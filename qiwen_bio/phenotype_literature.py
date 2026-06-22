"""Phenotype-oriented literature retrieval and conservative claim-support
classification for stage 3F.

This module connects cellular-process evidence that already has **direct**
support (UniProt GO annotation or direct KEGG membership) to phenotype-oriented
PubMed literature. Each retrieved article is classified at the
**abstract-metadata boundary only** into one of three support levels, and a
phenotype hypothesis is emitted **only** when a process link has at least one
article classified as ``supports``.

Biological-correctness guardrails (enforced here, not at the caller):

* Search retrieval is never "support" by itself.
* Classification uses title + abstract text only. No full-text appraisal and
  no causal inference is performed.
* Enrichment-only processes (no UniProt GO or direct KEGG support) are never
  hypothesis-eligible; they may still appear in ``literature_context``.
* Every hypothesis carries an explicit uncertainty string. Process association
  does not establish activity, direction, mechanism, causality, or phenotype.
"""

from __future__ import annotations

import re
from typing import Literal

from pydantic import BaseModel, Field

from qiwen_bio.cellular_processes import CellularProcessEvidence, ProcessSupport
from qiwen_bio.pubmed import LiteratureEvidence, PubMedArticle, PubMedServiceError


SupportLevel = Literal["supports", "mentions", "no_abstract"]


# Pinned, conservative outcome verbs. A claim is only "supports" when the
# abstract uses one of these relational/outcome verbs alongside both the gene
# and the process label. This list is intentionally short and clinical: each
# verb implies a stated relationship rather than mere co-occurrence. Editing
# it is a behaviour change and should be reviewed.
OUTCOME_VERBS: tuple[str, ...] = (
    "associated with",
    "correlates",
    "correlated",
    "correlation",
    "linked to",
    "linked with",
    "involved in",
    "involved",
    "promotes",
    "promote",
    "inhibits",
    "inhibit",
    "inhibition of",
    "activates",
    "activate",
    "activation of",
    "suppresses",
    "suppress",
    "suppression of",
    "induces",
    "induce",
    "induction of",
    "regulates",
    "regulate",
    "regulation of",
    "mediates",
    "mediate",
    "required for",
    "required",
    "disrupt",
    "disrupts",
    "disruption of",
    "loss of",
    "overexpression",
    "overexpressed",
    "drives",
    "drive",
    "modulates",
    "modulate",
)

# Process labels often contain generic stop words ("by", "of", "pathway",
# "signaling") that would produce false-positive matches. We require at least
# one *substantive* token from the label to appear in the abstract.
_PROCESS_STOPWORDS = frozenset(
    {
        "a",
        "an",
        "and",
        "by",
        "of",
        "or",
        "the",
        "to",
        "in",
        "for",
        "with",
        "via",
        "pathway",
        "pathways",
        "signaling",
        "signalling",
        "process",
        "response",
        "activity",
        "mediator",
        "class",
        "human",
        "homo",
        "sapiens",
    }
)


def _tokens(text: str) -> set[str]:
    return {match.group(0).lower() for match in re.finditer(r"[A-Za-z0-9_-]+", text)}


def _substantive_tokens(label: str) -> set[str]:
    tokens = _tokens(label)
    substantive = {token for token in tokens if token not in _PROCESS_STOPWORDS and len(token) > 2}
    # Fall back to the raw token set if the label is entirely stop words
    # (e.g. a one-word generic label) so we never claim a match on nothing.
    return substantive or tokens


def _gene_present(gene: str, article: PubMedArticle) -> bool:
    gene_token = gene.strip().upper()
    if not gene_token:
        return False
    haystack = f"{article.title} {article.abstract or ''}".upper()
    # Word-boundary match so "TP53" does not match "TP53BP1".
    return re.search(rf"(?<![A-Z0-9]){re.escape(gene_token)}(?![A-Z0-9])", haystack) is not None


def _process_present(process_label: str, article: PubMedArticle) -> bool:
    substantive = _substantive_tokens(process_label)
    if not substantive:
        return False
    text = f"{article.title} {article.abstract or ''}".lower()
    tokens = _tokens(text)
    # Require at least one substantive label token. Using "any" is deliberately
    # permissive for the *presence* test; the outcome-verb test is what gates
    # the stronger "supports" classification.
    return any(token in tokens for token in substantive)


def _matched_outcome_verbs(article: PubMedArticle) -> list[str]:
    text = f"{article.title} {article.abstract or ''}".lower()
    text = re.sub(r"\s+", " ", text)
    return [verb for verb in OUTCOME_VERBS if verb in text]


def classify_support(
    gene: str,
    process_label: str,
    article: PubMedArticle,
) -> tuple[SupportLevel, str, list[str]]:
    """Classify how strongly an article supports the gene->process link.

    Returns ``(support_level, basis, matched_phrases)``. ``basis`` is a short,
    auditable string describing which direct support fed the link and which
    phrase(s) matched, so a reader can trace why an article was labelled
    ``supports``. ``matched_phrases`` lists the pinned outcome verbs found in
    the title or abstract (empty unless the level is ``supports``).
    """
    if not article.abstract:
        return "no_abstract", "no abstract available at the abstract-metadata boundary", []

    if not _gene_present(gene, article):
        return "mentions", f"gene {gene} not co-mentioned in title or abstract", []

    if not _process_present(process_label, article):
        return "mentions", "process label tokens not co-mentioned in title or abstract", []

    verbs = _matched_outcome_verbs(article)
    if not verbs:
        return (
            "mentions",
            "gene and process co-mentioned; no outcome verb in the pinned conservative list",
            [],
        )

    # Deduplicate while preserving order, capping at three for readability.
    seen: set[str] = set()
    unique: list[str] = []
    for verb in verbs:
        if verb not in seen:
            seen.add(verb)
            unique.append(verb)
        if len(unique) == 3:
            break
    return (
        "supports",
        f"matched outcome verb(s): {', '.join(unique)}",
        unique,
    )


class PhenotypeClaimLink(BaseModel):
    process_id: str
    process_label: str
    gene: str
    support_level: SupportLevel
    pmid: str
    title: str
    evidence_basis: str
    matched_phrases: list[str] = Field(default_factory=list)
    uncertainty: str


class PhenotypeLiteratureEvidence(BaseModel):
    gene: str
    process_links: list[PhenotypeClaimLink]
    counts_by_level: dict[str, int]
    hypotheses: list[str]
    query_urls: list[str]
    disclaimer: str = (
        "Phenotype literature classification is performed at the abstract-metadata boundary "
        "only. It is not full-text appraisal, claim-level evidence grading, or causal inference."
    )
    boundary: str = (
        "A hypothesis is emitted only when a directly supported process has at least one "
        "article classified as 'supports'. Process association does not establish activity, "
        "effect direction, mechanism, causality, or a cellular or organismal phenotype."
    )


def _direct_support_sources(process_supports: list[ProcessSupport]) -> list[str]:
    return [
        support.source_name
        for support in process_supports
        if support.evidence_type in ("database_annotation", "database_membership")
    ]


def _is_directly_supported(process_supports: list[ProcessSupport]) -> bool:
    return any(
        support.evidence_type in ("database_annotation", "database_membership")
        for support in process_supports
    )


def _sort_articles_by_support(
    articles: list[PubMedArticle],
    gene: str,
    process_label: str,
) -> list[tuple[PubMedArticle, SupportLevel, str, list[str]]]:
    classified = [
        (article, *classify_support(gene, process_label, article))
        for article in articles
    ]
    priority = {"supports": 0, "mentions": 1, "no_abstract": 2}
    return sorted(classified, key=lambda item: (priority[item[1]], item[0].pmid))


def _build_hypothesis(
    gene: str,
    process_id: str,
    process_label: str,
    direct_sources: list[str],
    supports_links: list[PhenotypeClaimLink],
) -> str:
    pmids = ", ".join(f"PMID {link.pmid}" for link in supports_links)
    # Aggregate the matched outcome verbs across all supports articles,
    # deduplicated and capped, so the hypothesis cites concrete language.
    seen: set[str] = set()
    phrases: list[str] = []
    for link in supports_links:
        for phrase in link.matched_phrases:
            if phrase not in seen:
                seen.add(phrase)
                phrases.append(phrase)
    phrase_text = ", ".join(f"'{phrase}'" for phrase in phrases[:4]) or "outcome verbs"
    return (
        f"Because {', '.join(direct_sources)} links {gene} to {process_label} "
        f"({process_id}), and {pmids} use {phrase_text} in the abstract alongside "
        f"both {gene} and the process label, {gene} may be associated with "
        f"{process_label}. This is an abstract-level hypothesis, not a causal or "
        f"phenotype conclusion; it requires full-text appraisal and experimental "
        f"validation before use."
    )


def build_phenotype_literature(
    gene: str,
    cellular_processes: CellularProcessEvidence,
    pubmed_client,
    limit_per_process: int = 3,
) -> PhenotypeLiteratureEvidence:
    """Build phenotype-oriented literature evidence from cellular-process evidence.

    ``pubmed_client`` must expose ``search_with_abstracts(protein, context_terms,
    limit)``. Failures degrade gracefully: the process link is skipped and the
    boundary note is preserved, mirroring the optional-layer pattern used
    elsewhere in the pipeline.
    """
    normalized_gene = gene.strip().upper()
    eligible = [
        process
        for process in cellular_processes.processes
        if _is_directly_supported(process.supports)
    ]

    process_links: list[PhenotypeClaimLink] = []
    hypotheses: list[str] = []
    query_urls: list[str] = []

    for process in eligible:
        try:
            literature = pubmed_client.search_with_abstracts(
                protein=normalized_gene,
                context_terms=[process.label],
                limit=limit_per_process,
            )
        except PubMedServiceError:
            # Optional-layer degradation: skip this process, keep going.
            continue
        if literature.query:
            query_urls.append(literature.query)

        classified = _sort_articles_by_support(
            literature.articles, normalized_gene, process.label
        )[:limit_per_process]

        links_for_process: list[PhenotypeClaimLink] = []
        for article, level, basis, matched_phrases in classified:
            direct_sources = _direct_support_sources(process.supports)
            links_for_process.append(
                PhenotypeClaimLink(
                    process_id=process.canonical_id,
                    process_label=process.label,
                    gene=normalized_gene,
                    support_level=level,
                    pmid=article.pmid,
                    title=article.title,
                    evidence_basis=(
                        f"Direct support: {', '.join(direct_sources)}. {basis}"
                    ),
                    matched_phrases=list(matched_phrases),
                    uncertainty=(
                        "Abstract-level co-mention and outcome verb; not a causal, "
                        "directional, or phenotype conclusion."
                    ),
                )
            )
        process_links.extend(links_for_process)

        supports_links = [link for link in links_for_process if link.support_level == "supports"]
        if supports_links:
            direct_sources = _direct_support_sources(process.supports)
            hypotheses.append(
                _build_hypothesis(
                    normalized_gene,
                    process.canonical_id,
                    process.label,
                    direct_sources,
                    supports_links,
                )
            )

    if not eligible:
        boundary = (
            "No eligible directly supported processes were available for phenotype "
            "literature retrieval; enrichment-only processes are not hypothesis-eligible."
        )
    else:
        boundary = PhenotypeLiteratureEvidence.model_fields["boundary"].default

    counts: dict[str, int] = {}
    for link in process_links:
        counts[link.support_level] = counts.get(link.support_level, 0) + 1

    return PhenotypeLiteratureEvidence(
        gene=normalized_gene,
        process_links=process_links,
        counts_by_level=dict(sorted(counts.items())),
        hypotheses=hypotheses,
        query_urls=query_urls,
        boundary=boundary,
    )

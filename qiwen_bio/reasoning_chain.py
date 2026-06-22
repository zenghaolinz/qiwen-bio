"""Cross-scale reasoning chain (stage 3G).

This module assembles the already-retrieved evidence layers into a single
structured chain that connects a mutation (or, in its absence, the protein
itself) to structure, function, pathway, and phenotype evidence. It is a
**synthesis layer**: it performs no new network requests and introduces no new
biological claims of its own.

Biological-correctness guardrails (enforced here):

* The chain is a **pure evidence chain**. Each step states deterministic,
  source-bound facts first. A hypothesis sentence ("可能影响/可能相关/提示") is
  emitted only when explicit upstream evidence supports it, and every
  hypothesis carries an uncertainty clause. No step asserts causality,
  direction, mechanism, activity, or a phenotype outcome.
* For a mutation-impact chain, the wild-type residue is validated against the
  UniProt sequence. A mismatch suppresses downstream "可能影响" hypotheses:
  facts remain, but no functional-impact inference is drawn from a mutation
  that does not match the annotated sequence.
* Per ADR-0012/0013: process association and abstract-level literature do not
  establish activity, direction, mechanism, causality, or phenotype.
"""

from __future__ import annotations

import re
from typing import Literal

from pydantic import BaseModel

from qiwen_bio.alphafold import AlphaFoldAnalysis
from qiwen_bio.cellular_processes import CellularProcessEvidence
from qiwen_bio.interpro import DomainAnnotation
from qiwen_bio.kegg import KeggPathwayAnnotation
from qiwen_bio.models import AnalysisResponse
from qiwen_bio.phenotype_literature import PhenotypeLiteratureEvidence
from qiwen_bio.stringdb import EvidenceGraph
from qiwen_bio.uniprot import UniProtAnnotation


MUTATION_PATTERN = re.compile(r"^([A-Z])(\d+)([A-Z])$")

ChainType = Literal["mutation_impact", "protein_function"]
StepConfidence = Literal["high", "medium", "low", "insufficient"]


class ChainStep(BaseModel):
    step_id: str
    title: str
    evidence_facts: list[str]
    hypothesis: str | None = None
    evidence_sources: list[str]
    confidence: StepConfidence
    uncertainty: str
    available: bool


class ReasoningChain(BaseModel):
    chain_type: ChainType
    gene: str
    mutation: str | None
    steps: list[ChainStep]
    missing_layers: list[str]
    summary: str
    boundary: str = (
        "This chain organizes retrieved and calculated evidence. It does not establish "
        "activity, effect direction, mechanism, causality, pathogenicity, or a cellular "
        "or organismal phenotype. Hypothesis sentences are explicitly labelled and require "
        "full-text appraisal and experimental validation before use."
    )


def _parse_mutation(mutation: str | None) -> tuple[str, int, str] | None:
    """Parse a mutation string into (wild_type, position, mutant).

    Returns None if the mutation is missing or does not match the strict
    ``X123Y`` format. This is deliberately stricter than the legacy
    ``mutation[1:-1].isdigit()`` slice in synthesis.py and is the canonical
    parse for the reasoning chain.
    """
    if not mutation:
        return None
    match = MUTATION_PATTERN.fullmatch(mutation.strip().upper())
    if not match:
        return None
    wild_type, position_text, mutant = match.groups()
    return wild_type, int(position_text), mutant


def _wild_type_matches_sequence(
    sequence: str, wild_type: str, position: int
) -> bool:
    if position < 1 or position > len(sequence):
        return False
    return sequence[position - 1] == wild_type


def _confidence_from_plddt(plddt: float | None) -> StepConfidence:
    if plddt is None:
        return "insufficient"
    if plddt >= 90:
        return "high"
    if plddt >= 70:
        return "medium"
    if plddt >= 50:
        return "low"
    return "low"


def _direct_processes(cellular_processes: CellularProcessEvidence) -> list:
    """Processes with at least one non-enrichment support (direct annotation/membership)."""
    return [
        process
        for process in cellular_processes.processes
        if any(
            support.evidence_type in ("database_annotation", "database_membership")
            for support in process.supports
        )
    ]


def build_mutation_step(
    annotation: UniProtAnnotation,
    structure: AlphaFoldAnalysis | None,
    domains: DomainAnnotation | None,
    mutation: str | None,
) -> ChainStep:
    """Build the mutation step. Records a mismatch fact when the wild-type
    residue does not match the UniProt sequence, which downstream steps use to
    suppress functional-impact hypotheses."""
    parsed = _parse_mutation(mutation)
    if parsed is None:
        return ChainStep(
            step_id="mutation",
            title="Mutation input",
            evidence_facts=[f"Invalid or malformed mutation input: {mutation!r}."],
            evidence_sources=[],
            confidence="insufficient",
            uncertainty="No valid mutation could be parsed; downstream impact hypotheses are suppressed.",
            available=False,
        )

    wild_type, position, mutant = parsed
    facts: list[str] = [f"User-supplied mutation: {wild_type}{position}{mutant}."]
    sources: list[str] = ["User input"]

    matches = _wild_type_matches_sequence(annotation.sequence, wild_type, position)
    if matches:
        facts.append(
            f"Wild-type {wild_type} at position {position} matches the UniProt sequence."
        )
    else:
        actual = (
            annotation.sequence[position - 1]
            if 1 <= position <= len(annotation.sequence)
            else "absent"
        )
        facts.append(
            f"Wild-type mismatch: UniProt sequence has {actual} at position {position}, "
            f"not {wild_type}. The mutation may use different numbering or refer to "
            f"a different isoform; impact hypotheses are suppressed."
        )
    sources.append(annotation.source_url)

    domain_overlap_names: list[str] = []
    if domains and domains.mutation_overlaps:
        domain_overlap_names = [entry.accession for entry in domains.mutation_overlaps]
        facts.append(
            f"Mutation position overlaps {len(domains.mutation_overlaps)} annotated "
            f"domain entry/entries: {', '.join(domain_overlap_names)}."
        )
        sources.extend(entry.source_url for entry in domains.mutation_overlaps)
    elif domains:
        facts.append("Mutation position does not overlap any annotated InterPro/Pfam domain entry.")

    if structure and structure.mutation_site:
        site = structure.mutation_site
        facts.append(
            f"AlphaFold mutation site pLDDT: {site.plddt:.2f} ({site.confidence}); "
            f"{len(structure.mutation_neighborhood)} CA neighbours within 8 A."
        )
        sources.append(structure.structure_url)

    has_overlap = bool(domain_overlap_names)
    hypothesis: str | None = None
    if matches and has_overlap and structure and structure.mutation_site:
        hypothesis = (
            f"假设：因为 {wild_type}{position}{mutant} 位于结构域 "
            f"{', '.join(domain_overlap_names)} 内，该突变可能影响该结构域的局部结构或功能。"
            f"这是基于坐标重叠和模型置信度的假设，不是功能影响或致病性结论；"
            f"需实验验证。"
        )

    confidence: StepConfidence = "low"
    if not matches:
        confidence = "insufficient"
    elif has_overlap and structure and structure.mutation_site:
        confidence = _confidence_from_plddt(structure.mutation_site.plddt)

    return ChainStep(
        step_id="mutation",
        title="Mutation and structural context",
        evidence_facts=facts,
        hypothesis=hypothesis,
        evidence_sources=sources,
        confidence=confidence,
        uncertainty=(
            "Positional overlap is a coordinate containment observation, not a "
            "functional-effect or pathogenicity prediction."
        ),
        available=True,
    )


def build_structure_step(
    structure: AlphaFoldAnalysis | None, mutation: str | None
) -> ChainStep:
    if structure is None:
        return ChainStep(
            step_id="structure",
            title="Predicted structure",
            evidence_facts=["No AlphaFold structure is available for this entry."],
            evidence_sources=[],
            confidence="insufficient",
            uncertainty="Structure layer unavailable; geometric context cannot be assessed.",
            available=False,
        )
    facts = [
        f"AlphaFold model: {structure.residue_count} residues, mean pLDDT {structure.mean_plddt:.2f}.",
        f"Non-local CA contacts (8 A threshold): {structure.contact_map.total_contacts}.",
    ]
    if mutation and structure.mutation_site:
        site = structure.mutation_site
        facts.append(
            f"Mutation site {site.wild_type}{site.position}{site.mutant}: pLDDT {site.plddt:.2f} "
            f"({site.confidence}), {len(structure.mutation_neighborhood)} neighbours within 8 A."
        )
    return ChainStep(
        step_id="structure",
        title="Predicted structure",
        evidence_facts=facts,
        hypothesis=None,
        evidence_sources=[structure.structure_url],
        confidence=_confidence_from_plddt(structure.mean_plddt),
        uncertainty=(
            "pLDDT and CA proximity describe model confidence and geometry, not "
            "pathogenicity, stability, or functional effect."
        ),
        available=True,
    )


def build_function_step(
    annotation: UniProtAnnotation,
    domains: DomainAnnotation | None,
    mutation_step: ChainStep | None,
) -> ChainStep:
    facts: list[str] = []
    sources: list[str] = [annotation.source_url]
    if annotation.functions:
        facts.append(f"UniProt FUNCTION: {annotation.functions[0]}")
    else:
        facts.append("No UniProt FUNCTION annotation available.")
    bp_terms = [
        term for term in annotation.go_terms if term.aspect == "biological_process"
    ]
    facts.append(f"UniProt GO biological_process terms: {len(bp_terms)}.")
    if domains:
        names = ", ".join(entry.name for entry in domains.entries[:6])
        facts.append(f"InterPro/Pfam domain entries: {domains.entry_count} ({names}).")
        sources.extend(entry.source_url for entry in domains.entries[:6])
    else:
        facts.append("No InterPro/Pfam domain annotation available.")

    # A functional-impact hypothesis requires a matching mutation with a domain overlap.
    hypothesis: str | None = None
    if (
        mutation_step
        and mutation_step.available
        and mutation_step.hypothesis is not None
        and domains
        and domains.mutation_overlaps
    ):
        overlap_names = ", ".join(entry.accession for entry in domains.mutation_overlaps)
        hypothesis = (
            f"假设：结合结构域 {overlap_names} 的坐标重叠，该突变可能影响蛋白功能。"
            f"这是基于注释重叠的假设，不是功能丧失或获得结论；需实验验证。"
        )

    confidence: StepConfidence = "medium" if annotation.functions else "low"
    return ChainStep(
        step_id="function",
        title="Function and domain annotation",
        evidence_facts=facts,
        hypothesis=hypothesis,
        evidence_sources=sources,
        confidence=confidence,
        uncertainty=(
            "Database annotations describe known function; they do not predict how a "
            "specific mutation alters that function."
        ),
        available=True,
    )


def build_pathway_step(
    cellular_processes: CellularProcessEvidence,
    kegg: KeggPathwayAnnotation | None,
    graph: EvidenceGraph | None,
    function_step: ChainStep,
) -> ChainStep:
    facts: list[str] = []
    sources: list[str] = []
    direct = _direct_processes(cellular_processes)
    if direct:
        ids = ", ".join(process.canonical_id for process in direct[:8])
        facts.append(f"Directly supported processes/pathways: {len(direct)} ({ids}).")
        sources.extend(
            support.source_url
            for process in direct
            for support in process.supports
            if support.evidence_type in ("database_annotation", "database_membership")
        )
    else:
        facts.append("No directly supported process/pathway evidence (UniProt GO or direct KEGG).")

    if kegg:
        facts.append(
            f"Direct KEGG memberships: {kegg.pathway_count} pathway(s) for gene "
            f"{', '.join(kegg.gene_ids)}."
        )
        sources.append("https://www.kegg.jp/")

    if graph:
        interactions = sum(1 for edge in graph.edges if edge.type == "interacts_with")
        terms = [node for node in graph.nodes if node.type != "protein"]
        facts.append(f"STRING interaction edges: {interactions}; enrichment terms: {len(terms)}.")
        sources.append("https://string-db.org/")

    # A pathway-association hypothesis requires a functional-impact hypothesis upstream.
    hypothesis: str | None = None
    if function_step.hypothesis is not None and direct:
        hypothesis = (
            f"假设：如果该蛋白功能受影响，其参与的通路（如 {direct[0].label}）可能相关。"
            f"这是基于功能假设和数据库成员关系的二级假设，不是通路激活、方向性或表型因果结论；"
            f"需实验验证。"
        )

    confidence: StepConfidence = "medium" if direct else "low"
    return ChainStep(
        step_id="pathway",
        title="Pathway and interaction context",
        evidence_facts=facts,
        hypothesis=hypothesis,
        evidence_sources=sources,
        confidence=confidence,
        uncertainty=(
            "Pathway membership and interaction scores are database associations, not "
            "evidence of pathway activation, directionality, or phenotype causality."
        ),
        available=True,
    )


def build_phenotype_step(
    phenotype_literature: PhenotypeLiteratureEvidence | None,
    pathway_step: ChainStep,
) -> ChainStep:
    facts: list[str] = []
    sources: list[str] = []
    if phenotype_literature:
        counts = ", ".join(
            f"{level}={count}"
            for level, count in phenotype_literature.counts_by_level.items()
        ) or "none"
        facts.append(
            f"Phenotype literature claim-support: {counts}; "
            f"{len(phenotype_literature.process_links)} process-link classification(s)."
        )
        pmids = [link.pmid for link in phenotype_literature.process_links[:5]]
        if pmids:
            facts.append(f"Classified PMIDs: {', '.join(pmids)}.")
            sources.extend(
                f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" for pmid in pmids
            )
    else:
        facts.append("No phenotype literature layer available.")

    hypothesis: str | None = None
    if phenotype_literature and phenotype_literature.hypotheses:
        hypothesis = phenotype_literature.hypotheses[0]
    elif phenotype_literature:
        facts.append(
            "No phenotype hypothesis was emitted: no directly supported process had an "
            "article classified as 'supports' at the abstract-metadata boundary."
        )

    confidence: StepConfidence = "low"
    if phenotype_literature and phenotype_literature.counts_by_level.get("supports", 0) > 0:
        confidence = "medium"

    return ChainStep(
        step_id="phenotype",
        title="Phenotype literature evidence",
        evidence_facts=facts,
        hypothesis=hypothesis,
        evidence_sources=sources,
        confidence=confidence,
        uncertainty=(
            "Abstract-metadata classification is not full-text appraisal, claim-level "
            "grading, or causal inference."
        ),
        available=phenotype_literature is not None,
    )


def build_reasoning_chain(
    annotation: UniProtAnnotation,
    analysis: AnalysisResponse,
    structure: AlphaFoldAnalysis | None,
    domains: DomainAnnotation | None,
    kegg: KeggPathwayAnnotation | None,
    graph: EvidenceGraph | None,
    cellular_processes: CellularProcessEvidence,
    phenotype_literature: PhenotypeLiteratureEvidence | None,
    mutation: str | None,
) -> ReasoningChain:
    """Assemble the cross-scale reasoning chain from already-retrieved layers."""
    gene = annotation.gene_names[0] if annotation.gene_names else annotation.accession
    # chain_type reflects user intent: any supplied mutation string (even a
    # malformed one) triggers the mutation-impact chain type, so the chain
    # surfaces the invalid input rather than silently switching types.
    chain_type: ChainType = "mutation_impact" if mutation else "protein_function"

    steps: list[ChainStep] = []
    missing_layers: list[str] = []

    mutation_step: ChainStep | None = None
    if chain_type == "mutation_impact":
        mutation_step = build_mutation_step(annotation, structure, domains, mutation)
        steps.append(mutation_step)

    structure_step = build_structure_step(structure, mutation)
    if not structure_step.available:
        missing_layers.append("structure")
    steps.append(structure_step)

    function_step = build_function_step(annotation, domains, mutation_step)
    steps.append(function_step)

    pathway_step = build_pathway_step(cellular_processes, kegg, graph, function_step)
    steps.append(pathway_step)
    if not cellular_processes.processes:
        missing_layers.append("cellular_processes")

    phenotype_step = build_phenotype_step(phenotype_literature, pathway_step)
    if not phenotype_step.available:
        missing_layers.append("phenotype_literature")
    steps.append(phenotype_step)

    if kegg is None:
        missing_layers.append("kegg")
    if graph is None:
        missing_layers.append("string_graph")

    if chain_type == "mutation_impact":
        summary = (
            f"Mutation-impact reasoning chain for {gene} {mutation}: connects the mutation "
            f"to structural context, function, pathway, and phenotype literature evidence."
        )
    else:
        summary = (
            f"Protein-function reasoning chain for {gene}: connects domain/function "
            f"annotation, pathway, and phenotype literature evidence (no mutation supplied)."
        )

    return ReasoningChain(
        chain_type=chain_type,
        gene=gene,
        mutation=mutation if chain_type == "mutation_impact" else None,
        steps=steps,
        missing_layers=missing_layers,
        summary=summary,
    )

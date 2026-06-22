from qiwen_bio.cellular_processes import (
    CellularProcess,
    CellularProcessEvidence,
    ProcessSupport,
)
from qiwen_bio.phenotype_literature import (
    build_phenotype_literature,
    classify_support,
)
from qiwen_bio.pubmed import LiteratureEvidence, PubMedArticle, PubMedServiceError


def _process(
    canonical_id: str,
    label: str,
    supports: list[ProcessSupport],
) -> CellularProcess:
    return CellularProcess(
        canonical_id=canonical_id,
        label=label,
        process_type="biological_process",
        supports=supports,
    )


def _support(source_name: str, evidence_type: str) -> ProcessSupport:
    return ProcessSupport(
        source_name=source_name,
        evidence_type=evidence_type,
        source_url="https://example.test/source",
        detail=f"{source_name} support",
    )


def _evidence(processes: list[CellularProcess]) -> CellularProcessEvidence:
    return CellularProcessEvidence(
        protein_accession="P04637",
        processes=processes,
        counts_by_source={},
        literature_context=[process.label for process in processes],
    )


def _article(
    pmid: str,
    title: str,
    abstract: str | None,
) -> PubMedArticle:
    return PubMedArticle(
        pmid=pmid,
        title=title,
        authors=["Smith A"],
        journal="Evidence Journal",
        published="2025 Jan",
        doi="10.1000/example",
        url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
        abstract=abstract,
        abstract_sections=[(None, abstract)] if abstract else [],
    )


def test_classify_support_requires_gene_process_and_outcome_verb() -> None:
    article = _article(
        "1",
        "TP53 and p53 signaling pathway.",
        "TP53 regulates the p53 signaling pathway during cell cycle arrest.",
    )
    level, basis, phrases = classify_support("TP53", "p53 signaling pathway", article)
    assert level == "supports"
    assert "regulates" in basis
    assert "regulates" in phrases


def test_classify_support_returns_mentions_without_outcome_verb() -> None:
    article = _article(
        "2",
        "TP53 and p53 signaling pathway.",
        "TP53 and the p53 signaling pathway are discussed.",
    )
    level, _, phrases = classify_support("TP53", "p53 signaling pathway", article)
    assert level == "mentions"
    assert phrases == []


def test_classify_support_returns_no_abstract_when_missing() -> None:
    article = _article("3", "TP53 and p53 signaling pathway.", None)
    level, _, phrases = classify_support("TP53", "p53 signaling pathway", article)
    assert level == "no_abstract"
    assert phrases == []


def test_classify_support_returns_mentions_when_gene_absent() -> None:
    article = _article(
        "4",
        "MDM2 study.",
        "MDM2 regulates the p53 signaling pathway.",
    )
    level, _, phrases = classify_support("TP53", "p53 signaling pathway", article)
    assert level == "mentions"
    assert phrases == []


class StubPhenotypePubMedClient:
    """Returns a fixed LiteratureEvidence per (protein, context_terms) call.

    The first context_term is the process label, so the stub keys responses by it.
    """

    def __init__(self, responses: dict[str, LiteratureEvidence]) -> None:
        self.responses = responses
        self.calls: list[tuple[str, list[str], int]] = []

    def search_with_abstracts(
        self, protein: str, context_terms: list[str], limit: int
    ) -> LiteratureEvidence:
        self.calls.append((protein, list(context_terms), limit))
        key = context_terms[0] if context_terms else ""
        if key not in self.responses:
            return LiteratureEvidence(
                protein=protein,
                context_terms=context_terms,
                query='"PROTEIN"[Title/Abstract]',
                articles=[],
            )
        return self.responses[key]


def _direct_supported_process() -> CellularProcess:
    return _process(
        "GO:0072331",
        "signal transduction by p53 class mediator",
        [
            _support("UniProt GO", "database_annotation"),
            _support("STRING enrichment", "enrichment_statistic"),
        ],
    )


def test_builder_gates_hypothesis_on_direct_support_and_supports_article() -> None:
    process = _direct_supported_process()
    evidence = _evidence([process])
    article = _article(
        "111",
        "TP53 and p53 signaling.",
        "TP53 activates signal transduction by p53 class mediator in apoptosis.",
    )
    pubmed = StubPhenotypePubMedClient(
        {
            "signal transduction by p53 class mediator": LiteratureEvidence(
                protein="TP53",
                context_terms=["signal transduction by p53 class mediator"],
                query='"TP53"[Title/Abstract]',
                articles=[article],
            )
        }
    )

    result = build_phenotype_literature("TP53", evidence, pubmed, limit_per_process=3)

    assert result.gene == "TP53"
    assert len(result.process_links) == 1
    link = result.process_links[0]
    assert link.process_id == "GO:0072331"
    assert link.support_level == "supports"
    assert link.pmid == "111"
    assert "UniProt GO" in link.evidence_basis
    assert link.uncertainty
    assert result.counts_by_level == {"supports": 1}
    assert len(result.hypotheses) == 1
    assert "TP53" in result.hypotheses[0]
    assert "signal transduction by p53 class mediator" in result.hypotheses[0]
    assert "abstract-level hypothesis" in result.hypotheses[0]
    assert pubmed.calls[0][0] == "TP53"
    assert pubmed.calls[0][1] == ["signal transduction by p53 class mediator"]


def test_builder_emits_no_hypothesis_when_only_mentions_articles() -> None:
    process = _direct_supported_process()
    evidence = _evidence([process])
    article = _article(
        "222",
        "TP53 and p53 mediator.",
        "TP53 and signal transduction by p53 class mediator are described.",
    )
    pubmed = StubPhenotypePubMedClient(
        {
            "signal transduction by p53 class mediator": LiteratureEvidence(
                protein="TP53",
                context_terms=["signal transduction by p53 class mediator"],
                query='"TP53"[Title/Abstract]',
                articles=[article],
            )
        }
    )

    result = build_phenotype_literature("TP53", evidence, pubmed)

    assert result.counts_by_level == {"mentions": 1}
    assert result.hypotheses == []
    assert result.process_links[0].support_level == "mentions"
    assert result.boundary


def test_builder_skips_enrichment_only_processes() -> None:
    enrichment_only = _process(
        "GO:9999999",
        "neighbor-only process",
        [_support("STRING enrichment", "enrichment_statistic")],
    )
    evidence = _evidence([enrichment_only])
    pubmed = StubPhenotypePubMedClient({})

    result = build_phenotype_literature("TP53", evidence, pubmed)

    assert result.process_links == []
    assert result.hypotheses == []
    assert pubmed.calls == []
    assert "no eligible" in result.boundary.lower()


def test_builder_degrades_gracefully_on_pubmed_failure() -> None:
    process = _direct_supported_process()
    evidence = _evidence([process])

    class FailingPubMedClient:
        def search_with_abstracts(
            self, protein: str, context_terms: list[str], limit: int
        ) -> LiteratureEvidence:
            raise PubMedServiceError("efetch unavailable")

    result = build_phenotype_literature("TP53", evidence, FailingPubMedClient())

    assert result.process_links == []
    assert result.hypotheses == []
    assert result.boundary


def test_builder_caps_links_per_process_and_prioritizes_supports() -> None:
    process = _direct_supported_process()
    evidence = _evidence([process])
    supports_article = _article(
        "301",
        "TP53 activates p53 mediator.",
        "TP53 activates signal transduction by p53 class mediator.",
    )
    mention_article = _article(
        "302",
        "TP53 and p53 mediator.",
        "TP53 and signal transduction by p53 class mediator are mentioned.",
    )
    extra_supports = _article(
        "303",
        "TP53 promotes p53 mediator.",
        "TP53 promotes signal transduction by p53 class mediator.",
    )
    pubmed = StubPhenotypePubMedClient(
        {
            "signal transduction by p53 class mediator": LiteratureEvidence(
                protein="TP53",
                context_terms=["signal transduction by p53 class mediator"],
                query='"TP53"[Title/Abstract]',
                articles=[mention_article, supports_article, extra_supports],
            )
        }
    )

    result = build_phenotype_literature("TP53", evidence, pubmed, limit_per_process=2)

    levels = [link.support_level for link in result.process_links]
    assert levels.count("supports") == 2
    assert "mentions" not in levels
    pmids = [link.pmid for link in result.process_links]
    assert "301" in pmids and "303" in pmids

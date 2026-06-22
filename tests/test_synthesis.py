from qiwen_bio.alphafold import AlphaFoldServiceError, parse_alphafold_pdb
from qiwen_bio.pipeline import AnalysisPipeline
from qiwen_bio.pubmed import LiteratureEvidence, PubMedArticle, PubMedServiceError
from qiwen_bio.stringdb import StringServiceError, build_evidence_graph
from qiwen_bio.interpro import (
    DomainAnnotation,
    DomainEntry,
    DomainLocation,
    InterProServiceError,
)
from qiwen_bio.kegg import (
    KeggPathway,
    KeggPathwayAnnotation,
    KeggServiceError,
)
from qiwen_bio.synthesis import build_comprehensive_analysis
from qiwen_bio.uniprot import parse_uniprot_record
from tests.test_alphafold import PDB_TEXT
from tests.test_string_graph import ENRICHMENT_RECORDS, NETWORK_RECORDS
from tests.test_uniprot import UNIPROT_RECORD


class StubUniProtClient:
    def resolve(self, identifier: str, organism_id: int = 9606):
        return parse_uniprot_record(UNIPROT_RECORD)


class StubAlphaFoldClient:
    def analyze(self, accession: str, mutation: str | None = None):
        return parse_alphafold_pdb(accession, PDB_TEXT, "https://example.test/model.pdb", mutation)


class StubStringClient:
    def build_graph(self, identifier: str, species: int, limit: int, required_score: int):
        return build_evidence_graph(identifier, species, NETWORK_RECORDS, ENRICHMENT_RECORDS)


class StubPubMedClient:
    def __init__(self) -> None:
        self.context_terms: list[str] = []
        self.abstract_calls: list[tuple[str, list[str], int]] = []

    def search(self, protein: str, context_terms: list[str], limit: int):
        self.context_terms = context_terms
        return LiteratureEvidence(
            protein=protein,
            context_terms=context_terms,
            query='"TP53"[Title/Abstract]',
            articles=[
                PubMedArticle(
                    pmid="12345",
                    title="TP53 evidence study.",
                    authors=["Smith A"],
                    journal="Evidence Journal",
                    published="2025 Jan",
                    doi="10.1000/example",
                    url="https://pubmed.ncbi.nlm.nih.gov/12345/",
                )
            ],
        )

    def search_with_abstracts(self, protein: str, context_terms: list[str], limit: int):
        self.abstract_calls.append((protein, list(context_terms), limit))
        # Return a supports-level article for any process label so the
        # hypothesis gate fires and phenotype_hypotheses is populated.
        label = context_terms[0] if context_terms else "process"
        return LiteratureEvidence(
            protein=protein,
            context_terms=context_terms,
            query=f'"TP53"[Title/Abstract] AND "{label}"[Title/Abstract]',
            articles=[
                PubMedArticle(
                    pmid="99999",
                    title=f"TP53 regulates {label}.",
                    authors=["Smith A"],
                    journal="Evidence Journal",
                    published="2025 Jan",
                    doi="10.1000/example",
                    url="https://pubmed.ncbi.nlm.nih.gov/99999/",
                    abstract=f"TP53 regulates {label} in tumor suppression.",
                    abstract_sections=[(None, f"TP53 regulates {label} in tumor suppression.")],
                )
            ],
        )


class StubInterProClient:
    def fetch(self, accession: str, mutation_position: int | None = None):
        assert accession == "P04637"
        entry = DomainEntry(
            accession="PF00870",
            name="P53 DNA-binding domain",
            source_database="pfam",
            entry_type="domain",
            integrated_accession="IPR011615",
            source_url="https://www.ebi.ac.uk/interpro/entry/pfam/PF00870/",
            locations=[DomainLocation(start=1, end=8, status="CONTINUOUS")],
            go_terms=[],
            overlaps_mutation=mutation_position == 2,
        )
        return DomainAnnotation(
            protein_accession=accession,
            protein_length=10,
            mutation_position=mutation_position,
            entries=[entry],
            mutation_overlaps=[entry] if entry.overlaps_mutation else [],
            entry_count=1,
            location_count=1,
            source_urls=["https://www.ebi.ac.uk/interpro/api/"],
        )


class StubKeggClient:
    def fetch(self, accession: str, limit: int = 20):
        assert accession == "P04637"
        return KeggPathwayAnnotation(
            protein_accession=accession,
            gene_ids=["hsa:7157"],
            pathways=[
                KeggPathway(
                    pathway_id="hsa04115",
                    name="p53 signaling pathway - Homo sapiens (human)",
                    description="p53 stress response.",
                    classes=["Cellular Processes", "Cell growth and death"],
                    source_url="https://www.kegg.jp/entry/hsa04115",
                )
            ],
            pathway_count=1,
            linked_pathway_count=1,
            truncated=False,
            query_urls=["https://rest.kegg.jp/link/pathway/hsa:7157"],
        )


def test_comprehensive_analysis_combines_all_layers_and_scores_coverage() -> None:
    pubmed = StubPubMedClient()
    result = build_comprehensive_analysis(
        identifier="TP53",
        organism_id=9606,
        mutation="R2H",
        pipeline=AnalysisPipeline(),
        uniprot_client=StubUniProtClient(),
        alphafold_client=StubAlphaFoldClient(),
        string_client=StubStringClient(),
        pubmed_client=pubmed,
        interpro_client=StubInterProClient(),
        kegg_client=StubKeggClient(),
    )

    assert result.coverage.score == 100
    assert result.coverage.label == "comprehensive"
    assert all(component.available for component in result.coverage.components)
    assert result.structure.mutation_site.position == 2
    assert result.graph.seed == "TP53"
    assert result.literature.articles[0].pmid == "12345"
    assert result.domains.mutation_overlaps[0].accession == "PF00870"
    assert result.kegg.pathways[0].pathway_id == "hsa04115"
    assert [item.canonical_id for item in result.cellular_processes.processes] == [
        "GO:0072331",
        "hsa04115",
    ]
    normalized_pathway = result.cellular_processes.processes[1]
    assert [support.source_name for support in normalized_pathway.supports] == [
        "KEGG",
        "STRING enrichment",
    ]
    # Phenotype literature is populated for directly supported processes and
    # the hypothesis gate fires because the stub returns supports-level articles.
    # Per ADR-0012 the cellular-process layer keeps its own list empty; the
    # hypotheses live on the separate PhenotypeLiteratureEvidence object.
    assert result.phenotype_literature is not None
    assert result.phenotype_literature.gene == "TP53"
    assert result.phenotype_literature.counts_by_level.get("supports", 0) >= 1
    assert result.cellular_processes.phenotype_hypotheses == []
    assert len(result.phenotype_literature.hypotheses) >= 1
    assert "abstract-level hypothesis" in result.phenotype_literature.hypotheses[0]
    # Cross-scale reasoning chain (stage 3G): mutation-impact chain with 5 steps.
    assert result.reasoning_chain is not None
    assert result.reasoning_chain.chain_type == "mutation_impact"
    assert result.reasoning_chain.mutation == "R2H"
    assert [step.step_id for step in result.reasoning_chain.steps] == [
        "mutation", "structure", "function", "pathway", "phenotype"
    ]
    assert "## Reasoning chain" in result.report_markdown
    assert pubmed.context_terms[0].startswith("p53 signaling pathway")
    assert any(component.name == "domains" for component in result.coverage.components)
    assert "## Structure evidence" in result.report_markdown
    assert "## Interaction and pathway evidence" in result.report_markdown
    assert "[PMID 12345]" in result.report_markdown
    assert "## Domain evidence" in result.report_markdown
    assert "overlaps mutation position 2" in result.report_markdown
    assert "## Direct KEGG pathway evidence" in result.report_markdown
    assert "Returned pathway records: 1 of 1 linked" in result.report_markdown
    assert "## Cellular-process evidence" in result.report_markdown
    # Regression: every normalized process must render (the append was
    # previously outside the loop, so only the last process appeared).
    assert "GO:0072331" in result.report_markdown
    assert "hsa04115" in result.report_markdown
    assert "## Phenotype literature evidence" in result.report_markdown
    assert "abstract-level hypothesis" in result.report_markdown
    pathway_component = next(
        component for component in result.coverage.components if component.name == "pathways"
    )
    assert "direct KEGG" in pathway_component.detail


class FailingAlphaFoldClient:
    def analyze(self, accession: str, mutation: str | None = None):
        raise AlphaFoldServiceError("structure unavailable")


class FailingStringClient:
    def build_graph(self, identifier: str, species: int, limit: int, required_score: int):
        raise StringServiceError("graph unavailable")


class FailingPubMedClient:
    def search(self, protein: str, context_terms: list[str], limit: int):
        raise PubMedServiceError("literature unavailable")

    def search_with_abstracts(self, protein: str, context_terms: list[str], limit: int):
        raise PubMedServiceError("literature unavailable")


class FailingInterProClient:
    def fetch(self, accession: str, mutation_position: int | None = None):
        raise InterProServiceError("domains unavailable")


class FailingKeggClient:
    def fetch(self, accession: str, limit: int = 20):
        raise KeggServiceError("pathways unavailable")


def test_optional_service_failures_return_partial_report_with_warnings() -> None:
    result = build_comprehensive_analysis(
        identifier="TP53",
        organism_id=9606,
        mutation=None,
        pipeline=AnalysisPipeline(),
        uniprot_client=StubUniProtClient(),
        alphafold_client=FailingAlphaFoldClient(),
        string_client=FailingStringClient(),
        pubmed_client=FailingPubMedClient(),
        interpro_client=FailingInterProClient(),
        kegg_client=FailingKeggClient(),
    )

    assert result.coverage.score == 25
    assert result.coverage.label == "limited"
    assert result.structure is None
    assert result.graph is None
    assert result.literature is None
    assert result.domains is None
    assert result.kegg is None
    assert result.cellular_processes.processes == []
    assert len(result.warnings) == 5
    assert "## Unavailable evidence layers" in result.report_markdown
    # No mutation supplied → protein_function chain with no mutation step.
    assert result.reasoning_chain is not None
    assert result.reasoning_chain.chain_type == "protein_function"
    assert "mutation" not in {step.step_id for step in result.reasoning_chain.steps}
    assert "structure" in result.reasoning_chain.missing_layers


def test_kegg_failure_uses_explicit_string_enrichment_fallback() -> None:
    result = build_comprehensive_analysis(
        identifier="TP53",
        organism_id=9606,
        mutation="R2H",
        pipeline=AnalysisPipeline(),
        uniprot_client=StubUniProtClient(),
        alphafold_client=StubAlphaFoldClient(),
        string_client=StubStringClient(),
        pubmed_client=StubPubMedClient(),
        interpro_client=StubInterProClient(),
        kegg_client=FailingKeggClient(),
    )

    pathway_component = next(
        component for component in result.coverage.components if component.name == "pathways"
    )
    assert pathway_component.available is True
    assert "STRING enrichment fallback" in pathway_component.detail
    assert any(warning.startswith("KEGG:") for warning in result.warnings)
    assert len(result.cellular_processes.processes) == 2


class MentionsOnlyPubMedClient(StubPubMedClient):
    """Returns articles that mention the gene and process but use no outcome
    verb, so the hypothesis gate stays closed."""

    def search_with_abstracts(self, protein: str, context_terms: list[str], limit: int):
        label = context_terms[0] if context_terms else "process"
        return LiteratureEvidence(
            protein=protein,
            context_terms=context_terms,
            query=f'"TP53"[Title/Abstract] AND "{label}"[Title/Abstract]',
            articles=[
                PubMedArticle(
                    pmid="55555",
                    title=f"TP53 and {label}.",
                    authors=["Smith A"],
                    journal="Evidence Journal",
                    published="2025 Jan",
                    doi=None,
                    url="https://pubmed.ncbi.nlm.nih.gov/55555/",
                    abstract=f"TP53 and {label} are described in this review.",
                    abstract_sections=[(None, f"TP53 and {label} are described in this review.")],
                )
            ],
        )


def test_phenotype_literature_gate_stays_closed_without_supports_articles() -> None:
    result = build_comprehensive_analysis(
        identifier="TP53",
        organism_id=9606,
        mutation="R2H",
        pipeline=AnalysisPipeline(),
        uniprot_client=StubUniProtClient(),
        alphafold_client=StubAlphaFoldClient(),
        string_client=StubStringClient(),
        pubmed_client=MentionsOnlyPubMedClient(),
        interpro_client=StubInterProClient(),
        kegg_client=StubKeggClient(),
    )

    assert result.phenotype_literature is not None
    assert result.phenotype_literature.counts_by_level.get("mentions", 0) >= 1
    assert result.phenotype_literature.counts_by_level.get("supports", 0) == 0
    assert result.phenotype_literature.hypotheses == []
    assert result.cellular_processes.phenotype_hypotheses == []
    assert "no phenotype hypothesis" in result.report_markdown.lower()


def test_phenotype_literature_degrades_gracefully_when_pubmed_fails() -> None:
    result = build_comprehensive_analysis(
        identifier="TP53",
        organism_id=9606,
        mutation="R2H",
        pipeline=AnalysisPipeline(),
        uniprot_client=StubUniProtClient(),
        alphafold_client=StubAlphaFoldClient(),
        string_client=StubStringClient(),
        pubmed_client=FailingPubMedClient(),
        interpro_client=StubInterProClient(),
        kegg_client=StubKeggClient(),
    )

    # Direct-supported processes exist, but PubMed failures are swallowed by
    # the optional-layer pattern: no phenotype links, no hypotheses, and the
    # main literature layer surfaces its own warning.
    assert result.phenotype_literature is not None
    assert result.phenotype_literature.process_links == []
    assert result.phenotype_literature.hypotheses == []
    assert result.cellular_processes.phenotype_hypotheses == []
    assert any(warning.startswith("PubMed:") for warning in result.warnings)

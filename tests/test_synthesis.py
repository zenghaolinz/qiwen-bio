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
    assert result.cellular_processes.phenotype_hypotheses == []
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
    assert "No phenotype hypotheses were generated automatically" in result.report_markdown
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

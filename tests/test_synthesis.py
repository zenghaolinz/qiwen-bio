from qiwen_bio.alphafold import AlphaFoldServiceError, parse_alphafold_pdb
from qiwen_bio.pipeline import AnalysisPipeline
from qiwen_bio.pubmed import LiteratureEvidence, PubMedArticle, PubMedServiceError
from qiwen_bio.stringdb import StringServiceError, build_evidence_graph
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
    def search(self, protein: str, context_terms: list[str], limit: int):
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


def test_comprehensive_analysis_combines_all_layers_and_scores_coverage() -> None:
    result = build_comprehensive_analysis(
        identifier="TP53",
        organism_id=9606,
        mutation="R2H",
        pipeline=AnalysisPipeline(),
        uniprot_client=StubUniProtClient(),
        alphafold_client=StubAlphaFoldClient(),
        string_client=StubStringClient(),
        pubmed_client=StubPubMedClient(),
    )

    assert result.coverage.score == 100
    assert result.coverage.label == "comprehensive"
    assert all(component.available for component in result.coverage.components)
    assert result.structure.mutation_site.position == 2
    assert result.graph.seed == "TP53"
    assert result.literature.articles[0].pmid == "12345"
    assert "## Structure evidence" in result.report_markdown
    assert "## Interaction and pathway evidence" in result.report_markdown
    assert "[PMID 12345]" in result.report_markdown


class FailingAlphaFoldClient:
    def analyze(self, accession: str, mutation: str | None = None):
        raise AlphaFoldServiceError("structure unavailable")


class FailingStringClient:
    def build_graph(self, identifier: str, species: int, limit: int, required_score: int):
        raise StringServiceError("graph unavailable")


class FailingPubMedClient:
    def search(self, protein: str, context_terms: list[str], limit: int):
        raise PubMedServiceError("literature unavailable")


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
    )

    assert result.coverage.score == 35
    assert result.coverage.label == "limited"
    assert result.structure is None
    assert result.graph is None
    assert result.literature is None
    assert len(result.warnings) == 3
    assert "## Unavailable evidence layers" in result.report_markdown

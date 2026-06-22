import httpx
from fastapi.testclient import TestClient

from qiwen_bio.api import (
    app,
    get_alphafold_client,
    get_embedding_service,
    get_interpro_client,
    get_kegg_client,
    get_pubmed_client,
    get_string_client,
)
from qiwen_bio.api import get_uniprot_client
from qiwen_bio.uniprot import parse_uniprot_record
from tests.test_uniprot import UNIPROT_RECORD
from tests.test_alphafold import PDB_TEXT
from qiwen_bio.alphafold import parse_alphafold_pdb
from qiwen_bio.stringdb import build_evidence_graph
from tests.test_string_graph import ENRICHMENT_RECORDS, NETWORK_RECORDS
from qiwen_bio.pubmed import PubMedClient
from tests.test_pubmed import SUMMARY_PAYLOAD
from tests.test_synthesis import (
    StubAlphaFoldClient as SynthesisAlphaFoldClient,
    StubPubMedClient as SynthesisPubMedClient,
    StubStringClient as SynthesisStringClient,
    StubUniProtClient as SynthesisUniProtClient,
    StubInterProClient as SynthesisInterProClient,
    StubKeggClient as SynthesisKeggClient,
)
from qiwen_bio.embedding import EmbeddingCache, EmbeddingService, ModelLoadError
from tests.test_embedding import FakeProvider
from qiwen_bio.interpro import (
    DomainAnnotation,
    DomainEntry,
    DomainLocation,
    InterProNotFoundError,
    InterProServiceError,
)
from qiwen_bio.kegg import (
    KeggNotFoundError,
    KeggPathway,
    KeggPathwayAnnotation,
    KeggServiceError,
)


client = TestClient(app)


def test_health() -> None:
    response = client.get("/api/health")
    assert response.status_code == 200
    assert response.json() == {"status": "ok", "version": "0.1.0"}


def test_analyze_endpoint() -> None:
    response = client.post(
        "/api/v1/analyze",
        json={"name": "Peptide", "sequence": "KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK"},
    )
    assert response.status_code == 200
    payload = response.json()
    assert payload["name"] == "Peptide"
    assert payload["features"]["length"] == 37
    assert payload["prediction"]["task"] == "antimicrobial_peptide_demo"
    assert len(payload["evidence_chain"]) == 3


def test_analyze_endpoint_rejects_invalid_sequence() -> None:
    response = client.post("/api/v1/analyze", json={"sequence": "DNA?"})
    assert response.status_code == 422


def test_analyze_uniprot_uses_resolved_sequence_and_adds_provenance() -> None:
    class StubUniProtClient:
        def resolve(self, identifier: str, organism_id: int = 9606):
            assert identifier == "TP53"
            assert organism_id == 9606
            return parse_uniprot_record(UNIPROT_RECORD)

    app.dependency_overrides[get_uniprot_client] = lambda: StubUniProtClient()
    try:
        response = client.post(
            "/api/v1/analyze/uniprot",
            json={"identifier": "TP53", "organism_id": 9606, "mutation": "R5H"},
        )
    finally:
        app.dependency_overrides.clear()

    assert response.status_code == 200
    payload = response.json()
    assert payload["annotation"]["accession"] == "P04637"
    assert payload["analysis"]["features"]["length"] == 10
    assert payload["analysis"]["evidence_chain"][0]["source"].endswith("/P04637")
    assert payload["analysis"]["evidence_chain"][0]["evidence_type"] == "database_record"


def test_alphafold_endpoint_returns_plddt_and_mutation_context() -> None:
    class StubAlphaFoldClient:
        def analyze(self, accession: str, mutation: str | None = None):
            assert accession == "PTEST1"
            return parse_alphafold_pdb(
                accession,
                PDB_TEXT,
                "https://example.test/model_v6.pdb",
                mutation,
            )

    app.dependency_overrides[get_alphafold_client] = lambda: StubAlphaFoldClient()
    try:
        response = client.post(
            "/api/v1/structure/alphafold",
            json={"accession": "PTEST1", "mutation": "R2H"},
        )
    finally:
        app.dependency_overrides.clear()

    assert response.status_code == 200
    payload = response.json()
    assert payload["structure"]["mean_plddt"] == 71.75
    assert payload["structure"]["mutation_site"]["confidence"] == "low"
    assert "not a pathogenicity" in payload["interpretation"]
    assert "geometric proximity" in payload["interpretation"]


def test_string_graph_endpoint_returns_typed_evidence_graph() -> None:
    class StubStringClient:
        def build_graph(
            self,
            identifier: str,
            species: int = 9606,
            limit: int = 10,
            required_score: int = 700,
        ):
            assert (identifier, species, limit, required_score) == ("TP53", 9606, 5, 800)
            return build_evidence_graph(identifier, species, NETWORK_RECORDS, ENRICHMENT_RECORDS)

    app.dependency_overrides[get_string_client] = lambda: StubStringClient()
    try:
        response = client.post(
            "/api/v1/graph/string",
            json={"identifier": "TP53", "organism_id": 9606, "limit": 5, "required_score": 800},
        )
    finally:
        app.dependency_overrides.clear()

    assert response.status_code == 200
    payload = response.json()
    assert payload["seed"] == "TP53"
    assert {node["type"] for node in payload["nodes"]} == {"protein", "process", "pathway"}
    assert any(edge["type"] == "interacts_with" for edge in payload["edges"])


def test_pubmed_endpoint_returns_articles_and_report_section() -> None:
    pubmed = PubMedClient(
        transport=httpx.MockTransport(
            lambda request: (
                httpx.Response(200, json={"esearchresult": {"idlist": ["12345", "67890"]}})
                if request.url.path.endswith("/esearch.fcgi")
                else httpx.Response(200, json=SUMMARY_PAYLOAD)
            )
        )
    )
    app.dependency_overrides[get_pubmed_client] = lambda: pubmed
    try:
        response = client.post(
            "/api/v1/literature/pubmed",
            json={"protein": "TP53", "context_terms": ["Cell cycle"], "limit": 2},
        )
    finally:
        app.dependency_overrides.clear()

    assert response.status_code == 200
    payload = response.json()
    assert len(payload["evidence"]["articles"]) == 2
    assert "[PMID 12345]" in payload["markdown_section"]
    assert "does not by itself validate" in payload["evidence"]["disclaimer"]


def test_comprehensive_report_endpoint_returns_server_generated_bundle() -> None:
    app.dependency_overrides[get_uniprot_client] = lambda: SynthesisUniProtClient()
    app.dependency_overrides[get_alphafold_client] = lambda: SynthesisAlphaFoldClient()
    app.dependency_overrides[get_string_client] = lambda: SynthesisStringClient()
    app.dependency_overrides[get_pubmed_client] = lambda: SynthesisPubMedClient()
    app.dependency_overrides[get_interpro_client] = lambda: SynthesisInterProClient()
    app.dependency_overrides[get_kegg_client] = lambda: SynthesisKeggClient()
    try:
        response = client.post(
            "/api/v1/report/comprehensive",
            json={"identifier": "TP53", "organism_id": 9606, "mutation": "R2H"},
        )
    finally:
        app.dependency_overrides.clear()

    assert response.status_code == 200
    payload = response.json()
    assert payload["coverage"]["score"] == 100
    assert payload["annotation"]["accession"] == "P04637"
    assert payload["structure"]["mutation_site"]["position"] == 2
    assert payload["kegg"]["pathways"][0]["pathway_id"] == "hsa04115"
    assert payload["cellular_processes"]["processes"][1]["canonical_id"] == "hsa04115"
    # Per ADR-0012 the cellular-process layer keeps its own list empty;
    # hypotheses live on the separate phenotype_literature object.
    assert payload["cellular_processes"]["phenotype_hypotheses"] == []
    assert payload["phenotype_literature"]["gene"] == "TP53"
    assert payload["phenotype_literature"]["counts_by_level"]["supports"] >= 1
    assert len(payload["phenotype_literature"]["hypotheses"]) >= 1
    assert payload["reasoning_chain"]["chain_type"] == "mutation_impact"
    assert payload["reasoning_chain"]["mutation"] == "R2H"
    assert payload["literature"]["articles"][0]["pmid"] == "12345"
    assert payload["report_markdown"].startswith("# Qiwen Bio comprehensive report")


def test_embedding_endpoint_uses_versioned_cache(tmp_path) -> None:
    service = EmbeddingService(FakeProvider(), EmbeddingCache(tmp_path))
    app.dependency_overrides[get_embedding_service] = lambda: service
    try:
        first = client.post("/api/v1/embedding/esm2", json={"sequence": "ACDE"})
        second = client.post("/api/v1/embedding/esm2", json={"sequence": "ACDE"})
    finally:
        app.dependency_overrides.clear()

    assert first.status_code == 200
    assert first.json()["cached"] is False
    assert second.json()["cached"] is True
    assert second.json()["vector"] == [4.0, 2.0, 3.0]


def test_embedding_endpoint_maps_model_loading_failure_to_503() -> None:
    class FailingEmbeddingService:
        def embed(self, sequence: str):
            raise ModelLoadError("model files unavailable")

    app.dependency_overrides[get_embedding_service] = lambda: FailingEmbeddingService()
    try:
        response = client.post("/api/v1/embedding/esm2", json={"sequence": "ACDE"})
    finally:
        app.dependency_overrides.clear()

    assert response.status_code == 503
    assert response.json()["detail"] == "model files unavailable"


def _domain_annotation() -> DomainAnnotation:
    entry = DomainEntry(
        accession="PF00870",
        name="P53 DNA-binding domain",
        source_database="pfam",
        entry_type="domain",
        integrated_accession="IPR011615",
        source_url="https://www.ebi.ac.uk/interpro/entry/pfam/PF00870/",
        locations=[DomainLocation(start=100, end=288, status="CONTINUOUS")],
        go_terms=[],
        overlaps_mutation=True,
    )
    return DomainAnnotation(
        protein_accession="P04637",
        protein_length=393,
        mutation_position=175,
        entries=[entry],
        mutation_overlaps=[entry],
        entry_count=1,
        location_count=1,
        source_urls=["https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/P04637/"],
    )


def test_interpro_endpoint_returns_domain_and_mutation_overlap() -> None:
    class StubInterProClient:
        def fetch(self, accession: str, mutation_position: int | None = None):
            assert (accession, mutation_position) == ("P04637", 175)
            return _domain_annotation()

    app.dependency_overrides[get_interpro_client] = lambda: StubInterProClient()
    try:
        response = client.post(
            "/api/v1/domains/interpro",
            json={"accession": "P04637", "mutation": "R175H"},
        )
    finally:
        app.dependency_overrides.clear()

    assert response.status_code == 200
    assert response.json()["entries"][0]["accession"] == "PF00870"
    assert response.json()["mutation_overlaps"][0]["locations"][0]["start"] == 100
    assert "does not by itself" in response.json()["disclaimer"]


def test_interpro_endpoint_maps_not_found_and_service_errors() -> None:
    class FailingInterProClient:
        def __init__(self, error: Exception) -> None:
            self.error = error

        def fetch(self, accession: str, mutation_position: int | None = None):
            raise self.error

    for error, expected_status in (
        (InterProNotFoundError("no domains"), 404),
        (InterProServiceError("service unavailable"), 502),
    ):
        app.dependency_overrides[get_interpro_client] = lambda error=error: FailingInterProClient(
            error
        )
        try:
            response = client.post(
                "/api/v1/domains/interpro", json={"accession": "P04637"}
            )
        finally:
            app.dependency_overrides.clear()
        assert response.status_code == expected_status


def _kegg_annotation() -> KeggPathwayAnnotation:
    return KeggPathwayAnnotation(
        protein_accession="P04637",
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


def test_kegg_endpoint_returns_direct_pathway_records() -> None:
    class StubKeggClient:
        def fetch(self, accession: str, limit: int = 20):
            assert (accession, limit) == ("P04637", 5)
            return _kegg_annotation()

    app.dependency_overrides[get_kegg_client] = lambda: StubKeggClient()
    try:
        response = client.post(
            "/api/v1/pathways/kegg", json={"accession": "P04637", "limit": 5}
        )
    finally:
        app.dependency_overrides.clear()

    assert response.status_code == 200
    assert response.json()["pathways"][0]["pathway_id"] == "hsa04115"
    assert response.json()["pathways"][0]["source_url"].endswith("/hsa04115")
    assert "not evidence" in response.json()["disclaimer"]


def test_kegg_endpoint_maps_not_found_and_service_errors() -> None:
    class FailingKeggClient:
        def __init__(self, error: Exception) -> None:
            self.error = error

        def fetch(self, accession: str, limit: int = 20):
            raise self.error

    for error, expected_status in (
        (KeggNotFoundError("no pathways"), 404),
        (KeggServiceError("service unavailable"), 502),
    ):
        app.dependency_overrides[get_kegg_client] = lambda error=error: FailingKeggClient(
            error
        )
        try:
            response = client.post(
                "/api/v1/pathways/kegg", json={"accession": "P04637"}
            )
        finally:
            app.dependency_overrides.clear()
        assert response.status_code == expected_status


def test_phenotype_literature_endpoint_returns_gated_evidence() -> None:
    app.dependency_overrides[get_uniprot_client] = lambda: SynthesisUniProtClient()
    app.dependency_overrides[get_string_client] = lambda: SynthesisStringClient()
    app.dependency_overrides[get_pubmed_client] = lambda: SynthesisPubMedClient()
    app.dependency_overrides[get_kegg_client] = lambda: SynthesisKeggClient()
    try:
        response = client.post(
            "/api/v1/literature/phenotype",
            json={"identifier": "TP53", "organism_id": 9606, "limit_per_process": 2},
        )
    finally:
        app.dependency_overrides.clear()

    assert response.status_code == 200
    payload = response.json()
    evidence = payload["evidence"]
    assert evidence["gene"] == "TP53"
    assert evidence["counts_by_level"]["supports"] >= 1
    assert len(evidence["hypotheses"]) >= 1
    assert "abstract-level hypothesis" in evidence["hypotheses"][0]
    assert "## Phenotype literature evidence" in payload["markdown_section"]
    assert "abstract-metadata boundary" in payload["markdown_section"]


def test_phenotype_literature_endpoint_maps_uniprot_failure_to_404() -> None:
    from qiwen_bio.uniprot import ProteinNotFoundError

    class FailingUniProtClient:
        def resolve(self, identifier: str, organism_id: int = 9606):
            raise ProteinNotFoundError("not found")

    app.dependency_overrides[get_uniprot_client] = lambda: FailingUniProtClient()
    try:
        response = client.post(
            "/api/v1/literature/phenotype",
            json={"identifier": "NOSUCH", "organism_id": 9606},
        )
    finally:
        app.dependency_overrides.clear()

    assert response.status_code == 404


def test_reasoning_chain_endpoint_returns_mutation_impact_chain() -> None:
    app.dependency_overrides[get_uniprot_client] = lambda: SynthesisUniProtClient()
    app.dependency_overrides[get_alphafold_client] = lambda: SynthesisAlphaFoldClient()
    app.dependency_overrides[get_string_client] = lambda: SynthesisStringClient()
    app.dependency_overrides[get_pubmed_client] = lambda: SynthesisPubMedClient()
    app.dependency_overrides[get_interpro_client] = lambda: SynthesisInterProClient()
    app.dependency_overrides[get_kegg_client] = lambda: SynthesisKeggClient()
    try:
        response = client.post(
            "/api/v1/reasoning/chain",
            json={"identifier": "TP53", "organism_id": 9606, "mutation": "R2H"},
        )
    finally:
        app.dependency_overrides.clear()

    assert response.status_code == 200
    payload = response.json()
    chain = payload["chain"]
    assert chain["chain_type"] == "mutation_impact"
    assert chain["gene"] == "TP53"
    assert chain["mutation"] == "R2H"
    assert [step["step_id"] for step in chain["steps"]] == [
        "mutation", "structure", "function", "pathway", "phenotype"
    ]
    assert "## Reasoning chain" in payload["markdown_section"]
    assert "假设" in payload["markdown_section"] or "no hypothesis" in payload["markdown_section"]


def test_reasoning_chain_endpoint_maps_uniprot_failure_to_404() -> None:
    from qiwen_bio.uniprot import ProteinNotFoundError

    class FailingUniProtClient:
        def resolve(self, identifier: str, organism_id: int = 9606):
            raise ProteinNotFoundError("not found")

    app.dependency_overrides[get_uniprot_client] = lambda: FailingUniProtClient()
    try:
        response = client.post(
            "/api/v1/reasoning/chain",
            json={"identifier": "NOSUCH", "organism_id": 9606},
        )
    finally:
        app.dependency_overrides.clear()

    assert response.status_code == 404

import httpx
from fastapi.testclient import TestClient

from qiwen_bio.api import (
    app,
    get_alphafold_client,
    get_embedding_service,
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
)
from qiwen_bio.embedding import EmbeddingCache, EmbeddingService, ModelLoadError
from tests.test_embedding import FakeProvider


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

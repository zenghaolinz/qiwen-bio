from fastapi.testclient import TestClient

from qiwen_bio.api import app
from qiwen_bio.api import get_uniprot_client
from qiwen_bio.uniprot import parse_uniprot_record
from tests.test_uniprot import UNIPROT_RECORD


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

from fastapi.testclient import TestClient

from qiwen_bio.api import app


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

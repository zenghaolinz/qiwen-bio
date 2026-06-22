import httpx
import pytest

from qiwen_bio.interpro import InterProClient, InterProServiceError


def _result(
    accession: str,
    source: str,
    name: str,
    start: int,
    end: int,
    *,
    integrated: str | None = None,
    model: str | None = None,
    score: float | None = None,
) -> dict:
    return {
        "metadata": {
            "accession": accession,
            "name": name,
            "source_database": source,
            "type": "domain" if source == "pfam" else "family",
            "integrated": integrated,
            "go_terms": [],
        },
        "proteins": [
            {
                "accession": "p04637",
                "protein_length": 393,
                "entry_protein_locations": [
                    {
                        "fragments": [
                            {"start": start, "end": end, "dc-status": "CONTINUOUS"}
                        ],
                        "model": model,
                        "score": score,
                    }
                ],
            }
        ],
    }


def test_client_combines_cursor_paginated_interpro_and_pfam_entries() -> None:
    requested: list[str] = []

    def handler(request: httpx.Request) -> httpx.Response:
        requested.append(str(request.url))
        if "/entry/interpro/" in request.url.path and "cursor" not in request.url.params:
            return httpx.Response(
                200,
                json={
                    "count": 2,
                    "next": "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/P04637/?cursor=page2",
                    "results": [
                        _result("IPR0001", "interpro", "p53 family", 100, 200)
                    ],
                },
            )
        if "/entry/interpro/" in request.url.path:
            return httpx.Response(
                200,
                json={
                    "count": 2,
                    "next": None,
                    "results": [
                        _result("IPR0002", "interpro", "N terminus", 1, 50)
                    ],
                },
            )
        return httpx.Response(
            200,
            json={
                "count": 1,
                "next": None,
                "results": [
                    _result(
                        "PF00870",
                        "pfam",
                        "P53 DNA-binding domain",
                        120,
                        180,
                        integrated="IPR0001",
                        model="PF00870",
                        score=1.1e-59,
                    )
                ],
            },
        )

    with httpx.Client(transport=httpx.MockTransport(handler)) as http_client:
        annotation = InterProClient(http_client=http_client).fetch(
            "P04637", mutation_position=175
        )

    assert len(requested) == 3
    assert any("/entry/interpro/" in url for url in requested)
    assert any("/entry/pfam/" in url for url in requested)
    assert annotation.protein_length == 393
    assert {entry.source_database for entry in annotation.entries} == {"interpro", "pfam"}
    assert annotation.entry_count == 3
    assert annotation.location_count == 3
    assert {entry.accession for entry in annotation.mutation_overlaps} == {
        "IPR0001",
        "PF00870",
    }
    pfam = next(entry for entry in annotation.entries if entry.accession == "PF00870")
    assert pfam.integrated_accession == "IPR0001"
    assert pfam.locations[0].model == "PF00870"
    assert pfam.locations[0].score == 1.1e-59


def test_client_rejects_cursor_pagination_loop() -> None:
    loop_url = "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/P04637/"

    def handler(request: httpx.Request) -> httpx.Response:
        if "/entry/interpro/" in request.url.path:
            return httpx.Response(200, json={"count": 0, "next": loop_url, "results": []})
        return httpx.Response(200, json={"count": 0, "next": None, "results": []})

    with httpx.Client(transport=httpx.MockTransport(handler)) as http_client:
        with pytest.raises(InterProServiceError, match="pagination loop"):
            InterProClient(http_client=http_client).fetch("P04637")


def test_client_rejects_pagination_url_outside_official_api() -> None:
    def handler(request: httpx.Request) -> httpx.Response:
        if request.url.host == "www.ebi.ac.uk" and "/entry/interpro/" in request.url.path:
            return httpx.Response(
                200,
                json={
                    "count": 0,
                    "next": "http://127.0.0.1:8080/private",
                    "results": [],
                },
            )
        return httpx.Response(200, json={"count": 0, "next": None, "results": []})

    with httpx.Client(transport=httpx.MockTransport(handler)) as http_client:
        with pytest.raises(InterProServiceError, match="unsafe pagination URL"):
            InterProClient(http_client=http_client).fetch("P04637")

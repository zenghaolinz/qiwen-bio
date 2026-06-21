import httpx

from qiwen_bio.pubmed import PubMedClient
from qiwen_bio.reporting import render_literature_section


SUMMARY_PAYLOAD = {
    "result": {
        "uids": ["12345", "67890"],
        "12345": {
            "uid": "12345",
            "title": "TP53 and cell-cycle control.",
            "fulljournalname": "Journal of Example Biology",
            "pubdate": "2025 Jan",
            "authors": [{"name": "Smith A"}, {"name": "Chen B"}],
            "articleids": [
                {"idtype": "pubmed", "value": "12345"},
                {"idtype": "doi", "value": "10.1000/example.1"},
            ],
        },
        "67890": {
            "uid": "67890",
            "title": "A second p53 study.",
            "fulljournalname": "Evidence Reports",
            "pubdate": "2024 Dec",
            "authors": [{"name": "Garcia C"}],
            "articleids": [{"idtype": "pubmed", "value": "67890"}],
        },
    }
}


def test_client_builds_context_query_and_parses_article_metadata() -> None:
    requests: list[httpx.Request] = []

    def handler(request: httpx.Request) -> httpx.Response:
        requests.append(request)
        if request.url.path.endswith("/esearch.fcgi"):
            term = request.url.params["term"]
            assert '"TP53"[Title/Abstract]' in term
            assert '"Cell cycle"[Title/Abstract]' in term
            return httpx.Response(200, json={"esearchresult": {"idlist": ["12345", "67890"]}})
        if request.url.path.endswith("/esummary.fcgi"):
            assert request.url.params["id"] == "12345,67890"
            return httpx.Response(200, json=SUMMARY_PAYLOAD)
        return httpx.Response(404)

    client = PubMedClient(transport=httpx.MockTransport(handler))
    evidence = client.search("TP53", context_terms=["Cell cycle"], limit=5)

    assert len(requests) == 2
    assert evidence.protein == "TP53"
    assert evidence.context_terms == ["Cell cycle"]
    assert evidence.articles[0].pmid == "12345"
    assert evidence.articles[0].authors == ["Smith A", "Chen B"]
    assert evidence.articles[0].doi == "10.1000/example.1"
    assert evidence.articles[0].url == "https://pubmed.ncbi.nlm.nih.gov/12345/"


def test_empty_search_returns_auditable_empty_evidence() -> None:
    client = PubMedClient(
        transport=httpx.MockTransport(
            lambda _: httpx.Response(200, json={"esearchresult": {"idlist": []}})
        )
    )

    evidence = client.search("NORESULT", limit=3)

    assert evidence.articles == []
    assert "NORESULT" in evidence.query


def test_literature_markdown_links_each_pubmed_record() -> None:
    client = PubMedClient(
        transport=httpx.MockTransport(
            lambda request: (
                httpx.Response(200, json={"esearchresult": {"idlist": ["12345", "67890"]}})
                if request.url.path.endswith("/esearch.fcgi")
                else httpx.Response(200, json=SUMMARY_PAYLOAD)
            )
        )
    )
    evidence = client.search("TP53", limit=2)

    markdown = render_literature_section(evidence)

    assert "[PMID 12345](https://pubmed.ncbi.nlm.nih.gov/12345/)" in markdown
    assert "TP53 and cell-cycle control" in markdown
    assert "does not by itself validate" in markdown

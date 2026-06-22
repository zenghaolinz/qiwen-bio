import httpx

from qiwen_bio.pubmed import PubMedClient, PubMedServiceError
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


EFETCH_XML = """<?xml version="1.0"?>
<PubmedArticleSet>
<PubmedArticle><MedlineCitation><PMID>12345</PMID><Article>
<Journal><Title>Journal of Example Biology</Title></Journal>
<ArticleTitle>TP53 and cell-cycle control.</ArticleTitle>
<Abstract>
<AbstractText Label="BACKGROUND">TP53 is a tumor suppressor.</AbstractText>
<AbstractText Label="RESULTS">TP53 regulates the p53 signaling pathway in cell cycle arrest.</AbstractText>
</Abstract>
</Article></MedlineCitation></PubmedArticle>
<PubmedArticle><MedlineCitation><PMID>67890</PMID><Article>
<Journal><Title>Evidence Reports</Title></Journal>
<ArticleTitle>A second p53 study.</ArticleTitle>
</Article></MedlineCitation></PubmedArticle>
</PubmedArticleSet>"""


def test_search_with_abstracts_fetches_efetch_xml_and_parses_sections() -> None:
    requests: list[httpx.Request] = []

    def handler(request: httpx.Request) -> httpx.Response:
        requests.append(request)
        path = request.url.path
        if path.endswith("/esearch.fcgi"):
            return httpx.Response(200, json={"esearchresult": {"idlist": ["12345", "67890"]}})
        if path.endswith("/esummary.fcgi"):
            return httpx.Response(200, json=SUMMARY_PAYLOAD)
        if path.endswith("/efetch.fcgi"):
            assert request.url.params["db"] == "pubmed"
            assert request.url.params["retmode"] == "xml"
            assert request.url.params["id"] == "12345,67890"
            return httpx.Response(200, text=EFETCH_XML, headers={"Content-Type": "application/xml"})
        return httpx.Response(404)

    client = PubMedClient(transport=httpx.MockTransport(handler))
    evidence = client.search_with_abstracts(
        "TP53", context_terms=["Cell cycle"], limit=2
    )

    assert [req.url.path for req in requests] == [
        "/entrez/eutils/esearch.fcgi",
        "/entrez/eutils/esummary.fcgi",
        "/entrez/eutils/efetch.fcgi",
    ]
    by_pmid = {article.pmid: article for article in evidence.articles}
    assert by_pmid["12345"].abstract is not None
    assert "regulates the p53 signaling pathway" in by_pmid["12345"].abstract
    assert by_pmid["12345"].abstract_sections[0] == ("BACKGROUND", "TP53 is a tumor suppressor.")
    assert by_pmid["12345"].abstract_sections[1][0] == "RESULTS"
    assert "regulates" in by_pmid["12345"].abstract_sections[1][1]
    assert by_pmid["67890"].abstract is None
    assert by_pmid["67890"].abstract_sections == []


def test_search_with_abstracts_maps_efetch_failure_to_service_error() -> None:
    def handler(request: httpx.Request) -> httpx.Response:
        path = request.url.path
        if path.endswith("/esearch.fcgi"):
            return httpx.Response(200, json={"esearchresult": {"idlist": ["12345"]}})
        if path.endswith("/esummary.fcgi"):
            return httpx.Response(200, json=SUMMARY_PAYLOAD)
        if path.endswith("/efetch.fcgi"):
            return httpx.Response(500)
        return httpx.Response(404)

    client = PubMedClient(transport=httpx.MockTransport(handler))
    try:
        client.search_with_abstracts("TP53", limit=1)
    except PubMedServiceError:
        return
    raise AssertionError("expected PubMedServiceError on efetch failure")

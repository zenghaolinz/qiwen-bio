import httpx

from qiwen_bio.stringdb import StringClient, build_evidence_graph


NETWORK_RECORDS = [
    {
        "stringId_A": "9606.ENSP_TP53",
        "stringId_B": "9606.ENSP_MDM2",
        "preferredName_A": "TP53",
        "preferredName_B": "MDM2",
        "ncbiTaxonId": "9606",
        "score": 0.982,
    },
    {
        "stringId_A": "9606.ENSP_TP53",
        "stringId_B": "9606.ENSP_ATM",
        "preferredName_A": "TP53",
        "preferredName_B": "ATM",
        "ncbiTaxonId": "9606",
        "score": 0.941,
    },
    {
        "stringId_A": "9606.ENSP_MDM2",
        "stringId_B": "9606.ENSP_TP53",
        "preferredName_A": "MDM2",
        "preferredName_B": "TP53",
        "ncbiTaxonId": "9606",
        "score": 0.982,
    },
]

ENRICHMENT_RECORDS = [
    {
        "category": "Process",
        "term": "GO:0072331",
        "description": "Signal transduction by p53 class mediator",
        "fdr": 1.2e-8,
        "preferredNames": ["TP53", "MDM2", "ATM"],
    },
    {
        "category": "KEGG",
        "term": "hsa04115",
        "description": "p53 signaling pathway",
        "fdr": 2.4e-7,
        "preferredNames": ["TP53", "MDM2", "ATM"],
    },
    {
        "category": "COMPARTMENTS",
        "term": "GOCC:0005634",
        "description": "Nucleus",
        "fdr": 1e-10,
        "preferredNames": ["TP53"],
    },
]


def test_graph_normalizes_proteins_interactions_and_supported_terms() -> None:
    graph = build_evidence_graph("TP53", 9606, NETWORK_RECORDS, ENRICHMENT_RECORDS)

    assert graph.seed == "TP53"
    assert len([node for node in graph.nodes if node.type == "protein"]) == 3
    assert len([edge for edge in graph.edges if edge.type == "interacts_with"]) == 2
    assert {node.type for node in graph.nodes} == {"protein", "process", "pathway"}
    interaction = next(edge for edge in graph.edges if edge.type == "interacts_with")
    assert interaction.score in {0.982, 0.941}
    assert interaction.source_url == "https://string-db.org/api/json/network"
    pathway = next(node for node in graph.nodes if node.external_id == "hsa04115")
    assert pathway.label == "p53 signaling pathway"
    assert pathway.fdr == 2.4e-7


def test_client_fetches_network_then_enriches_resolved_proteins() -> None:
    requests: list[httpx.Request] = []

    def handler(request: httpx.Request) -> httpx.Response:
        requests.append(request)
        if request.url.path == "/api/json/network":
            return httpx.Response(200, json=NETWORK_RECORDS)
        if request.url.path == "/api/json/enrichment":
            assert "TP53" in str(request.url) and "MDM2" in str(request.url)
            return httpx.Response(200, json=ENRICHMENT_RECORDS)
        return httpx.Response(404)

    client = StringClient(transport=httpx.MockTransport(handler))
    graph = client.build_graph("TP53", species=9606, limit=10, required_score=700)

    assert [request.url.path for request in requests] == [
        "/api/json/network",
        "/api/json/enrichment",
    ]
    assert len(graph.nodes) == 5

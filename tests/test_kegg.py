import httpx
import pytest

from qiwen_bio.kegg import KeggClient, KeggNotFoundError, KeggServiceError


FLAT_FILE = """ENTRY       hsa04115                    Pathway
NAME        p53 signaling pathway - Homo sapiens (human)
DESCRIPTION p53 activation is induced by stress signals and DNA damage.
            It can lead to cell cycle arrest or apoptosis.
CLASS       Cellular Processes; Cell growth and death
///
ENTRY       hsa05200                    Pathway
NAME        Pathways in cancer - Homo sapiens (human)
DESCRIPTION Overview of cancer-associated pathways.
CLASS       Human Diseases; Cancer: overview
///
"""


def test_client_maps_uniprot_and_parses_direct_pathway_records() -> None:
    requested: list[str] = []

    def handler(request: httpx.Request) -> httpx.Response:
        requested.append(request.url.path)
        if request.url.path.startswith("/conv/genes/"):
            return httpx.Response(200, text="up:P04637\thsa:7157\n")
        if request.url.path.startswith("/link/pathway/"):
            return httpx.Response(
                200,
                text=(
                    "hsa:7157\tpath:hsa04115\n"
                    "hsa:7157\tpath:hsa05200\n"
                    "hsa:7157\tpath:hsa04115\n"
                ),
            )
        return httpx.Response(200, text=FLAT_FILE)

    with httpx.Client(transport=httpx.MockTransport(handler)) as http_client:
        annotation = KeggClient(http_client=http_client).fetch("P04637", limit=10)

    assert annotation.gene_ids == ["hsa:7157"]
    assert annotation.pathway_count == 2
    assert annotation.linked_pathway_count == 2
    assert annotation.truncated is False
    assert [pathway.pathway_id for pathway in annotation.pathways] == [
        "hsa04115",
        "hsa05200",
    ]
    assert annotation.pathways[0].name.startswith("p53 signaling pathway")
    assert "cell cycle arrest or apoptosis" in annotation.pathways[0].description
    assert annotation.pathways[0].classes == [
        "Cellular Processes",
        "Cell growth and death",
    ]
    assert annotation.pathways[0].source_url == "https://www.kegg.jp/entry/hsa04115"
    assert requested == [
        "/conv/genes/uniprot:P04637",
        "/link/pathway/hsa:7157",
        "/get/hsa04115+hsa05200",
    ]


def test_client_batches_pathway_get_requests_at_ten_entries() -> None:
    get_batches: list[list[str]] = []
    pathway_ids = [f"hsa{index:05d}" for index in range(1, 13)]

    def handler(request: httpx.Request) -> httpx.Response:
        if request.url.path.startswith("/conv/genes/"):
            return httpx.Response(200, text="up:P04637\thsa:7157\n")
        if request.url.path.startswith("/link/pathway/"):
            return httpx.Response(
                200,
                text="".join(f"hsa:7157\tpath:{item}\n" for item in pathway_ids),
            )
        batch = request.url.path.removeprefix("/get/").split("+")
        get_batches.append(batch)
        return httpx.Response(
            200,
            text="".join(
                f"ENTRY       {item}                    Pathway\nNAME        {item} name\n///\n"
                for item in batch
            ),
        )

    with httpx.Client(transport=httpx.MockTransport(handler)) as http_client:
        annotation = KeggClient(http_client=http_client).fetch("P04637", limit=11)

    assert [len(batch) for batch in get_batches] == [10, 1]
    assert annotation.pathway_count == 11
    assert annotation.linked_pathway_count == 12
    assert annotation.truncated is True


def test_client_distinguishes_not_found_from_service_failure() -> None:
    with httpx.Client(
        transport=httpx.MockTransport(lambda request: httpx.Response(200, text=""))
    ) as http_client:
        with pytest.raises(KeggNotFoundError, match="No KEGG gene mapping"):
            KeggClient(http_client=http_client).fetch("P00000")

    with httpx.Client(
        transport=httpx.MockTransport(lambda request: httpx.Response(503, text="down"))
    ) as http_client:
        with pytest.raises(KeggServiceError, match="request failed"):
            KeggClient(http_client=http_client).fetch("P04637")

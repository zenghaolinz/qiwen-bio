import httpx

from qiwen_bio.dataset import UniProtAmpDatasetClient, write_prepared_dataset


TSV = """Entry\tEntry Name\tSequence\tLength\tKeywords
P12345\tAMP_ONE\tACDEFG\t6\tAntimicrobial; Secreted
Q98765\tAMP_TWO\tKKLL\t4\tAntimicrobial
"""


def test_uniprot_positive_query_and_tsv_parsing() -> None:
    captured: dict[str, str] = {}

    def handler(request: httpx.Request) -> httpx.Response:
        captured["query"] = request.url.params["query"]
        captured["fields"] = request.url.params["fields"]
        return httpx.Response(200, text=TSV, headers={"content-type": "text/tab-separated-values"})

    with httpx.Client(transport=httpx.MockTransport(handler)) as http_client:
        records = UniProtAmpDatasetClient(http_client=http_client).fetch_records(
            label=1, limit=2, seed=9
        )

    assert "reviewed:true" in captured["query"]
    assert "keyword:KW-0929" in captured["query"]
    assert captured["fields"] == "accession,id,sequence,length,keyword"
    assert [record.record_id for record in records] == ["P12345", "Q98765"]
    assert all(record.label == 1 for record in records)
    assert records[0].source_url == "https://www.uniprot.org/uniprotkb/P12345/entry"


def test_uniprot_negative_query_is_explicitly_a_proxy() -> None:
    captured: dict[str, str] = {}

    def handler(request: httpx.Request) -> httpx.Response:
        captured["query"] = request.url.params["query"]
        return httpx.Response(200, text=TSV)

    with httpx.Client(transport=httpx.MockTransport(handler)) as http_client:
        records = UniProtAmpDatasetClient(http_client=http_client).fetch_records(
            label=0, limit=1, seed=3
        )

    assert "NOT (keyword:KW-0929)" in captured["query"]
    assert len(records) == 1
    assert records[0].label == 0
    assert records[0].source == "UniProtKB reviewed proxy-negative"


def test_write_prepared_dataset_includes_policy_and_auditable_rows(tmp_path) -> None:
    responses = [TSV, TSV.replace("P12345", "N12345").replace("Q98765", "N98765")]

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, text=responses.pop(0))

    with httpx.Client(transport=httpx.MockTransport(handler)) as http_client:
        client = UniProtAmpDatasetClient(http_client=http_client)
        prepared = client.prepare_dataset(positive_limit=2, negative_limit=2, seed=5)

    csv_path = tmp_path / "samples.csv"
    manifest_path = tmp_path / "manifest.json"
    write_prepared_dataset(prepared, csv_path=csv_path, manifest_path=manifest_path)

    csv_text = csv_path.read_text(encoding="utf-8")
    manifest_text = manifest_path.read_text(encoding="utf-8")
    assert "sample_id,sequence,label,source_ids,sources,cluster_id,split" in csv_text
    assert '"positive_policy"' in manifest_text
    assert '"negative_policy"' in manifest_text
    assert "proxy-negative" in manifest_text
    assert '"source_url": "https://rest.uniprot.org/uniprotkb/search"' in manifest_text
    assert '"license_url": "https://www.uniprot.org/help/license"' in manifest_text
    assert "reviewed:true" in manifest_text

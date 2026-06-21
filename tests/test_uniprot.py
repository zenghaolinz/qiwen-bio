import httpx
import pytest

from qiwen_bio.uniprot import (
    AmbiguousProteinError,
    ProteinNotFoundError,
    UniProtClient,
    parse_uniprot_record,
)


UNIPROT_RECORD = {
    "primaryAccession": "P04637",
    "proteinDescription": {
        "recommendedName": {"fullName": {"value": "Cellular tumor antigen p53"}}
    },
    "genes": [{"geneName": {"value": "TP53"}, "synonyms": [{"value": "P53"}]}],
    "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
    "sequence": {"value": "MEEPQSDPSV", "length": 10},
    "comments": [
        {"commentType": "FUNCTION", "texts": [{"value": "Acts as a tumor suppressor."}]}
    ],
    "uniProtKBCrossReferences": [
        {
            "database": "GO",
            "id": "GO:0003677",
            "properties": [{"key": "GoTerm", "value": "F:DNA binding"}],
        },
        {"database": "AlphaFoldDB", "id": "P04637"},
    ],
}


def test_parse_uniprot_record_extracts_traceable_annotation() -> None:
    annotation = parse_uniprot_record(UNIPROT_RECORD)

    assert annotation.accession == "P04637"
    assert annotation.gene_names == ["TP53", "P53"]
    assert annotation.sequence == "MEEPQSDPSV"
    assert annotation.functions == ["Acts as a tumor suppressor."]
    assert annotation.go_terms[0].aspect == "molecular_function"
    assert annotation.alphafold_url == "https://alphafold.ebi.ac.uk/entry/P04637"
    assert annotation.source_url == "https://rest.uniprot.org/uniprotkb/P04637"


def test_client_resolves_exact_gene_for_selected_organism() -> None:
    def handler(request: httpx.Request) -> httpx.Response:
        assert request.url.path == "/uniprotkb/search"
        assert "gene_exact%3ATP53" in str(request.url)
        assert "organism_id%3A9606" in str(request.url)
        return httpx.Response(200, json={"results": [UNIPROT_RECORD]})

    client = UniProtClient(transport=httpx.MockTransport(handler))
    annotation = client.resolve("TP53", organism_id=9606)

    assert annotation.accession == "P04637"


def test_client_distinguishes_missing_and_ambiguous_results() -> None:
    missing = UniProtClient(
        transport=httpx.MockTransport(lambda _: httpx.Response(200, json={"results": []}))
    )
    with pytest.raises(ProteinNotFoundError):
        missing.resolve("NOTAGENE")

    ambiguous_payload = {"results": [UNIPROT_RECORD, {**UNIPROT_RECORD, "primaryAccession": "Q00001"}]}
    ambiguous = UniProtClient(
        transport=httpx.MockTransport(lambda _: httpx.Response(200, json=ambiguous_payload))
    )
    with pytest.raises(AmbiguousProteinError):
        ambiguous.resolve("TP53")

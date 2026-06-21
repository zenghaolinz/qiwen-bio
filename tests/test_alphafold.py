import httpx
import pytest

from qiwen_bio.alphafold import (
    AlphaFoldClient,
    MutationMismatchError,
    parse_alphafold_pdb,
)


PDB_TEXT = """\
ATOM      1  CA  MET A   1      10.000  10.000  10.000  1.00 42.00           C
ATOM      2  CA  ARG A   2      11.000  10.000  10.000  1.00 68.00           C
ATOM      3  CA  GLY A   3      12.000  10.000  10.000  1.00 82.00           C
ATOM      4  CA  LYS A   4      13.000  10.000  10.000  1.00 95.00           C
END
"""


def test_parse_pdb_summarizes_plddt_and_mutation_site() -> None:
    analysis = parse_alphafold_pdb(
        accession="PTEST1",
        pdb_text=PDB_TEXT,
        structure_url="https://example.test/model.pdb",
        mutation="R2H",
    )

    assert analysis.residue_count == 4
    assert analysis.mean_plddt == 71.75
    assert analysis.confidence_distribution.model_dump() == {
        "very_high": 0.25,
        "confident": 0.25,
        "low": 0.25,
        "very_low": 0.25,
    }
    assert analysis.mutation_site.position == 2
    assert analysis.mutation_site.plddt == 68.0
    assert analysis.mutation_site.confidence == "low"
    assert [(item.position, item.distance_angstrom) for item in analysis.mutation_neighborhood] == [
        (1, 1.0),
        (3, 1.0),
        (4, 2.0),
    ]
    assert [(item.residue_a, item.residue_b) for item in analysis.contact_map.contacts] == [(1, 4)]
    assert analysis.contact_map.threshold_angstrom == 8.0
    assert len(analysis.coordinates) == 4


def test_structure_without_mutation_has_contact_map_but_no_neighborhood() -> None:
    analysis = parse_alphafold_pdb("PTEST1", PDB_TEXT, "https://example.test/model.pdb")

    assert analysis.mutation_neighborhood == []
    assert analysis.contact_map.total_contacts == 1


def test_parse_pdb_rejects_mismatched_mutation_residue() -> None:
    with pytest.raises(MutationMismatchError, match="expected R at position 3, found G"):
        parse_alphafold_pdb("PTEST1", PDB_TEXT, "https://example.test/model.pdb", "R3H")


def test_client_uses_versioned_url_from_prediction_metadata() -> None:
    def handler(request: httpx.Request) -> httpx.Response:
        if request.url.path == "/api/prediction/PTEST1":
            return httpx.Response(
                200,
                json=[{
                    "entryId": "AF-PTEST1-F1",
                    "uniprotAccession": "PTEST1",
                    "pdbUrl": "https://files.example.test/AF-PTEST1-F1-model_v6.pdb",
                    "cifUrl": "https://files.example.test/AF-PTEST1-F1-model_v6.cif",
                    "latestVersion": 6,
                }],
            )
        if request.url.host == "files.example.test":
            return httpx.Response(200, text=PDB_TEXT)
        return httpx.Response(404)

    client = AlphaFoldClient(transport=httpx.MockTransport(handler))
    analysis = client.analyze("PTEST1", mutation="R2H")

    assert analysis.model_version == 6
    assert analysis.structure_url.endswith("model_v6.pdb")
    assert analysis.mutation_site.mutant == "H"

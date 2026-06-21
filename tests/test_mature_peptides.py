import httpx

from qiwen_bio.dataset import RawAmpRecord, prepare_amp_dataset, write_prepared_dataset
from qiwen_bio.mature_peptides import (
    MaturePeptide,
    MaturePeptideEntry,
    UniProtMaturePeptideClient,
    curate_mature_positive_dataset,
    parse_mature_peptides,
)
from qiwen_bio.mature_peptide_cli import build_parser


ENTRY = {
    "primaryAccession": "P12345",
    "sequence": {"value": "MACDEFGKKLL"},
    "features": [
        {
            "type": "Signal",
            "location": {
                "start": {"value": 1, "modifier": "EXACT"},
                "end": {"value": 2, "modifier": "EXACT"},
            },
        },
        {
            "type": "Peptide",
            "description": "Mature one",
            "location": {
                "start": {"value": 3, "modifier": "EXACT"},
                "end": {"value": 6, "modifier": "EXACT"},
            },
        },
        {
            "type": "Peptide",
            "description": "Mature two",
            "location": {
                "start": {"value": 8, "modifier": "EXACT"},
                "end": {"value": 11, "modifier": "EXACT"},
            },
        },
    ],
}


def test_parser_extracts_exact_one_based_closed_peptide_features() -> None:
    result = parse_mature_peptides(ENTRY, expected_accession="P12345")

    assert [peptide.sequence for peptide in result.peptides] == ["CDEF", "KKLL"]
    assert [(peptide.start, peptide.end) for peptide in result.peptides] == [(3, 6), (8, 11)]
    assert [peptide.description for peptide in result.peptides] == ["Mature one", "Mature two"]
    assert result.exclusions == []


def test_parser_audits_fuzzy_and_out_of_range_features() -> None:
    entry = {
        **ENTRY,
        "features": [
            {
                "type": "Peptide",
                "location": {
                    "start": {"value": 3, "modifier": "UNCERTAIN"},
                    "end": {"value": 6, "modifier": "EXACT"},
                },
            },
            {
                "type": "Peptide",
                "location": {
                    "start": {"value": 8, "modifier": "EXACT"},
                    "end": {"value": 99, "modifier": "EXACT"},
                },
            },
        ],
    }

    result = parse_mature_peptides(entry, expected_accession="P12345")

    assert result.peptides == []
    assert [exclusion.reason for exclusion in result.exclusions] == [
        "non_exact_coordinate",
        "out_of_range_coordinate",
    ]


def test_client_uses_official_accession_json_endpoint() -> None:
    captured: dict[str, str] = {}

    def handler(request: httpx.Request) -> httpx.Response:
        captured["url"] = str(request.url)
        return httpx.Response(200, json=ENTRY)

    with httpx.Client(transport=httpx.MockTransport(handler)) as http_client:
        result = UniProtMaturePeptideClient(http_client=http_client).fetch("P12345")

    assert captured["url"] == "https://rest.uniprot.org/uniprotkb/P12345.json"
    assert len(result.peptides) == 2


def _peptide(accession: str, sequence: str, feature_index: int = 1) -> MaturePeptide:
    return MaturePeptide(
        accession=accession,
        feature_index=feature_index,
        description=f"peptide {accession}",
        start=2,
        end=len(sequence) + 1,
        sequence=sequence,
        parent_sequence_sha256=f"parent-{accession}",
    )


class FakeMatureClient:
    def __init__(self) -> None:
        self.requested: list[str] = []

    def fetch(self, accession: str) -> MaturePeptideEntry:
        self.requested.append(accession)
        peptides = {
            "p1": [_peptide("p1", "CDEF")],
            "p2": [_peptide("p2", "CDEF"), _peptide("p2", "KKLL", 2)],
            "p3": [],
        }[accession]
        return MaturePeptideEntry(
            accession=accession,
            parent_length=20,
            peptides=peptides,
            exclusions=[],
        )


def test_batch_curation_merges_duplicates_and_keeps_readiness_gates_closed(tmp_path) -> None:
    prepared = prepare_amp_dataset(
        [
            RawAmpRecord("p1", "ACDE", 1, "test"),
            RawAmpRecord("p2", "KKKK", 1, "test"),
            RawAmpRecord("p3", "LLLL", 1, "test"),
            RawAmpRecord("n1", "VVVV", 0, "test"),
        ],
        seed=5,
    )
    baseline_csv = tmp_path / "baseline.csv"
    baseline_manifest = tmp_path / "baseline.json"
    write_prepared_dataset(
        prepared, csv_path=baseline_csv, manifest_path=baseline_manifest
    )
    output_csv = tmp_path / "mature.csv"
    output_manifest = tmp_path / "mature-manifest.json"
    client = FakeMatureClient()

    result = curate_mature_positive_dataset(
        baseline_csv=baseline_csv,
        baseline_manifest=baseline_manifest,
        output_csv=output_csv,
        output_manifest=output_manifest,
        client=client,
    )

    assert client.requested == ["p1", "p2", "p3"]
    assert result.manifest.input_positive_records == 3
    assert result.manifest.entries_with_exact_peptide == 2
    assert result.manifest.entries_without_exact_peptide == 1
    assert result.manifest.extracted_instances == 3
    assert result.manifest.unique_mature_sequences == 2
    assert result.manifest.duplicate_instances == 1
    assert result.manifest.exact_feature_entry_coverage == 2 / 3
    assert result.manifest.min_mature_length == 4
    assert result.manifest.max_mature_length == 4
    assert result.manifest.mean_mature_length == 4.0
    assert result.manifest.has_defensible_negatives is False
    assert result.manifest.has_independent_external_benchmark is False
    assert result.manifest.ready_for_model_replacement is False
    merged = next(item for item in result.sequences if item.sequence == "CDEF")
    assert merged.source_accessions == ["p1", "p2"]
    assert "mature_id,sequence,source_accessions,descriptions,instance_count" in output_csv.read_text(
        encoding="utf-8"
    )


def test_mature_peptide_cli_has_auditable_default_paths() -> None:
    args = build_parser().parse_args([])

    assert args.baseline_csv.name == "uniprot_amp_samples.csv"
    assert args.baseline_manifest.name == "uniprot_amp_baseline.json"
    assert args.output_csv.name == "uniprot_mature_amp_sequences.csv"
    assert args.output_manifest.name == "uniprot_mature_amp_audit.json"

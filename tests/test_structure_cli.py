"""Tests for the structure-feature export (Stage 2, task 6).

The export produces a flat record for a future controlled structure-enhanced
benchmark. It must carry structure features only — never pathogenicity or
functional-effect labels — and must be deterministic.
"""

import json

from qiwen_bio.structure_cli import to_export_record
from qiwen_bio.structure_features import build_structure_feature_summary
from tests.test_structure_features import _annotation, _domains, _structure


def _summary(mutation="R175H", site_plddt=92.0, overlap=True):
    return build_structure_feature_summary(
        annotation=_annotation(),
        structure=_structure(mutation_site_plddt=site_plddt, mutation=mutation),
        domains=_domains(overlap=overlap, mutation_position=175),
        mutation=mutation,
    )


def test_export_record_has_all_required_fields() -> None:
    record = to_export_record(_summary(), sequence_length=393)
    required = {
        "accession", "sequence_length", "mutation", "mean_plddt",
        "mutation_site_plddt", "mutation_confidence_band",
        "contact_count_8a", "neighbor_count_8a", "domain_overlap",
        "overlapping_domains", "label", "source",
    }
    assert required.issubset(record.keys())
    assert record["accession"] == "P04637"
    assert record["sequence_length"] == 393
    assert record["mutation"] == "R175H"
    assert record["mean_plddt"] == 78.2
    assert record["mutation_site_plddt"] == 92.0
    assert record["mutation_confidence_band"] == "very_high"
    assert record["contact_count_8a"] == 12
    assert record["neighbor_count_8a"] == 2
    assert record["domain_overlap"] is True
    assert record["label"] is None
    assert record["source"] == "alphafold_predicted"


def test_export_record_label_is_settable() -> None:
    record = to_export_record(_summary(), sequence_length=393, label=1)
    assert record["label"] == 1


def test_export_record_handles_no_structure() -> None:
    summary = build_structure_feature_summary(
        annotation=_annotation(),
        structure=None,
        domains=_domains(overlap=False, mutation_position=175),
        mutation="R175H",
    )
    record = to_export_record(summary, sequence_length=393)
    assert record["mean_plddt"] is None
    assert record["mutation_site_plddt"] is None
    assert record["mutation_confidence_band"] == "unavailable"
    assert record["contact_count_8a"] is None
    assert record["source"] == "none"


def test_export_record_is_json_serializable() -> None:
    record = to_export_record(_summary(), sequence_length=393)
    # Must round-trip through JSON for JSONL export.
    serialized = json.dumps(record)
    assert json.loads(serialized) == record


def test_export_record_carries_no_pathogenicity_or_functional_fields() -> None:
    record = to_export_record(_summary(), sequence_length=393)
    # The export must not carry pathogenicity/functional-effect labels.
    forbidden = {"pathogenicity", "functional_effect", "stability_change", "impact"}
    assert not (forbidden & record.keys())

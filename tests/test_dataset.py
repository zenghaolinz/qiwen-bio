from qiwen_bio.dataset import RawAmpRecord, global_sequence_identity, prepare_amp_dataset


def test_global_sequence_identity_handles_substitution_and_indel() -> None:
    assert global_sequence_identity("ACDE", "ACDE") == 1.0
    assert global_sequence_identity("ACDE", "ACDF") == 0.75
    assert global_sequence_identity("ACDE", "ACDDE") == 0.8


def test_preparation_deduplicates_same_label_and_removes_label_conflicts() -> None:
    records = [
        RawAmpRecord(record_id="p1", sequence="ACDEFG", label=1, source="test"),
        RawAmpRecord(record_id="p1-dup", sequence="ac defg", label=1, source="test"),
        RawAmpRecord(record_id="conflict-pos", sequence="KKLL", label=1, source="test"),
        RawAmpRecord(record_id="conflict-neg", sequence="KKLL", label=0, source="test"),
        RawAmpRecord(record_id="n1", sequence="VVVV", label=0, source="test"),
    ]

    prepared = prepare_amp_dataset(records, identity_threshold=0.8, seed=7)

    assert len(prepared.samples) == 2
    positive = next(sample for sample in prepared.samples if sample.label == 1)
    assert positive.source_ids == ["p1", "p1-dup"]
    assert prepared.manifest.input_records == 5
    assert prepared.manifest.duplicate_records == 1
    assert prepared.manifest.conflicting_sequences == 1
    assert prepared.conflicts[0].sequence == "KKLL"


def test_similarity_clusters_never_cross_splits_and_are_deterministic() -> None:
    records = [
        RawAmpRecord(record_id="p1", sequence="ACDEFG", label=1, source="test"),
        RawAmpRecord(record_id="p2", sequence="ACDEFA", label=1, source="test"),
        RawAmpRecord(record_id="p3", sequence="KKKKKK", label=1, source="test"),
        RawAmpRecord(record_id="p4", sequence="HHHHHH", label=1, source="test"),
        RawAmpRecord(record_id="n1", sequence="VVVVVV", label=0, source="test"),
        RawAmpRecord(record_id="n2", sequence="VVVVVA", label=0, source="test"),
        RawAmpRecord(record_id="n3", sequence="GGGGGG", label=0, source="test"),
        RawAmpRecord(record_id="n4", sequence="TTTTTT", label=0, source="test"),
    ]

    first = prepare_amp_dataset(records, identity_threshold=0.8, seed=11)
    second = prepare_amp_dataset(records, identity_threshold=0.8, seed=11)

    assignments = {sample.sample_id: sample.split for sample in first.samples}
    assert assignments == {sample.sample_id: sample.split for sample in second.samples}
    cluster_splits: dict[str, set[str]] = {}
    for sample in first.samples:
        cluster_splits.setdefault(sample.cluster_id, set()).add(sample.split)
    assert all(len(splits) == 1 for splits in cluster_splits.values())
    p1 = next(sample for sample in first.samples if sample.sample_id == "p1")
    p2 = next(sample for sample in first.samples if sample.sample_id == "p2")
    assert p1.cluster_id == p2.cluster_id
    assert first.manifest.split_counts["train"] + first.manifest.split_counts["validation"] + first.manifest.split_counts["test"] == 8
    for split in ("train", "validation", "test"):
        assert {sample.label for sample in first.samples if sample.split == split} == {0, 1}

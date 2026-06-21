import json
from hashlib import sha256

import pytest

from qiwen_bio.dataset import RawAmpRecord, prepare_amp_dataset, write_prepared_dataset
from qiwen_bio.embedding import ProteinEmbedding
from qiwen_bio.training import (
    CalibratedAmpModel,
    EmbeddingDataset,
    EmbeddingSample,
    load_embedding_dataset,
    precompute_embedding_dataset,
    render_model_card,
    train_calibrated_baseline,
)
from qiwen_bio.training_cli import build_parser


class FakeEmbeddingService:
    def __init__(self) -> None:
        self.sequences: list[str] = []

    def embed(self, sequence: str) -> ProteinEmbedding:
        self.sequences.append(sequence)
        vector = [float(len(sequence)), float(sequence.count("A"))]
        return ProteinEmbedding(
            model_id="fake/esm",
            model_revision="revision-1",
            pooling="mean",
            sequence_sha256=sha256(sequence.encode()).hexdigest(),
            sequence_length=len(sequence),
            dimension=2,
            vector=vector,
            cached=False,
        )


def _write_dataset(tmp_path):
    records = [
        RawAmpRecord("p1", "ACDE", 1, "test"),
        RawAmpRecord("p2", "KKKK", 1, "test"),
        RawAmpRecord("n1", "VVVV", 0, "test"),
        RawAmpRecord("n2", "GGGG", 0, "test"),
    ]
    prepared = prepare_amp_dataset(records, identity_threshold=0.8, seed=7)
    csv_path = tmp_path / "samples.csv"
    manifest_path = tmp_path / "manifest.json"
    write_prepared_dataset(prepared, csv_path=csv_path, manifest_path=manifest_path)
    return csv_path, manifest_path, prepared


def test_precompute_binds_vectors_to_dataset_and_model_metadata(tmp_path) -> None:
    csv_path, manifest_path, prepared = _write_dataset(tmp_path)
    output_path = tmp_path / "embeddings.json"
    service = FakeEmbeddingService()

    artifact = precompute_embedding_dataset(
        csv_path=csv_path,
        manifest_path=manifest_path,
        output_path=output_path,
        embedding_service=service,
    )

    assert service.sequences == [sample.sequence for sample in prepared.samples]
    assert artifact.dataset_sha256 == prepared.manifest.dataset_sha256
    assert artifact.model_id == "fake/esm"
    assert artifact.model_revision == "revision-1"
    assert artifact.dimension == 2
    assert [row.split for row in artifact.samples] == [sample.split for sample in prepared.samples]
    assert output_path.exists()
    assert json.loads(output_path.read_text(encoding="utf-8"))["embedding_sha256"] == artifact.embedding_sha256


def test_precompute_rejects_csv_that_no_longer_matches_manifest(tmp_path) -> None:
    csv_path, manifest_path, _ = _write_dataset(tmp_path)
    csv_path.write_text(
        csv_path.read_text(encoding="utf-8").replace("ACDE", "ACDF"),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="dataset SHA-256"):
        precompute_embedding_dataset(
            csv_path=csv_path,
            manifest_path=manifest_path,
            output_path=tmp_path / "embeddings.json",
            embedding_service=FakeEmbeddingService(),
        )


def _synthetic_embeddings(*, test_shift: float = 0.0) -> EmbeddingDataset:
    samples: list[EmbeddingSample] = []
    split_sizes = {"train": 12, "validation": 6, "test": 6}
    for split, size in split_sizes.items():
        for index in range(size):
            label = index % 2
            center = 2.0 if label else -2.0
            shift = test_shift if split == "test" else 0.0
            samples.append(
                EmbeddingSample(
                    sample_id=f"{split}-{index}",
                    sequence_sha256=f"sha-{split}-{index}",
                    label=label,
                    cluster_id=f"cluster-{split}-{index}",
                    split=split,
                    vector=[center + index * 0.01 + shift, center * 0.5],
                )
            )
    return EmbeddingDataset(
        dataset_sha256="dataset-sha",
        model_id="fake/esm",
        model_revision="revision-1",
        pooling="mean",
        dimension=2,
        embedding_sha256=f"embedding-sha-{test_shift}",
        samples=samples,
    )


def test_training_exports_calibrated_json_model_with_metrics(tmp_path) -> None:
    output_path = tmp_path / "model.json"

    result = train_calibrated_baseline(
        _synthetic_embeddings(), output_path=output_path, seed=13
    )
    loaded = CalibratedAmpModel.model_validate_json(output_path.read_text(encoding="utf-8"))

    assert result.model.dataset_sha256 == "dataset-sha"
    assert result.model.embedding_model_revision == "revision-1"
    assert 0 < result.model.threshold < 1
    assert result.validation_metrics.sample_count == 6
    assert result.test_metrics.sample_count == 6
    assert result.test_metrics.roc_auc == 1.0
    assert loaded.predict_probability([2.0, 1.0]) == pytest.approx(
        result.model.predict_probability([2.0, 1.0])
    )
    assert loaded.predict_probability([2.0, 1.0]) > loaded.predict_probability([-2.0, -1.0])


def test_test_split_cannot_change_fitted_parameters(tmp_path) -> None:
    original = train_calibrated_baseline(
        _synthetic_embeddings(), output_path=tmp_path / "first.json", seed=17
    ).model
    shifted = train_calibrated_baseline(
        _synthetic_embeddings(test_shift=100.0),
        output_path=tmp_path / "second.json",
        seed=17,
    ).model

    assert shifted.scaler_mean == pytest.approx(original.scaler_mean)
    assert shifted.scaler_scale == pytest.approx(original.scaler_scale)
    assert shifted.classifier_coefficients == pytest.approx(original.classifier_coefficients)
    assert shifted.classifier_intercept == pytest.approx(original.classifier_intercept)
    assert shifted.calibration_coefficient == pytest.approx(original.calibration_coefficient)
    assert shifted.calibration_intercept == pytest.approx(original.calibration_intercept)
    assert shifted.threshold == pytest.approx(original.threshold)


def test_embedding_loader_rejects_tampered_vector(tmp_path) -> None:
    csv_path, manifest_path, _ = _write_dataset(tmp_path)
    output_path = tmp_path / "embeddings.json"
    precompute_embedding_dataset(
        csv_path=csv_path,
        manifest_path=manifest_path,
        output_path=output_path,
        embedding_service=FakeEmbeddingService(),
    )
    payload = json.loads(output_path.read_text(encoding="utf-8"))
    payload["samples"][0]["vector"][0] += 1
    output_path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="embedding SHA-256"):
        load_embedding_dataset(output_path)


def test_model_card_reports_metrics_and_non_deployment_decision(tmp_path) -> None:
    result = train_calibrated_baseline(
        _synthetic_embeddings(), output_path=tmp_path / "model.json", seed=23
    )

    card = render_model_card(result.model)

    assert "# AMP ESM-2 Logistic Baseline Model Card" in card
    assert "Test ROC AUC" in card
    assert "Proxy-negative" in card
    assert "precursor proteins" in card
    assert "same validation split" in card
    assert "Do not replace the demo predictor" in card


def test_training_cli_exposes_precompute_and_train_commands() -> None:
    parser = build_parser()

    precompute = parser.parse_args(["precompute"])
    train = parser.parse_args(["train"])

    assert precompute.command == "precompute"
    assert precompute.output.name == "amp_esm2_embeddings.json"
    assert train.command == "train"
    assert train.model_output.name == "amp_esm2_logistic.json"
    assert train.model_card.name == "amp-esm2-logistic-v0.1.md"

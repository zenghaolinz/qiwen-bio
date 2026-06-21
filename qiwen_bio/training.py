from __future__ import annotations

import json
import math
from csv import DictReader
from hashlib import sha256
from pathlib import Path
from typing import Protocol

from pydantic import BaseModel

from qiwen_bio.dataset import PreparedAmpSample, Split, prepared_samples_sha256
from qiwen_bio.embedding import ProteinEmbedding


class EmbeddingServiceLike(Protocol):
    def embed(self, sequence: str) -> ProteinEmbedding: ...


class EmbeddingSample(BaseModel):
    sample_id: str
    sequence_sha256: str
    label: int
    cluster_id: str
    split: Split
    vector: list[float]


class EmbeddingDataset(BaseModel):
    dataset_sha256: str
    model_id: str
    model_revision: str
    pooling: str
    dimension: int
    embedding_sha256: str
    samples: list[EmbeddingSample]


class BinaryMetrics(BaseModel):
    sample_count: int
    positive_count: int
    negative_count: int
    roc_auc: float
    average_precision: float
    accuracy: float
    balanced_accuracy: float
    precision: float
    recall: float
    f1: float
    brier_score: float
    expected_calibration_error: float


class CalibratedAmpModel(BaseModel):
    model_name: str = "esm2_logistic_platt"
    model_version: str = "0.1.0-proxy"
    dataset_sha256: str
    embedding_sha256: str
    embedding_model_id: str
    embedding_model_revision: str
    embedding_pooling: str
    dimension: int
    seed: int
    scaler_mean: list[float]
    scaler_scale: list[float]
    classifier_coefficients: list[float]
    classifier_intercept: float
    calibration_coefficient: float
    calibration_intercept: float
    threshold: float
    validation_metrics: BinaryMetrics
    test_metrics: BinaryMetrics
    limitations: list[str]

    def predict_probability(self, vector: list[float]) -> float:
        if len(vector) != self.dimension:
            raise ValueError(f"expected embedding dimension {self.dimension}")
        standardized = [
            (value - mean) / scale
            for value, mean, scale in zip(vector, self.scaler_mean, self.scaler_scale)
        ]
        base_score = self.classifier_intercept + sum(
            coefficient * value
            for coefficient, value in zip(self.classifier_coefficients, standardized)
        )
        calibrated_score = (
            self.calibration_coefficient * base_score + self.calibration_intercept
        )
        if calibrated_score >= 0:
            return 1.0 / (1.0 + math.exp(-calibrated_score))
        exponential = math.exp(calibrated_score)
        return exponential / (1.0 + exponential)


class TrainingResult(BaseModel):
    model: CalibratedAmpModel
    validation_metrics: BinaryMetrics
    test_metrics: BinaryMetrics


def _read_prepared_samples(csv_path: Path) -> list[PreparedAmpSample]:
    with csv_path.open(encoding="utf-8", newline="") as source:
        return [
            PreparedAmpSample(
                sample_id=row["sample_id"],
                sequence=row["sequence"],
                label=int(row["label"]),
                source_ids=row["source_ids"].split("|") if row["source_ids"] else [],
                sources=row["sources"].split("|") if row["sources"] else [],
                cluster_id=row["cluster_id"],
                split=row["split"],
            )
            for row in DictReader(source)
        ]


def _embedding_digest(
    *,
    dataset_sha256: str,
    model_id: str,
    model_revision: str,
    pooling: str,
    dimension: int,
    samples: list[EmbeddingSample],
) -> str:
    payload = {
        "dataset_sha256": dataset_sha256,
        "model_id": model_id,
        "model_revision": model_revision,
        "pooling": pooling,
        "dimension": dimension,
        "samples": [sample.model_dump(mode="json") for sample in samples],
    }
    return sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def load_embedding_dataset(path: Path) -> EmbeddingDataset:
    artifact = EmbeddingDataset.model_validate_json(path.read_text(encoding="utf-8"))
    expected = _embedding_digest(
        dataset_sha256=artifact.dataset_sha256,
        model_id=artifact.model_id,
        model_revision=artifact.model_revision,
        pooling=artifact.pooling,
        dimension=artifact.dimension,
        samples=artifact.samples,
    )
    if artifact.embedding_sha256 != expected:
        raise ValueError("embedding artifact does not match embedding SHA-256")
    return artifact


def precompute_embedding_dataset(
    *,
    csv_path: Path,
    manifest_path: Path,
    output_path: Path,
    embedding_service: EmbeddingServiceLike,
) -> EmbeddingDataset:
    prepared_samples = _read_prepared_samples(csv_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    dataset_sha = prepared_samples_sha256(prepared_samples)
    if dataset_sha != manifest.get("dataset_sha256"):
        raise ValueError("prepared CSV does not match manifest dataset SHA-256")
    if not prepared_samples:
        raise ValueError("embedding dataset cannot be empty")

    embedded: list[tuple[PreparedAmpSample, ProteinEmbedding]] = []
    for sample in prepared_samples:
        embedded.append((sample, embedding_service.embed(sample.sequence)))
    first_result = embedded[0][1]
    for _, result in embedded:
        if (
            result.model_id != first_result.model_id
            or result.model_revision != first_result.model_revision
            or result.pooling != first_result.pooling
            or result.dimension != first_result.dimension
            or len(result.vector) != first_result.dimension
        ):
            raise ValueError("embedding provider returned inconsistent metadata")

    samples = [
        EmbeddingSample(
            sample_id=sample.sample_id,
            sequence_sha256=result.sequence_sha256,
            label=sample.label,
            cluster_id=sample.cluster_id,
            split=sample.split,
            vector=result.vector,
        )
        for sample, result in embedded
    ]
    embedding_sha = _embedding_digest(
        dataset_sha256=dataset_sha,
        model_id=first_result.model_id,
        model_revision=first_result.model_revision,
        pooling=first_result.pooling,
        dimension=first_result.dimension,
        samples=samples,
    )
    artifact = EmbeddingDataset(
        dataset_sha256=dataset_sha,
        model_id=first_result.model_id,
        model_revision=first_result.model_revision,
        pooling=first_result.pooling,
        dimension=first_result.dimension,
        embedding_sha256=embedding_sha,
        samples=samples,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    temporary.write_text(artifact.model_dump_json(indent=2), encoding="utf-8")
    temporary.replace(output_path)
    return artifact


def _require_valid_splits(artifact: EmbeddingDataset) -> None:
    seen_clusters: dict[str, Split] = {}
    for sample in artifact.samples:
        if len(sample.vector) != artifact.dimension:
            raise ValueError("embedding vector dimension mismatch")
        previous = seen_clusters.setdefault(sample.cluster_id, sample.split)
        if previous != sample.split:
            raise ValueError("homology cluster crosses dataset splits")
    for split in ("train", "validation", "test"):
        labels = {sample.label for sample in artifact.samples if sample.split == split}
        if labels != {0, 1}:
            raise ValueError(f"{split} split must contain both labels")


def _binary_metrics(labels, probabilities, threshold: float) -> BinaryMetrics:
    import numpy as np
    from sklearn.metrics import (
        accuracy_score,
        average_precision_score,
        balanced_accuracy_score,
        brier_score_loss,
        f1_score,
        precision_score,
        recall_score,
        roc_auc_score,
    )

    labels_array = np.asarray(labels, dtype=int)
    probabilities_array = np.asarray(probabilities, dtype=float)
    predictions = (probabilities_array >= threshold).astype(int)
    bin_ids = np.minimum((probabilities_array * 10).astype(int), 9)
    calibration_error = 0.0
    for bin_id in range(10):
        mask = bin_ids == bin_id
        if mask.any():
            calibration_error += float(mask.mean()) * abs(
                float(probabilities_array[mask].mean()) - float(labels_array[mask].mean())
            )
    return BinaryMetrics(
        sample_count=len(labels_array),
        positive_count=int(labels_array.sum()),
        negative_count=int(len(labels_array) - labels_array.sum()),
        roc_auc=float(roc_auc_score(labels_array, probabilities_array)),
        average_precision=float(average_precision_score(labels_array, probabilities_array)),
        accuracy=float(accuracy_score(labels_array, predictions)),
        balanced_accuracy=float(balanced_accuracy_score(labels_array, predictions)),
        precision=float(precision_score(labels_array, predictions, zero_division=0)),
        recall=float(recall_score(labels_array, predictions, zero_division=0)),
        f1=float(f1_score(labels_array, predictions, zero_division=0)),
        brier_score=float(brier_score_loss(labels_array, probabilities_array)),
        expected_calibration_error=calibration_error,
    )


def _select_threshold(labels, probabilities) -> float:
    import numpy as np
    from sklearn.metrics import balanced_accuracy_score, f1_score

    candidates = sorted({float(value) for value in probabilities} | {0.5})
    return max(
        candidates,
        key=lambda threshold: (
            balanced_accuracy_score(labels, np.asarray(probabilities) >= threshold),
            f1_score(labels, np.asarray(probabilities) >= threshold, zero_division=0),
            -abs(threshold - 0.5),
        ),
    )


def train_calibrated_baseline(
    artifact: EmbeddingDataset, *, output_path: Path, seed: int = 42
) -> TrainingResult:
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    _require_valid_splits(artifact)

    def split_arrays(split: Split):
        selected = [sample for sample in artifact.samples if sample.split == split]
        return (
            np.asarray([sample.vector for sample in selected], dtype=float),
            np.asarray([sample.label for sample in selected], dtype=int),
        )

    train_vectors, train_labels = split_arrays("train")
    validation_vectors, validation_labels = split_arrays("validation")
    test_vectors, test_labels = split_arrays("test")
    scaler = StandardScaler().fit(train_vectors)
    classifier = LogisticRegression(random_state=seed, max_iter=2000).fit(
        scaler.transform(train_vectors), train_labels
    )
    validation_scores = classifier.decision_function(scaler.transform(validation_vectors))
    calibrator = LogisticRegression(random_state=seed, max_iter=2000).fit(
        validation_scores.reshape(-1, 1), validation_labels
    )
    validation_probabilities = calibrator.predict_proba(
        validation_scores.reshape(-1, 1)
    )[:, 1]
    threshold = _select_threshold(validation_labels, validation_probabilities)
    test_scores = classifier.decision_function(scaler.transform(test_vectors))
    test_probabilities = calibrator.predict_proba(test_scores.reshape(-1, 1))[:, 1]
    validation_metrics = _binary_metrics(
        validation_labels, validation_probabilities, threshold
    )
    test_metrics = _binary_metrics(test_labels, test_probabilities, threshold)
    model = CalibratedAmpModel(
        dataset_sha256=artifact.dataset_sha256,
        embedding_sha256=artifact.embedding_sha256,
        embedding_model_id=artifact.model_id,
        embedding_model_revision=artifact.model_revision,
        embedding_pooling=artifact.pooling,
        dimension=artifact.dimension,
        seed=seed,
        scaler_mean=scaler.mean_.tolist(),
        scaler_scale=scaler.scale_.tolist(),
        classifier_coefficients=classifier.coef_[0].tolist(),
        classifier_intercept=float(classifier.intercept_[0]),
        calibration_coefficient=float(calibrator.coef_[0][0]),
        calibration_intercept=float(calibrator.intercept_[0]),
        threshold=threshold,
        validation_metrics=validation_metrics,
        test_metrics=test_metrics,
        limitations=[
            "Positive labels are UniProt KW-0929 annotations, not assay-level activity measurements.",
            "Positive records may be precursor proteins rather than mature antimicrobial peptides.",
            "Negative labels are proxy negatives and may include unannotated antimicrobial proteins.",
            "Calibration and threshold selection reuse the same small validation split.",
            "The dataset is too small for clinical, safety-critical, or experimental decisions.",
        ],
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    temporary.write_text(model.model_dump_json(indent=2), encoding="utf-8")
    temporary.replace(output_path)
    return TrainingResult(
        model=model,
        validation_metrics=validation_metrics,
        test_metrics=test_metrics,
    )


def render_model_card(model: CalibratedAmpModel) -> str:
    validation = model.validation_metrics
    test = model.test_metrics
    return f"""# AMP ESM-2 Logistic Baseline Model Card

## Status

Offline research baseline only. **Do not replace the demo predictor** with this artifact.

## Model

- Classifier: standardized logistic regression with Platt calibration
- Embedding: `{model.embedding_model_id}` at `{model.embedding_model_revision}`
- Pooling/dimension: `{model.embedding_pooling}` / {model.dimension}
- Dataset SHA-256: `{model.dataset_sha256}`
- Embedding SHA-256: `{model.embedding_sha256}`
- Validation-selected threshold: {model.threshold:.6f}
- Random seed: {model.seed}

## Evaluation

- Test ROC AUC: {test.roc_auc:.4f}
- Test balanced accuracy: {test.balanced_accuracy:.4f}

| Metric | Validation (n={validation.sample_count}) | Test (n={test.sample_count}) |
| --- | ---: | ---: |
| ROC AUC | {validation.roc_auc:.4f} | {test.roc_auc:.4f} |
| Average precision | {validation.average_precision:.4f} | {test.average_precision:.4f} |
| Balanced accuracy | {validation.balanced_accuracy:.4f} | {test.balanced_accuracy:.4f} |
| F1 | {validation.f1:.4f} | {test.f1:.4f} |
| Brier score | {validation.brier_score:.4f} | {test.brier_score:.4f} |
| Expected calibration error | {validation.expected_calibration_error:.4f} | {test.expected_calibration_error:.4f} |

The scaler and base classifier were fitted on train only. Platt calibration and the decision threshold were fitted on validation only. Test was used only for the metrics above.

## Label Policy and Limitations

- Positives are reviewed UniProt entries annotated with antimicrobial keyword `KW-0929`, not assay-level measurements.
- Positive records may be precursor proteins rather than mature antimicrobial peptides.
- Proxy-negative examples are reviewed entries without `KW-0929`; missing annotation does not prove missing antimicrobial activity.
- Platt calibration and threshold selection reuse the same validation split, so both may be optimistic.
- The 100-record dataset and small test split cannot support clinical, safety-critical, or experimental decisions.
- UniProt contents can change; reproduce against the recorded dataset and embedding digests.

## Deployment Decision

Rejected for product replacement. Keep the transparent demo predictor labelled as a heuristic until a larger curated dataset with defensible negatives and external evaluation is available.
"""

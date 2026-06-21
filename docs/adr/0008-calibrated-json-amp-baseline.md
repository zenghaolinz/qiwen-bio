# ADR-0008: Keep the Calibrated AMP Baseline Offline and Export Pure JSON

## Status

Accepted

## Context

Stage 1C2a produced a frozen homology-aware split with 100 UniProt records. The project needs to test the complete embedding-to-classifier workflow without turning a small proxy-labelled experiment into a product claim. Pickled scikit-learn objects couple inference to library versions and are unsafe to load from untrusted sources.

## Decision

Precompute mean-pooled 320-dimensional vectors with the pinned ESM-2 revision and bind the artifact to both dataset and embedding SHA-256 digests. Fit `StandardScaler` and logistic regression on train only. Fit a one-dimensional Platt calibrator and choose the operating threshold on validation only. Read test only for the final reported metrics.

Export means, scales, classifier coefficients, calibration coefficients, threshold, provenance, metrics, and limitations as validated JSON. Implement inference math without scikit-learn or pickle. Keep the model offline and do not replace the transparent demo predictor.

## Evidence

The frozen 13-record test split produced ROC AUC 0.8571, average precision 0.8819, balanced accuracy 0.7619, F1 0.7273, Brier score 0.1831, and expected calibration error 0.3439. These point estimates are exploratory because the test set is very small.

## Consequences

### Positive

- Train, validation, and test responsibilities are explicit and regression-tested.
- Dataset or vector tampering is detected before training.
- JSON inference is inspectable, portable, and avoids pickle loading.
- The model card is generated from the exact trained artifact.

### Negative

- ECE remains high and the small validation split serves both calibration and threshold selection.
- Keyword positives may be precursor proteins rather than mature antimicrobial peptides.
- Proxy negatives may contain unannotated antimicrobial proteins.
- The point metrics have high uncertainty and no external validation.

## Alternatives Considered

- Deploy when ROC AUC exceeds a fixed threshold: rejected because label validity and test size dominate the point estimate.
- Train a neural classification head: rejected because 69 training records do not support the added capacity.
- Serialize a scikit-learn pipeline with pickle/joblib: rejected for portability and loading-safety reasons.
- Fit calibration on test: rejected because it would destroy the held-out evaluation.

## Next Decision Gate

Reconsider API integration only after curating mature-peptide positives and defensible negatives, using scalable homology clustering, separating calibration from threshold selection, and evaluating an external dataset.

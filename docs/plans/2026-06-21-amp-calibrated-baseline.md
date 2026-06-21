# AMP Calibrated Baseline Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Train and audit a calibrated AMP baseline on the frozen homology-aware UniProt proxy dataset without replacing the demo predictor prematurely.

**Architecture:** Precompute pinned ESM-2 vectors for the exact dataset digest, then train a standardized logistic regression on the train split. Fit a one-dimensional Platt calibrator and select the operating threshold on validation only; evaluate the frozen test split once. Export inference parameters as validated JSON rather than pickle and publish a model card that foregrounds proxy-negative bias and sample size.

**Tech Stack:** Python 3.11+, Pydantic, pinned ESM-2/Transformers/PyTorch, NumPy, scikit-learn, pytest.

---

### Task 1: Frozen Embedding Dataset

**Files:**
- Create: `qiwen_bio/training.py`
- Create: `tests/test_training.py`
- Modify: `.gitignore`

1. Write failing tests for CSV parsing, dataset digest matching, provider calls, output metadata, and refusal of stale embedding artifacts.
2. Run `python -m pytest -q tests/test_training.py` and confirm missing imports fail.
3. Implement a provider-injected precompute function that preserves sample IDs, labels, clusters, and splits and writes an atomic JSON artifact.
4. Run the focused tests and commit only after green.

### Task 2: Calibrated Logistic Baseline

**Files:**
- Modify: `qiwen_bio/training.py`
- Modify: `tests/test_training.py`
- Modify: `pyproject.toml`

1. Write failing synthetic-data tests proving train-only fitting, validation-only calibration and threshold selection, JSON round-trip inference, deterministic output, and metric schema.
2. Implement standard scaling, logistic regression, Platt calibration, threshold selection, metrics, and Pydantic JSON model schemas.
3. Add a `training` optional dependency group with NumPy and scikit-learn.
4. Run focused and full tests.

### Task 3: Offline Training CLI and Real Run

**Files:**
- Create: `qiwen_bio/training_cli.py`
- Create: `models/amp_esm2_logistic.json`
- Create: `docs/model-cards/amp-esm2-logistic-v0.1.md`

1. Add `precompute` and `train` CLI subcommands.
2. Precompute all 100 pinned ESM-2 embeddings into an ignored local artifact.
3. Train once and generate the JSON model plus Markdown model card.
4. Check that model metadata references the exact dataset and embedding digests.
5. Do not wire the model into the API unless the model card and test evidence justify replacement.

### Task 4: Stage Documentation and Verification

**Files:**
- Modify: `AGENTS.md`
- Modify: `README.md`
- Create: `docs/adr/0008-calibrated-json-amp-baseline.md`

1. Record actual validation/test metrics and limitations.
2. Mark 1C2b complete only if all artifacts and tests are present; separately state whether production replacement was accepted or rejected.
3. Run `python -m pytest -q`, compilation, model load/inference smoke tests, `git diff --check`, and inspect ignored training data.
4. Commit the stage and push `main`.

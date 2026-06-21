# AMP ESM-2 Logistic Baseline Model Card

## Status

Offline research baseline only. **Do not replace the demo predictor** with this artifact.

## Model

- Classifier: standardized logistic regression with Platt calibration
- Embedding: `facebook/esm2_t6_8M_UR50D` at `c731040fcd8d73dceaa04b0a8e6329b345b0f5df`
- Pooling/dimension: `mean` / 320
- Dataset SHA-256: `6727e736c6fa8373a955a5976fa288b0aeb7b0f25f16752bf11f870f10492730`
- Embedding SHA-256: `829cb54f499bb6d1a8c2293d5ccf94a7118482e1aeb2aaab1389248c3e637975`
- Validation-selected threshold: 0.620453
- Random seed: 42

## Evaluation

- Test ROC AUC: 0.8571
- Test balanced accuracy: 0.7619

| Metric | Validation (n=18) | Test (n=13) |
| --- | ---: | ---: |
| ROC AUC | 0.7750 | 0.8571 |
| Average precision | 0.7426 | 0.8819 |
| Balanced accuracy | 0.7875 | 0.7619 |
| F1 | 0.7778 | 0.7273 |
| Brier score | 0.1951 | 0.1831 |
| Expected calibration error | 0.2005 | 0.3439 |

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

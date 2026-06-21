# Qiwen Bio Agent Guide

## Mission

Build an evidence-oriented system that connects protein sequence and structure to function, pathways, cellular state, and phenotype without presenting unsupported model output as biological fact.

## Required Workflow

1. Use test-driven development for behavior changes: add a failing test, verify the failure, implement, then run the full suite.
2. Keep calculated facts, database records, user input, and model output distinct in `EvidenceItem.evidence_type`.
3. Every biological conclusion must expose its source and confidence or be labelled as a hypothesis.
4. Do not describe the demo AMP heuristic as a trained or validated predictor.
5. External clients must use dependency injection so tests do not require live network access.
6. Verify important integrations once against their official public endpoint before completing a stage.
7. After every major project stage, update this file, `README.md`, tests, and the stage history below in the same commit.

## Current Architecture

- `qiwen_bio/api.py`: FastAPI boundary and external-service error mapping.
- `qiwen_bio/pipeline.py`: deterministic analysis orchestration and evidence assembly.
- `qiwen_bio/features.py`: deterministic sequence-property calculations.
- `qiwen_bio/predictors.py`: replaceable prediction protocol and demo baseline.
- `qiwen_bio/uniprot.py`: UniProt lookup and annotation parsing.
- `qiwen_bio/alphafold.py`: AlphaFold metadata/PDB retrieval and pLDDT analysis.
- `qiwen_bio/stringdb.py`: STRING network/enrichment retrieval and typed evidence-graph normalization.
- `qiwen_bio/pubmed.py`: NCBI E-utilities retrieval and auditable PubMed metadata parsing.
- `qiwen_bio/synthesis.py`: optional-layer orchestration, graceful degradation, and evidence coverage scoring.
- `qiwen_bio/embedding.py`: pinned embedding providers, content-addressed disk cache, and lazy ESM-2 loading.
- `qiwen_bio/dataset.py`: AMP curation, conflict auditing, global-identity clustering, grouped splits, and UniProt retrieval.
- `qiwen_bio/dataset_cli.py`: reproducible dataset preparation and manifest export.
- `qiwen_bio/training.py`: frozen embedding artifacts, split validation, calibrated training, pure-JSON inference, metrics, and model-card rendering.
- `qiwen_bio/training_cli.py`: offline ESM precomputation and baseline training commands.
- `qiwen_bio/mature_peptides.py`: exact UniProt Peptide-feature extraction, positive curation, exclusions, and readiness gates.
- `qiwen_bio/mature_peptide_cli.py`: reproducible mature-positive audit against the frozen AMP baseline.
- `qiwen_bio/reporting.py`: Markdown report rendering.
- `qiwen_bio/static/`: dependency-free Web Demo.

## Stage History

| Stage | Status | Delivered | Explicitly not delivered |
| --- | --- | --- | --- |
| 1A Engineering MVP | Complete | Sequence validation, deterministic features, API, Web Demo, evidence chain, Markdown report | ESM embeddings and trained classifier |
| 1B UniProt annotation | Complete | Reviewed entry lookup, sequence/function/GO extraction, provenance, AlphaFold entry discovery | InterPro, KEGG, STRING, PubMed |
| 1C1 Embedding foundation | Complete | Pinned ESM-2 8M provider, 320D mean pooling, versioned disk cache, API and Web trigger | Larger ESM variants, batching, GPU verification |
| 1C2a Dataset foundation | Complete | Reviewed UniProt AMP-keyword positives, explicit proxy negatives, deduplication, conflict removal, 80% identity clusters, leakage-safe grouped split targeting 70/15/15, CSV/manifest export | Confirmed negative set, scalable MMseqs2 clustering, trained classifier |
| 1C2b Calibrated offline baseline | Complete | Frozen pinned embeddings, train-only scaler/classifier, validation-only Platt calibration and threshold, pure-JSON model, model card, held-out metrics | Product/API replacement rejected because of proxy labels, precursor ambiguity, and n=13 test set |
| 1C2c Mature-positive curation | Complete | Exact UniProt Peptide-feature extraction, 50-accession audit, 26 unique mature sequences, coverage/length manifest, explicit readiness gates | Defensible negatives and independent external benchmark |
| 1C2d External dataset validation | Planned | Defensible negatives, scalable clustering, independent benchmark, separate calibration/threshold sets | Not started |
| 2A Structure confidence | Complete | Dynamic AlphaFold model URL lookup, PDB parsing, mean/distributed pLDDT, mutation-site confidence | Contact maps, secondary structure, pockets, 3D viewer, SaProt |
| 2B Structure context | Complete | CA contact map, 8 A mutation neighborhood, interactive dependency-free backbone viewer | Secondary structure, pockets, all-atom contacts, SaProt |
| 3A Interaction graph | Complete | STRING interactions, process/pathway enrichment, typed source-linked graph, deterministic Canvas view | Direct KEGG API, PubMed, phenotype reasoning |
| 3B Literature retrieval | Complete | Context-bound PubMed search, structured citations, PMID links, Markdown report section | Abstract/full-text appraisal, claim-level support classification |
| 3C Evidence synthesis | Complete | One server-generated report, six-layer coverage score, optional-service degradation, unified Web result | Claim correctness scoring, phenotype prediction |
| 4 Experimental imaging | Planned | Gel/PCR/microscopy analysis | Not started |

## Next Priority

Implement stage 1C2d: qualify a redistributable source with defensible negatives and an independent benchmark, scale homology clustering, and separate calibration from threshold selection before reconsidering product integration.

## Dataset Baseline

- Reproduce with `python -m qiwen_bio.dataset_cli --positive-limit 50 --negative-limit 50`.
- The committed manifest is `data/manifests/uniprot_amp_baseline.json`; sequence rows remain under ignored `data/datasets/`.
- The frozen baseline retains 100 samples in 95 clusters with balanced labels; indivisible clusters produce an actual 69/18/13 split against the 70/15/15 target.
- Positive labels mean reviewed UniProtKB entries carrying keyword `KW-0929`.
- Negative labels are proxy negatives: absence of `KW-0929` is not evidence of absent antimicrobial activity.
- Clustering uses deterministic Needleman-Wunsch identity and connected components. It is explainable but O(n^2); use MMseqs2 or equivalent before expanding beyond small baselines.

## Offline Model Baseline

- Reproduce embeddings with `python -m qiwen_bio.training_cli precompute` and training with `python -m qiwen_bio.training_cli train`.
- Model artifact: `models/amp_esm2_logistic.json`; model card: `docs/model-cards/amp-esm2-logistic-v0.1.md`.
- Frozen test metrics (n=13): ROC AUC 0.8571, balanced accuracy 0.7619, F1 0.7273, Brier 0.1831, ECE 0.3439.
- The artifact is an offline research baseline. Do not wire it into `ProteinPredictor` or the API without a new stage decision based on curated labels and external evaluation.

## Mature-Positive Audit

- Reproduce with `python -m qiwen_bio.mature_peptide_cli`.
- The committed manifest is `data/manifests/uniprot_mature_amp_audit.json`; sequence rows remain under ignored `data/datasets/`.
- All 50 positive accessions were fetched successfully; 20 (40%) contain exact UniProt `Peptide` features.
- The audit extracted 26 unique mature sequences of length 9-56 residues (mean 31.12).
- `has_defensible_negatives` and `has_independent_external_benchmark` remain false, so `ready_for_model_replacement` must remain false.

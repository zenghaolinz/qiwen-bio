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
- `qiwen_bio/interpro.py`: direct cursor-paginated InterPro/Pfam domain retrieval, coordinate normalization, and mutation overlap.
- `qiwen_bio/kegg.py`: fixed-host UniProt-to-KEGG mapping, direct pathway links, batched flat-file metadata, and truncation audit.
- `qiwen_bio/cellular_processes.py`: stable-ID GO/KEGG/STRING process normalization, seed-edge filtering, multi-source supports, and empty phenotype-hypothesis boundary.
- `qiwen_bio/phenotype_literature.py`: phenotype-oriented PubMed abstract retrieval, conservative abstract-metadata claim-support classification, gated hypotheses with explicit uncertainty.
- `qiwen_bio/reasoning_chain.py`: cross-scale reasoning chain assembly, layered mutation/function triggering, wild-type validation, gated "可能影响" hypotheses.
- `qiwen_bio/reporting.py`: Markdown report rendering.
- `qiwen_bio/static/`: dependency-free Web Demo.

## Stage History

| Stage | Status | Delivered | Explicitly not delivered |
| --- | --- | --- | --- |
| 1A Engineering MVP | Complete | Sequence validation, deterministic features, API, Web Demo, evidence chain, Markdown report | ESM embeddings and trained classifier |
| 1B UniProt annotation | Complete | Reviewed entry lookup, sequence/function/GO extraction, provenance, AlphaFold entry discovery | STRING and PubMed are separate stages |
| 1C1 Embedding foundation | Complete | Pinned ESM-2 8M provider, 320D mean pooling, versioned disk cache, API and Web trigger | Larger ESM variants, batching, GPU verification |
| 1C2a Dataset foundation | Complete | Reviewed UniProt AMP-keyword positives, explicit proxy negatives, deduplication, conflict removal, 80% identity clusters, leakage-safe grouped split targeting 70/15/15, CSV/manifest export | Confirmed negative set, scalable MMseqs2 clustering, trained classifier |
| 1C2b Calibrated offline baseline | Complete | Frozen pinned embeddings, train-only scaler/classifier, validation-only Platt calibration and threshold, pure-JSON model, model card, held-out metrics | Product/API replacement rejected because of proxy labels, precursor ambiguity, and n=13 test set |
| 1C2c Mature-positive curation | Complete | Exact UniProt Peptide-feature extraction, 50-accession audit, 26 unique mature sequences, coverage/length manifest, explicit readiness gates | Defensible negatives and independent external benchmark |
| 1C2d External dataset validation | Planned | Defensible negatives, scalable clustering, independent benchmark, separate calibration/threshold sets | Not started |
| 2A Structure confidence | Complete | Dynamic AlphaFold model URL lookup, PDB parsing, mean/distributed pLDDT, mutation-site confidence | Contact maps, secondary structure, pockets, 3D viewer, SaProt |
| 2B Structure context | Complete | CA contact map, 8 A mutation neighborhood, interactive dependency-free backbone viewer | Secondary structure, pockets, all-atom contacts, SaProt |
| 2C Domain context | Complete | Direct InterPro and Pfam entries, cursor pagination, fragment coordinates, mutation-position overlap, API and report layer | Domain-impact prediction, conservation scoring |
| 3A Interaction graph | Complete | STRING interactions, process/pathway enrichment, typed source-linked graph, deterministic Canvas view | PubMed and phenotype reasoning are separate stages |
| 3B Literature retrieval | Complete | Context-bound PubMed search, structured citations, PMID links, Markdown report section | Abstract/full-text appraisal, claim-level support classification |
| 3C Evidence synthesis | Complete | One server-generated report, seven-layer coverage score, optional-service degradation, unified Web result | Claim correctness scoring, phenotype prediction |
| 3D Direct pathway evidence | Complete | UniProt-to-KEGG gene conversion, direct pathway links, batched metadata, truncation audit, API/report integration, explicit STRING fallback | Pathway activity, directionality, causality, phenotype prediction |
| 3E Cellular-process evidence | Complete | Stable-ID GO/KEGG/STRING merge, seed-linked enrichment filtering, typed multi-source supports, direct-first literature context, API/report bundle | Process activity, effect direction, phenotype hypotheses, causality |
| 3F Phenotype literature | Complete | Phenotype-oriented PubMed abstract retrieval (efetch XML), conservative abstract-metadata claim-support classification (supports/mentions/no_abstract), gated hypotheses with explicit uncertainty, API endpoint, report section, Web panel | Full-text appraisal, claim-level support grading, causal/directional inference, phenotype ontology, LLM reasoning |
| 3G Cross-scale reasoning chain | Complete | Layered mutation-impact / protein-function chain (mutation→structure→function→pathway→phenotype), wild-type validation against UniProt sequence, gated "可能影响" hypotheses, missing-layer tracking, API endpoint, report section, Web panel | LLM-based reasoning, causal/directional inference, cellular-state measurement, full-text appraisal |
| 4 Experimental imaging | Planned | Gel/PCR/microscopy analysis | Not started |

## Next Priority

Stages 3F and 3G close the protein-process-phenotype literature loop and assemble the cross-scale reasoning chain. Remaining Stage 3 work is the configurable LLM reasoning layer and measured cellular-state data. In parallel, model deployment still requires defensible AMP negatives and an independent benchmark. Other major gaps are SaProt, controlled structure-enhancement evaluation, and experimental imaging.

## Cellular-Process Evidence

- The comprehensive response exposes `cellular_processes.processes` with stable IDs and typed supports.
- UniProt GO biological-process annotations, direct KEGG memberships, and seed-linked STRING enrichment merge without losing provenance.
- STRING terms are included only when an `annotated_to` edge originates from the seed protein, preventing neighbor-only enrichment from becoming seed evidence.
- TP53 live verification produced 92 normalized records with support counts UniProt GO 67, KEGG 20, STRING enrichment 7.
- `phenotype_hypotheses` on `CellularProcessEvidence` is intentionally empty (ADR-0012): the cellular-process layer organizes process evidence but never infers activity, direction, mechanism, causality, or phenotype. Stage 3F phenotype hypotheses live on the separate `PhenotypeLiteratureEvidence` object and are gated on explicit abstract-level `supports` evidence.

## Phenotype Literature Evidence

- Endpoint: `POST /api/v1/literature/phenotype` with an identifier, organism, and 1-5 `limit_per_process`.
- Only processes with **direct** support (UniProt GO annotation or direct KEGG membership) are hypothesis-eligible; enrichment-only STRING processes are literature-context only.
- For each eligible process the builder runs a phenotype-oriented PubMed query `("<GENE>"[Title/Abstract]) AND ("<process label>"[Title/Abstract])`, fetches abstracts via efetch (`retmode=xml`), and parses structured `AbstractText` sections including their `Label` attribute.
- The classifier labels each article at the **abstract-metadata boundary only** as `supports`, `mentions`, or `no_abstract`. `supports` requires all of: gene token present, substantive process-label token present, and a pinned conservative outcome verb (e.g. `regulates`, `associated with`, `activates`, `suppresses`, `required for`) in the title or abstract.
- A phenotype hypothesis is emitted only when a directly supported process has at least one `supports` article. Each hypothesis cites the direct support sources, the PMID(s), the matched verb(s), and an explicit uncertainty statement.
- TP53 live verification: `p53 signaling pathway` query returned PMIDs 36859359 and 31365877, both classified `supports` (verbs `activate`, `suppress`, `regulates`, `associated with`, `correlation`), so hypotheses were emitted.
- PubMed failures degrade gracefully: the process link is skipped and no hypothesis is emitted; the main literature layer surfaces its own warning.
- This is not full-text appraisal, claim-level evidence grading, causal inference, or a phenotype prediction.

## Reasoning Chain Evidence

- Endpoint: `POST /api/v1/reasoning/chain` with an identifier, organism, optional mutation, and interaction/literature limits.
- The chain is a **synthesis layer**: it performs no new network requests and introduces no new biological claims. It assembles the already-retrieved UniProt, AlphaFold, InterPro/Pfam, STRING, KEGG, cellular-process, and phenotype-literature layers into a single structured chain.
- Layered triggering: a supplied mutation (any non-empty string) produces a `mutation_impact` chain (mutation→structure→function→pathway→phenotype); no mutation produces a `protein_function` chain (structure→function→pathway→phenotype). A malformed mutation keeps the `mutation_impact` type but marks the mutation step unavailable.
- Wild-type residue is validated against the UniProt sequence. A mismatch is recorded as a fact and suppresses downstream "可能影响" hypotheses (facts remain, but no functional-impact inference is drawn).
- Hypothesis gating: each "假设：可能影响/可能相关" sentence requires explicit upstream evidence (domain overlap for function, function hypothesis for pathway, ADR-0013 `supports` article for phenotype). Every hypothesis carries an uncertainty clause.
- Missing layers (structure, kegg, string_graph, phenotype_literature) are listed in `missing_layers` and the step is marked `available=False`.
- Per ADR-0014: the chain organizes evidence; it does not establish activity, direction, mechanism, causality, pathogenicity, or a phenotype.

## Direct KEGG Evidence

- Endpoint: `POST /api/v1/pathways/kegg` with UniProt accession and a 1-50 result limit.
- Query chain: UniProt accession -> KEGG gene ID -> direct pathway links -> pathway flat-file metadata in batches of ten.
- TP53 live verification mapped `P04637` to `hsa:7157`, found 51 linked pathways, and returned the default first 20 with `truncated=true`.
- KEGG content is retrieved on demand for academic use and is not bundled, cached, or redistributed by the repository.
- Membership is an association, not evidence of pathway activation, directionality, or phenotype causality.

## Domain Evidence

- Endpoint: `POST /api/v1/domains/interpro` with UniProt accession and optional mutation.
- Both InterPro and Pfam source endpoints are queried directly with cursor pagination and location deduplication.
- TP53 `R175H` live verification returned 13 entries (9 InterPro, 4 Pfam) and five position-overlapping entries.
- Mutation overlap means only that the supplied residue number lies inside an annotated fragment; it is not a functional-effect or pathogenicity prediction.

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

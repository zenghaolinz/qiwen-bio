# Qiwen Bio Project Status

Updated: 2026-06-22

## Scoring Method

Progress is measured against the 24 tasks in the original four-stage proposal. Complete tasks score 1, partial tasks score 0.5, and not-started tasks score 0. The fourth-stage experimental-image work is part of the full vision but optional to the current core bioinformatics platform.

- Full four-stage vision: **15 / 24 = 62.5%**
- Core software scope (stages 1-3): **15 / 18 = 83.33%**

These percentages measure delivered scope, not biological correctness or model accuracy.

## Stage 1: Minimum Viable Version

Score: **6.5 / 7 = 92.86%**

| Proposal task | Status | Evidence |
| --- | --- | --- |
| Project backend | Complete | FastAPI, OpenAPI, tests, Windows launcher |
| Integrate ESM-C or ESM-2 | Complete | Pinned ESM-2 provider |
| Protein embedding extraction | Complete | 320D mean pooling, cache, API, offline batch artifact |
| Select a narrow task | Complete | AMP demonstration task |
| Train a lightweight classifier | Complete | Calibrated logistic offline baseline and JSON artifact |
| Prediction result and explanation | Partial | Transparent demo is integrated; trained model deployment is correctly rejected |
| Markdown report | Complete | Sequence and comprehensive evidence reports |

## Stage 2: Structure-Enhanced Version

Score: **2.5 / 5 = 50%**

| Proposal task | Status | Evidence / gap |
| --- | --- | --- |
| Fetch AlphaFold structure | Complete | Dynamic AlphaFold model/PDB retrieval |
| Extract structural features | Partial | pLDDT, CA contacts, mutation neighborhood complete; secondary structure and pockets missing |
| Structure visualization | Complete | Dependency-free interactive backbone/contact views |
| SaProt structure embedding | Not started | Foldseek/SaProt pipeline absent |
| Compare sequence and structure models | Not started | No controlled benchmark |

## Stage 3: Multiscale Reasoning

Score: **5.5 / 6 = 91.67%**

| Proposal task | Status | Evidence / gap |
| --- | --- | --- |
| UniProt, GO, KEGG, STRING integration | Complete | UniProt/GO/STRING, direct InterPro/Pfam, and direct KEGG pathway records are integrated |
| Local knowledge graph | Complete | Typed interaction/process/pathway graph and visualization |
| Protein-pathway-phenotype chain | Partial | A deterministic cross-scale evidence chain (stage 3G) assembles mutation→structure→function→pathway→phenotype with gated hypotheses and wild-type validation. This is evidence-chain synthesis, not a validated phenotype prediction, causal inference, or Graph RAG. Process activity and cellular-state measurement remain out of scope. |
| PubMed evidence | Complete | Context retrieval, auditable metadata/citations, and abstract retrieval with claim-support classification |
| LLM structured report | Partial | Deterministic structured synthesis and cross-scale chain exist; configurable LLM reasoning layer absent |
| Confidence and uncertainty | Complete | Evidence coverage, provenance, optional-layer warnings, model limitations |

### What the Stage 3 system is and is not

The current Stage 3 system is a **deterministic evidence-chain synthesis** layer. It is explicitly **not**:

- a Graph RAG system;
- a full LLM reasoning system (no configurable LLM reasoning layer is implemented);
- a biological causal inference engine (hypotheses are gated and uncertainty-bound, never causal claims);
- a validated phenotype prediction system;
- a trained deployable biological model (the AMP predictor remains a heuristic / offline baseline);
- a clinical-grade evidence system.

There is no LoRA / Adapter fine-tuning and no SaProt structure-aware embedding in the current branch.

## Stage 4: Experimental and Imaging Extensions

Score: **0 / 6 = 0%**

Gel analysis, PCR interpretation, microscopy, colony counting, experimental-report automation, and image-driven follow-up suggestions have not started.

## Current Model Gate

The engineering training path is complete, but production model integration remains blocked by data validity:

- Exact mature-positive annotations: available for 20/50 source accessions (26 unique sequences).
- Defensible negatives: missing.
- Independent external benchmark: missing.
- Scalable homology clustering: not yet integrated.
- Separate calibration and threshold-selection sets: missing.

Accordingly, the trained AMP artifact remains offline and the product continues to expose the clearly labelled transparent heuristic.

## Next Priority

Stages 3F and 3G close the protein-process-phenotype literature loop and assemble a deterministic cross-scale evidence chain. The current system is an evidence-chain milestone, not a Graph RAG or full LLM reasoning system. Remaining Stage 3 work is the configurable LLM reasoning layer and measured cellular-state data. In parallel, model deployment still requires defensible AMP negatives and an independent benchmark. Other major gaps are SaProt, controlled structure-enhancement evaluation, and experimental imaging.

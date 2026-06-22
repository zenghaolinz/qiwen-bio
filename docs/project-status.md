# Qiwen Bio Project Status

Updated: 2026-06-22

## Scoring Method

Progress is measured against the 24 tasks in the original four-stage proposal. Complete tasks score 1, partial tasks score 0.5, and not-started tasks score 0. The fourth-stage experimental-image work is part of the full vision but optional to the current core bioinformatics platform.

- Full four-stage vision: **15.5 / 24 = 64.58%**
- Core software scope (stages 1-3): **15.5 / 18 = 86.11%**

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

Score: **3 / 5 = 60%**

| Proposal task | Status | Evidence / gap |
| --- | --- | --- |
| Fetch AlphaFold structure | Complete | Dynamic AlphaFold model/PDB retrieval |
| Extract structural features | Partial | Unified StructureFeatureSummary aggregates pLDDT, CA contacts, mutation neighborhood, and domain overlap; pLDDT band logic is single-sourced. Secondary structure and pockets missing |
| Structure visualization | Complete | Dependency-free interactive backbone/contact views |
| SaProt structure embedding | Not started | Foldseek/SaProt pipeline absent |
| Compare sequence and structure models | Not started | Structure-feature export (`python -m qiwen_bio.structure_cli`) prepares JSONL records for a future controlled benchmark; no benchmark or trained structure-enhanced model yet |

### Structure evidence scope

Structure evidence summary completed; controlled structure-enhanced prediction benchmark not started. The structure layer is an **evidence summary**, not a functional-effect prediction. AlphaFold pLDDT is local model confidence, not pathogenicity or functional-effect confidence. CA contacts within 8 A are geometric proximity, not confirmed biochemical interactions. There is no SaProt, no pocket prediction, no structure-enhanced model, and no LoRA/Adapter in the current branch.

The structure evidence summary now preserves source URLs, distinguishes missing/invalid/unmapped mutation contexts, and avoids functional-effect wording from structure-only evidence. Hypothesis sentences use "candidate site for follow-up functional-impact assessment" language rather than "may affect structure or function", reflecting that domain overlap and predicted-structure context are not functional-effect or pathogenicity evidence.

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

### Mutation-token interface policy

All mutation parsing routes through a single canonical parser,
`qiwen_bio.mutation.parse_mutation()`, which enforces the 20-canonical-amino-
acid alphabet and accepts both `R175H` and the HGVS short protein form
`p.R175H` (normalized to `R175H`). The request-model `Field(pattern=...)`
strings are not duplicated: `models.py` imports the
`MUTATION_TOKEN_PATTERN` constant from `qiwen_bio.mutation`, so the model
boundary and the parser cannot drift apart.

Two interface behaviours coexist **by design**, not by parser inconsistency:

1. **Synthesis / reasoning interfaces** (`ReasoningChainRequest`,
   `ComprehensiveReportRequest`, `AnalysisRequest`,
   `UniProtAnalysisRequest`) accept any string as `mutation` and do not 422 on
   a malformed token. A malformed or wild-type-mismatched mutation is
   surfaced through the evidence chain as an `invalid` / `unavailable`
   mutation step, and downstream impact hypotheses are suppressed. This is
   graceful degradation so a user still receives a usable report.

2. **Direct structure / domain interfaces** (`AlphaFoldAnalysisRequest`,
   `InterProAnnotationRequest`) apply the strict `MUTATION_TOKEN_PATTERN` at
   the request-model boundary and return HTTP 422 for an invalid mutation
   token (e.g. `Z175H`, where `Z` is not an amino acid). These are direct
   single-purpose annotation requests, so rejecting bad input early is more
   appropriate than returning a partial result.

Both behaviours use the same canonical alphabet; the difference is whether
the endpoint degrades gracefully or rejects early. Parser-consistency tests
in `tests/test_api.py` pin both behaviours.

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

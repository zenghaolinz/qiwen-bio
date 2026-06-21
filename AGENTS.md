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
- `qiwen_bio/reporting.py`: Markdown report rendering.
- `qiwen_bio/static/`: dependency-free Web Demo.

## Stage History

| Stage | Status | Delivered | Explicitly not delivered |
| --- | --- | --- | --- |
| 1A Engineering MVP | Complete | Sequence validation, deterministic features, API, Web Demo, evidence chain, Markdown report | ESM embeddings and trained classifier |
| 1B UniProt annotation | Complete | Reviewed entry lookup, sequence/function/GO extraction, provenance, AlphaFold entry discovery | InterPro, KEGG, STRING, PubMed |
| 2A Structure confidence | Complete | Dynamic AlphaFold model URL lookup, PDB parsing, mean/distributed pLDDT, mutation-site confidence | Contact maps, secondary structure, pockets, 3D viewer, SaProt |
| 2B Structure context | Complete | CA contact map, 8 A mutation neighborhood, interactive dependency-free backbone viewer | Secondary structure, pockets, all-atom contacts, SaProt |
| 3 Multiscale reasoning | Planned | Knowledge graph, pathways, literature evidence, phenotype chain | Not started |
| 4 Experimental imaging | Planned | Gel/PCR/microscopy analysis | Not started |

## Next Priority

Implement stage 3A: add STRING interactions and pathway annotations, normalize them into typed graph nodes/edges, and expose a local evidence graph without using an LLM to invent missing links.

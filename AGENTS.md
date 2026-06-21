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
- `qiwen_bio/reporting.py`: Markdown report rendering.
- `qiwen_bio/static/`: dependency-free Web Demo.

## Stage History

| Stage | Status | Delivered | Explicitly not delivered |
| --- | --- | --- | --- |
| 1A Engineering MVP | Complete | Sequence validation, deterministic features, API, Web Demo, evidence chain, Markdown report | ESM embeddings and trained classifier |
| 1B UniProt annotation | Complete | Reviewed entry lookup, sequence/function/GO extraction, provenance, AlphaFold entry discovery | InterPro, KEGG, STRING, PubMed |
| 1C1 Embedding foundation | Complete | Pinned ESM-2 8M provider, 320D mean pooling, versioned disk cache, API and Web trigger | Larger ESM variants, batching, GPU verification |
| 1C2 Trained classifier | Planned | Curated AMP data, homology-aware splits, calibrated classifier and model card | Not started |
| 2A Structure confidence | Complete | Dynamic AlphaFold model URL lookup, PDB parsing, mean/distributed pLDDT, mutation-site confidence | Contact maps, secondary structure, pockets, 3D viewer, SaProt |
| 2B Structure context | Complete | CA contact map, 8 A mutation neighborhood, interactive dependency-free backbone viewer | Secondary structure, pockets, all-atom contacts, SaProt |
| 3A Interaction graph | Complete | STRING interactions, process/pathway enrichment, typed source-linked graph, deterministic Canvas view | Direct KEGG API, PubMed, phenotype reasoning |
| 3B Literature retrieval | Complete | Context-bound PubMed search, structured citations, PMID links, Markdown report section | Abstract/full-text appraisal, claim-level support classification |
| 3C Evidence synthesis | Complete | One server-generated report, six-layer coverage score, optional-service degradation, unified Web result | Claim correctness scoring, phenotype prediction |
| 4 Experimental imaging | Planned | Gel/PCR/microscopy analysis | Not started |

## Next Priority

Implement stage 1C2: select and document a redistributable AMP dataset, remove duplicates, create sequence-similarity-aware train/validation/test splits, precompute pinned embeddings, and train a calibrated baseline before replacing the demo heuristic.

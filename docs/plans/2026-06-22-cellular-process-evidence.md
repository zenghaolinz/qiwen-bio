# Cellular-Process Evidence Layer Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Normalize protein-associated biological processes across UniProt GO, direct KEGG, and seed-linked STRING enrichment without generating phenotype or causal claims.

**Architecture:** A deterministic, network-free builder merges records by stable GO/KEGG identifier and attaches multiple source-specific support objects. UniProt biological-process terms are database annotations, KEGG pathways are direct memberships, and STRING terms are enrichment statistics included only when an `annotated_to` edge originates from the seed protein. The derived layer reuses existing pathway coverage points rather than double-counting evidence.

**Tech Stack:** Python 3.11+, Pydantic, existing UniProt/KEGG/STRING models, FastAPI comprehensive response, pytest.

---

### Task 1: Normalized Process Builder

**Files:**
- Create: `qiwen_bio/cellular_processes.py`
- Create: `tests/test_cellular_processes.py`

1. Write failing tests for GO filtering, stable-ID merging, source precedence, seed-edge filtering, FDR retention, deterministic ordering, and empty inputs.
2. Verify the missing-module failure.
3. Implement process/support/result models and the pure builder.
4. Run focused tests until green.

### Task 2: Comprehensive Synthesis and Report

**Files:**
- Modify: `qiwen_bio/synthesis.py`
- Modify: `qiwen_bio/reporting.py`
- Modify: `tests/test_synthesis.py`
- Modify: `tests/test_api.py`

1. Add failing tests for cellular-process response content, source supports, PubMed context, Markdown rendering, and no phenotype claims.
2. Build the layer after KEGG/STRING retrieval and use its labels for literature context.
3. Keep the seven-layer 100-point score unchanged and update the pathway detail with normalized process count.
4. Run focused tests until green.

### Task 3: Evidence Graph Extension

**Files:**
- Modify: `qiwen_bio/stringdb.py` only if the normalized layer exposes a graph adapter.
- Modify: `tests/test_cellular_processes.py`

1. Add a failing test for typed process nodes/support edges if an adapter is needed.
2. Implement only the minimal adapter required by the Web/API payload; avoid a new graph database.
3. Verify deterministic IDs and provenance.

### Task 4: Verification and Documentation

**Files:**
- Modify: `AGENTS.md`
- Modify: `README.md`
- Modify: `docs/project-status.md`
- Create: `docs/adr/0012-normalized-cellular-process-evidence.md`

1. Run a real TP53 comprehensive build or source-layer smoke test and record counts by source.
2. Document that process association is not activity, directionality, mechanism, or phenotype.
3. Update original-proposal scoring only if the protein-to-pathway/phenotype task meaningfully advances.
4. Run full verification, commit, and push.

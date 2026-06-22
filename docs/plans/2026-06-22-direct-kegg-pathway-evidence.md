# Direct KEGG Pathway Evidence Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add direct, source-distinguished KEGG pathway records for UniProt proteins without redistributing KEGG content or conflating it with STRING enrichment.

**Architecture:** A fixed-host KEGG REST client maps UniProt accession to KEGG gene IDs, links genes to pathways, and retrieves pathway flat files in batches of at most ten. It normalizes only identifiers, names, descriptions, and classes for per-request API/report use. Synthesis prefers direct KEGG for pathway coverage and PubMed context, but can explicitly fall back to STRING enrichment when KEGG is unavailable.

**Tech Stack:** Python 3.11+, httpx, Pydantic, FastAPI dependency injection, KEGG text formats, pytest.

---

### Task 1: KEGG REST Client

**Files:**
- Create: `qiwen_bio/kegg.py`
- Create: `tests/test_kegg.py`

1. Write failing tests for conversion/link TSV parsing, flat-file continuation parsing, deduplication, ten-entry batching, malformed responses, and 404/not-found behavior.
2. Verify the missing-module failure.
3. Implement normalized models and an injected fixed-host client with typed errors.
4. Run focused tests until green.

### Task 2: API Boundary

**Files:**
- Modify: `qiwen_bio/models.py`
- Modify: `qiwen_bio/api.py`
- Modify: `tests/test_api.py`

1. Write failing tests for `/api/v1/pathways/kegg`, request validation, 404, and 502.
2. Add the request model, dependency, endpoint, and error mapping.
3. Run API tests until green.

### Task 3: Comprehensive Synthesis

**Files:**
- Modify: `qiwen_bio/synthesis.py`
- Modify: `qiwen_bio/reporting.py`
- Modify: `tests/test_synthesis.py`
- Modify: `tests/test_api.py`

1. Add failing tests for direct KEGG report content and direct-vs-STRING coverage detail.
2. Pass KEGG through the comprehensive endpoint as an optional service.
3. Use direct pathway names as primary PubMed context and render a separate KEGG section.
4. Preserve the 100-point seven-layer score and explicitly label STRING-only fallback.

### Task 4: Official Verification and Documentation

**Files:**
- Modify: `AGENTS.md`
- Modify: `README.md`
- Modify: `docs/project-status.md`
- Create: `docs/adr/0011-direct-kegg-pathway-evidence.md`

1. Verify TP53 mapping, pathway count, and p53 pathway metadata against official KEGG REST.
2. Record academic-use, copyright, no-bulk-redistribution, and failure-mode decisions.
3. Update proposal progress based on direct KEGG completion.
4. Run full verification, commit, and push.

# InterPro and Pfam Integration Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add direct, auditable InterPro and Pfam domain evidence with mutation-coordinate overlap to the API and comprehensive report.

**Architecture:** A provider-specific client retrieves cursor-paginated InterPro API results for both `interpro` and `pfam`, normalizes entry metadata and protein fragments, and computes positional overlap as a database-coordinate fact. The API exposes the annotation independently, while synthesis treats it as an optional layer with explicit degradation. Domain evidence changes coverage composition but keeps the total at 100.

**Tech Stack:** Python 3.11+, httpx, Pydantic, FastAPI dependency injection, pytest.

---

### Task 1: Client, Pagination, and Mutation Overlap

**Files:**
- Create: `qiwen_bio/interpro.py`
- Create: `tests/test_interpro.py`

1. Write failing tests for InterPro/Pfam metadata, fragments, score/model fields, mutation overlap, cursor pagination, deduplication, and pagination-loop protection.
2. Run focused tests and confirm the missing-module failure.
3. Implement normalized models, injected transport/client, response validation, and typed service errors.
4. Run focused tests until green.

### Task 2: API Boundary

**Files:**
- Modify: `qiwen_bio/models.py`
- Modify: `qiwen_bio/api.py`
- Modify: `tests/test_api.py`

1. Write failing endpoint tests for successful annotation, request validation, 404, and 502 mapping.
2. Add `InterProAnnotationRequest`, dependency injection, and `/api/v1/domains/interpro`.
3. Run API tests until green.

### Task 3: Comprehensive Synthesis

**Files:**
- Modify: `qiwen_bio/synthesis.py`
- Modify: `qiwen_bio/reporting.py`
- Modify: `tests/test_synthesis.py`
- Modify: `tests/test_api.py`

1. Add failing tests for a domain layer, mutation overlap text, and optional-service degradation.
2. Add `domains` to the comprehensive model and pass the client through the API.
3. Rebalance coverage to sequence 10, annotation 15, domains 10, structure 20, interactions 15, pathways 15, literature 15.
4. Render domain coordinates and a strict interpretation boundary.
5. Run focused and full tests.

### Task 4: Official Endpoint and Documentation

**Files:**
- Modify: `AGENTS.md`
- Modify: `README.md`
- Modify: `docs/project-status.md`
- Create: `docs/adr/0010-direct-interpro-pfam-domain-evidence.md`

1. Verify TP53 against official InterPro and Pfam endpoints and record actual entry/overlap counts.
2. Update stage history and original-proposal progress scoring.
3. Run full tests, compilation, live endpoint smoke test, and Git checks.
4. Commit and push `main`.

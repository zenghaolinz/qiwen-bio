# UniProt Mature-Peptide Curation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace ambiguous precursor-level AMP positives with an auditable inventory of exact UniProt mature-peptide annotations and quantify what remains unsuitable for external validation.

**Architecture:** Read positive accessions from the frozen baseline CSV, retrieve each reviewed UniProt JSON entry through an injected HTTP client, and extract only canonical `Peptide` features with exact, in-range coordinates. Preserve every exclusion reason, deduplicate mature sequences, and write ignored sequence rows plus a committed provenance/readiness manifest. This stage does not invent negative labels or call same-source records an external benchmark.

**Tech Stack:** Python 3.11+, httpx, Pydantic, CSV/JSON, pytest.

---

### Task 1: Exact Feature Parser

**Files:**
- Create: `qiwen_bio/mature_peptides.py`
- Create: `tests/test_mature_peptides.py`

1. Write failing tests for exact 1-based coordinate slicing, multiple peptides, fuzzy-coordinate exclusion, out-of-range exclusion, and HTTP injection.
2. Run the focused test and confirm the module/import failure.
3. Implement Pydantic result models, parser, and UniProt client.
4. Run focused tests until green.

### Task 2: Batch Curation and Audit

**Files:**
- Modify: `qiwen_bio/mature_peptides.py`
- Modify: `tests/test_mature_peptides.py`
- Create: `qiwen_bio/mature_peptide_cli.py`

1. Write failing tests for positive-only accession selection, duplicate merging, missing-feature audit, atomic CSV/manifest output, and readiness gates.
2. Implement batch curation and CLI defaults against the frozen baseline.
3. Keep sequence output under ignored `data/datasets/`; commit only the aggregate manifest.

### Task 3: Official Endpoint Run

**Files:**
- Create: `data/manifests/uniprot_mature_amp_audit.json`

1. Fetch all 50 positive accessions from `https://rest.uniprot.org/uniprotkb/{accession}.json`.
2. Verify the official endpoint, exact feature coordinates, counts, duplicate handling, and manifest hash.
3. Record missing/fuzzy/error counts without silently dropping entries.

### Task 4: Documentation and Project Progress

**Files:**
- Modify: `AGENTS.md`
- Modify: `README.md`
- Create: `docs/adr/0009-exact-uniprot-mature-peptide-curation.md`

1. Mark mature-positive curation complete separately from defensible negatives and independent external evaluation.
2. Report the actual coverage and readiness decision.
3. Compare all original proposal stages with delivered artifacts and calculate transparent completion counts.
4. Run full tests, compilation, live artifact validation, and Git checks; commit and push.

"""Reproducible fixture/snapshot test for the TP53 R2H reasoning chain (Task 4).

This test pins the *structure* of the comprehensive-analysis output for the
TP53 R2H case (the stub-data analogue of R175H). It guards against silent
regressions in chain step layout, hypothesis gating, missing-layer tracking,
and report section presence.

The fixture at ``tests/fixtures/tp53_r2h_reasoning_chain.json`` is generated
deterministically from the stub clients (no network). If a deliberate
behaviour change alters the structure, regenerate the fixture with the same
stub inputs and commit the updated file.

Note: the fixture stores step *counts* and flags, not full text, so wording
improvements do not require a fixture update — only structural changes do.
"""

import json
from pathlib import Path

from qiwen_bio.pipeline import AnalysisPipeline
from qiwen_bio.synthesis import build_comprehensive_analysis
from tests.test_synthesis import (
    StubAlphaFoldClient,
    StubInterProClient,
    StubKeggClient,
    StubPubMedClient,
    StubStringClient,
    StubUniProtClient,
)

FIXTURE_PATH = Path(__file__).parent / "fixtures" / "tp53_r2h_reasoning_chain.json"


def _build_result():
    return build_comprehensive_analysis(
        identifier="TP53",
        organism_id=9606,
        mutation="R2H",
        pipeline=AnalysisPipeline(),
        uniprot_client=StubUniProtClient(),
        alphafold_client=StubAlphaFoldClient(),
        string_client=StubStringClient(),
        pubmed_client=StubPubMedClient(),
        interpro_client=StubInterProClient(),
        kegg_client=StubKeggClient(),
    )


def _actual_structure(result) -> dict:
    chain = result.reasoning_chain
    return {
        "chain_type": chain.chain_type,
        "gene": chain.gene,
        "mutation": chain.mutation,
        "missing_layers": chain.missing_layers,
        "summary": chain.summary,
        "steps": [
            {
                "step_id": s.step_id,
                "title": s.title,
                "available": s.available,
                "confidence": s.confidence,
                "evidence_facts_count": len(s.evidence_facts),
                "has_hypothesis": s.hypothesis is not None,
                "evidence_sources_count": len(s.evidence_sources),
            }
            for s in chain.steps
        ],
    }


def test_tp53_r2h_reasoning_chain_matches_snapshot() -> None:
    fixture = json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))
    result = _build_result()

    assert fixture["input"] == {"identifier": "TP53", "organism_id": 9606, "mutation": "R2H"}
    assert fixture["chain"] == _actual_structure(result)
    assert fixture["warnings"] == result.warnings
    assert fixture["report_section_present"] == {
        "coverage": "## Evidence coverage" in result.report_markdown,
        "structure": "## Structure evidence" in result.report_markdown,
        "domains": "## Domain evidence" in result.report_markdown,
        "kegg": "## Direct KEGG pathway evidence" in result.report_markdown,
        "cellular_processes": "## Cellular-process evidence" in result.report_markdown,
        "graph": "## Interaction and pathway evidence" in result.report_markdown,
        "literature": "## PubMed literature retrieval" in result.report_markdown,
        "phenotype_literature": "## Phenotype literature evidence" in result.report_markdown,
        "reasoning_chain": "## Reasoning chain" in result.report_markdown,
        "interpretation_boundary": "## Interpretation boundary" in result.report_markdown,
    }


def test_tp53_r2h_fixture_preserves_provenance_and_uncertainty() -> None:
    """The snapshot must retain evidence sources and uncertainty on every
    available step, and the report must carry the interpretation boundary."""
    result = _build_result()
    chain = result.reasoning_chain
    for step in chain.steps:
        if step.available:
            assert step.evidence_sources, f"step {step.step_id} lost provenance"
            assert step.uncertainty, f"step {step.step_id} lost uncertainty"
    # The interpretation boundary is the global non-causality disclaimer.
    assert "does not establish" in chain.boundary.lower()
    assert "## Interpretation boundary" in result.report_markdown

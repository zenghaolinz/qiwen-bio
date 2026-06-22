"""Comprehensive-report end-to-end tests for the structure evidence summary
(Stage 2 refinement, task 4).

These drive build_comprehensive_analysis with the stub clients and assert the
structure summary section renders correctly across: normal mutation, no
mutation, invalid mutation, no AlphaFold structure, low pLDDT, domain overlap,
and source_url preservation.
"""

from qiwen_bio.pipeline import AnalysisPipeline
from qiwen_bio.synthesis import build_comprehensive_analysis
from tests.test_synthesis import (
    FailingAlphaFoldClient,
    StubAlphaFoldClient,
    StubInterProClient,
    StubKeggClient,
    StubPubMedClient,
    StubStringClient,
    StubUniProtClient,
)


def _build(mutation, alphafold_client=None):
    return build_comprehensive_analysis(
        identifier="TP53",
        organism_id=9606,
        mutation=mutation,
        pipeline=AnalysisPipeline(),
        uniprot_client=StubUniProtClient(),
        alphafold_client=alphafold_client or StubAlphaFoldClient(),
        string_client=StubStringClient(),
        pubmed_client=StubPubMedClient(),
        interpro_client=StubInterProClient(),
        kegg_client=StubKeggClient(),
    )


def test_e2e_normal_mutation_renders_structure_summary_with_source_url() -> None:
    result = _build("R2H")
    chain = result.reasoning_chain
    assert chain is not None
    assert chain.structure_summary is not None
    assert chain.structure_summary.has_structure is True
    assert chain.structure_summary.source_url == "https://example.test/model.pdb"
    assert "## Structure evidence summary" in result.report_markdown
    assert "AlphaFold pLDDT is local model confidence" in result.report_markdown
    assert "https://example.test/model.pdb" in result.report_markdown
    # Mapped mutation: pLDDT value rendered, not "no mutation supplied".
    assert "Mutation site pLDDT: 68.00" in result.report_markdown
    assert "no mutation supplied" not in result.report_markdown


def test_e2e_no_mutation_protein_function_chain() -> None:
    result = _build(None)
    chain = result.reasoning_chain
    assert chain.chain_type == "protein_function"
    assert chain.structure_summary is not None
    assert chain.structure_summary.mutation_parse_status == "not_supplied"
    assert "no mutation supplied" in result.report_markdown
    assert "could not be parsed" not in result.report_markdown


def test_e2e_invalid_mutation_says_could_not_be_parsed_not_no_mutation() -> None:
    result = _build("bad")
    chain = result.reasoning_chain
    assert chain.chain_type == "mutation_impact"
    summary = chain.structure_summary
    assert summary.mutation_parse_status == "invalid"
    assert "could not be parsed" in result.report_markdown
    # Must NOT say "no mutation supplied" for the mutation-site line.
    structure_summary_idx = result.report_markdown.find("## Structure evidence summary")
    summary_section = result.report_markdown[structure_summary_idx:]
    assert "no mutation supplied" not in summary_section


def test_e2e_no_alphafold_structure_marks_layer_missing() -> None:
    result = _build("R2H", alphafold_client=FailingAlphaFoldClient())
    chain = result.reasoning_chain
    assert chain.structure_summary is not None
    assert chain.structure_summary.has_structure is False
    assert "structure" in chain.missing_layers
    assert "No AlphaFold structure is available" in result.report_markdown
    assert chain.structure_summary.source_url is None


def test_e2e_domain_overlap_reported_without_functional_impact_wording() -> None:
    result = _build("R2H")
    md = result.report_markdown
    assert "PF00870" in md
    assert "Domain overlap: yes" in md
    # No forbidden functional-impact phrasing from structure-only evidence.
    assert "可能影响局部结构或功能" not in md
    assert "可能影响蛋白功能" not in md
    assert "disrupts" not in md.lower()


def test_e2e_low_plddt_mutation_site_flags_limited_interpretation() -> None:
    # M1K: position 1 has pLDDT 42 (very_low, < 50). Wild-type M matches.
    result = _build("M1K")
    chain = result.reasoning_chain
    summary = chain.structure_summary
    assert summary.has_structure is True
    assert summary.mutation_site_plddt == 42.0
    assert summary.mutation_site_confidence_band == "very_low"
    assert summary.low_confidence_region is True
    assert "low-confidence predicted region" in result.report_markdown.lower()


def test_e2e_source_url_preserved_in_structure_step_evidence_sources() -> None:
    result = _build("R2H")
    structure_step = next(s for s in result.reasoning_chain.steps if s.step_id == "structure")
    assert "https://example.test/model.pdb" in structure_step.evidence_sources

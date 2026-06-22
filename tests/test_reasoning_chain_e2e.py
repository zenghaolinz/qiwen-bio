"""End-to-end reasoning-chain tests (Task 3 of the quality-convergence pass).

These tests drive the full comprehensive-analysis path with the stub clients
and assert the biological-correctness guardrails from the convergence guide:
no causal language, hypothesis gating on evidence, wild-type mismatch
suppression, and correct chain-type triggering.
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


FORBIDDEN_CAUSAL_TERMS = (
    "导致",
    "证明",
    "必然",
    "确定造成",
    "直接说明",
    "proves",
    "definitive",
    "establishes causality",
)

# "causes" is only forbidden when asserted positively, not when negated in a
# disclaimer (e.g. "not evidence that ... causes a phenotype"). We check it
# per-sentence and allow it inside negated clauses.
FORBIDDEN_POSITIVE_TERMS = ("causes",)


def _build(mutation):
    return build_comprehensive_analysis(
        identifier="TP53",
        organism_id=9606,
        mutation=mutation,
        pipeline=AnalysisPipeline(),
        uniprot_client=StubUniProtClient(),
        alphafold_client=StubAlphaFoldClient(),
        string_client=StubStringClient(),
        pubmed_client=StubPubMedClient(),
        interpro_client=StubInterProClient(),
        kegg_client=StubKeggClient(),
    )


def _assert_no_forbidden_causal_terms(markdown: str) -> None:
    # Absolutely forbidden terms never appear, in any context.
    lowered = markdown.lower()
    for term in FORBIDDEN_CAUSAL_TERMS:
        assert term.lower() not in lowered, (
            f"Report markdown contains forbidden causal term {term!r}"
        )
    # Positive causal verbs are forbidden unless the sentence negates them.
    negation_markers = ("not ", "does not ", "do not ", "cannot ", "no ")
    for sentence in markdown.split("."):
        s = sentence.lower()
        for term in FORBIDDEN_POSITIVE_TERMS:
            if term in s:
                idx = s.find(term)
                prefix = s[:idx]
                assert any(marker in prefix for marker in negation_markers), (
                    f"Report contains non-negated causal term {term!r}: "
                    f"{sentence.strip()[:120]}"
                )


def test_end_to_end_tp53_r175h_mutation_impact_chain() -> None:
    """Scenario 1: TP53 + R2H (stub's R175H analogue) mutation chain.

    The stub UniProt sequence is 'MEEPQSDPSV' (pos 2 = E), and the stub PDB
    has ARG at pos 2, so R2H matches AlphaFold but mismatches UniProt. This
    is a faithful test of the wild-type-validation guardrail.
    """
    result = _build("R2H")
    chain = result.reasoning_chain

    assert chain is not None
    assert chain.chain_type == "mutation_impact"
    assert chain.mutation == "R2H"
    assert [s.step_id for s in chain.steps] == [
        "mutation", "structure", "function", "pathway", "phenotype",
    ]
    # Mutation step present and records the mismatch (UniProt has E, not R).
    mutation_step = chain.steps[0]
    assert mutation_step.available is True
    assert any("mismatch" in f.lower() or "does not match" in f.lower()
               for f in mutation_step.evidence_facts)
    # Wild-type mismatch suppresses downstream impact hypotheses.
    assert mutation_step.hypothesis is None
    function_step = next(s for s in chain.steps if s.step_id == "function")
    assert function_step.hypothesis is None
    pathway_step = next(s for s in chain.steps if s.step_id == "pathway")
    assert pathway_step.hypothesis is None
    # Phenotype step carries the ADR-0013-gated literature hypothesis.
    phenotype_step = next(s for s in chain.steps if s.step_id == "phenotype")
    assert phenotype_step.hypothesis is not None
    assert "abstract-level hypothesis" in phenotype_step.hypothesis
    # Every step has at least one evidence source (provenance).
    for step in chain.steps:
        assert step.evidence_sources or not step.available, (
            f"Step {step.step_id} has no provenance"
        )
    # No forbidden causal language anywhere in the report.
    _assert_no_forbidden_causal_terms(result.report_markdown)
    assert "## Reasoning chain" in result.report_markdown


def test_end_to_end_no_mutation_protein_function_chain() -> None:
    """Scenario 2: TP53 with no mutation -> protein_function chain."""
    result = _build(None)
    chain = result.reasoning_chain

    assert chain is not None
    assert chain.chain_type == "protein_function"
    assert chain.mutation is None
    assert "mutation" not in {s.step_id for s in chain.steps}
    assert [s.step_id for s in chain.steps] == [
        "structure", "function", "pathway", "phenotype",
    ]
    # No mutation -> no mutation-specific impact hypotheses.
    for step in chain.steps:
        if step.step_id in ("function", "pathway"):
            assert step.hypothesis is None, (
                f"{step.step_id} step must not carry an impact hypothesis "
                f"without a mutation"
            )
    _assert_no_forbidden_causal_terms(result.report_markdown)


def test_end_to_end_wild_type_mismatch_suppresses_impact_hypotheses() -> None:
    """Scenario 3: a mutation whose wild type does not match the sequence.

    R2H against UniProt sequence 'MEEPQSDPSV' (pos 2 = E, not R) must:
    - record the mismatch as a fact;
    - suppress function/pathway impact hypotheses;
    - still return annotation facts (not an empty chain);
    - not emit forbidden causal language.
    """
    result = _build("R2H")
    chain = result.reasoning_chain

    mutation_step = chain.steps[0]
    assert any("mismatch" in f.lower() or "does not match" in f.lower()
               for f in mutation_step.evidence_facts)
    # Suppression: no functional-impact or pathway-impact hypothesis.
    function_step = next(s for s in chain.steps if s.step_id == "function")
    pathway_step = next(s for s in chain.steps if s.step_id == "pathway")
    assert function_step.hypothesis is None
    assert pathway_step.hypothesis is None
    # Facts still present (annotation is not destroyed by the mismatch).
    assert len(function_step.evidence_facts) > 0
    assert len(pathway_step.evidence_facts) > 0
    _assert_no_forbidden_causal_terms(result.report_markdown)


def test_end_to_end_malformed_mutation_keeps_chain_type_but_marks_unavailable() -> None:
    """A malformed mutation still produces a mutation_impact chain (user
    intent) but the mutation step is unavailable with an error fact."""
    result = _build("NOTAMUTATION")
    chain = result.reasoning_chain

    assert chain.chain_type == "mutation_impact"
    mutation_step = chain.steps[0]
    assert mutation_step.available is False
    assert mutation_step.confidence == "insufficient"
    assert any("invalid" in f.lower() or "malformed" in f.lower()
               for f in mutation_step.evidence_facts)
    # No impact hypotheses downstream of an unparseable mutation.
    for step in chain.steps[1:]:
        if step.step_id in ("function", "pathway"):
            assert step.hypothesis is None
    _assert_no_forbidden_causal_terms(result.report_markdown)


def test_end_to_end_missing_layers_are_explicit_not_fabricated() -> None:
    """When structure is unavailable, the chain must list it in missing_layers
    and mark the step unavailable rather than fabricate structural evidence."""
    result = build_comprehensive_analysis(
        identifier="TP53",
        organism_id=9606,
        mutation="R2H",
        pipeline=AnalysisPipeline(),
        uniprot_client=StubUniProtClient(),
        alphafold_client=FailingAlphaFoldClient(),
        string_client=StubStringClient(),
        pubmed_client=StubPubMedClient(),
        interpro_client=StubInterProClient(),
        kegg_client=StubKeggClient(),
    )
    chain = result.reasoning_chain
    assert "structure" in chain.missing_layers
    structure_step = next(s for s in chain.steps if s.step_id == "structure")
    assert structure_step.available is False
    assert structure_step.confidence == "insufficient"
    _assert_no_forbidden_causal_terms(result.report_markdown)

# ADR-0014: Cross-scale Reasoning Chain with Layered Triggering and Gated Hypotheses

## Status

Accepted

## Context

The project plan's core innovation ("创新点1: 从单点预测到多尺度推理") requires connecting a mutation (or, in its absence, the protein itself) to structure, function, pathway, and phenotype evidence in a single chain. Stage 3F delivered only the phenotype-literature tail of that chain. Stages 2A-2C, 3A-3E already retrieve the upstream layers (AlphaFold structure, InterPro/Pfam domains, UniProt function/GO, direct KEGG, STRING, cellular-process normalization, phenotype literature), but no object joins them into a cross-scale chain.

The chain must respect the biological-correctness guardrails established by ADR-0012 (cellular-process layer never infers activity/direction/mechanism/causality/phenotype) and ADR-0013 (phenotype hypotheses gated on abstract-level `supports` evidence). It must also close a loose mutation-parsing gap: `synthesis.py` parsed the mutation position with `mutation[1:-1].isdigit()`, which does not validate the wild-type residue against the UniProt sequence.

## Decision

Add a synthesis-only module `qiwen_bio/reasoning_chain.py` that assembles the already-retrieved layers into a `ReasoningChain` with no new network requests and no new biological claims.

### Layered triggering

- **Mutation supplied** (any non-empty string, even malformed) → `mutation_impact` chain with five steps: mutation → structure → function → pathway → phenotype.
- **No mutation** → `protein_function` chain with four steps (no mutation step): structure → function → pathway → phenotype.
- A malformed mutation produces a `mutation_impact` chain whose mutation step is `available=False` with an error fact, rather than silently switching chain types.

### Wild-type validation

The mutation step parses the mutation with a strict `^[A-Z]\d+[A-Z]$` regex and validates the wild-type residue against `annotation.sequence[position-1]`. A mismatch is recorded as a fact and **suppresses downstream "可能影响" hypotheses**: facts remain, but no functional-impact inference is drawn from a mutation that does not match the annotated sequence. This is stricter than the legacy `mutation[1:-1].isdigit()` slice and is the canonical parse for the chain.

### Hypothesis gating

Each step states deterministic, source-bound facts first. A hypothesis sentence ("可能影响/可能相关/提示", explicitly prefixed "假设：") is emitted only when explicit upstream evidence supports it:

- mutation step hypothesis requires: wild-type matches AND a domain overlap AND a structure mutation site;
- function step hypothesis requires: a mutation step hypothesis AND a domain overlap;
- pathway step hypothesis requires: a function step hypothesis AND a directly supported process;
- phenotype step hypothesis reuses the ADR-0013-gated `phenotype_literature.hypotheses[0]`.

Every hypothesis carries an uncertainty clause stating it is not a causal, directional, or phenotype conclusion and requires experimental validation.

### Missing layers

Absent layers (structure, kegg, string_graph, phenotype_literature) are recorded in `missing_layers` and the corresponding step is marked `available=False` with `confidence="insufficient"`. The chain is never empty: it always returns a structured object so the caller can see what was and was not available.

## Evidence

The unit test suite (`tests/test_reasoning_chain.py`) covers: mutation chain with matching wild-type (hypotheses fire at each gated step), mutation chain with wild-type mismatch (downstream hypotheses suppressed, facts only), no-mutation protein-function chain (no "可能影响" hypotheses), missing layers (step marked unavailable), malformed mutation (step unavailable, chain type preserved), and no-supports-articles phenotype step.

A pre-existing indentation bug in `reporting.py` (the `process_lines.append` call was outside the `for` loop, so only the last cellular process rendered) was fixed as part of this stage with a regression test, because the chain depends on correct process rendering.

## Consequences

### Positive

- The platform now delivers the plan's central cross-scale chain: mutation → structure → function → pathway → phenotype.
- Hypotheses are explicitly labelled and gated on upstream evidence; no unsupported inference is presented as biological fact.
- Wild-type validation closes the loose mutation-parsing gap and prevents impact hypotheses from a mutation that does not match the annotated sequence.
- The chain is a pure synthesis layer with no new network calls, so it adds negligible latency and cannot introduce new service failures.
- Missing layers are explicit, not silent.

### Negative

- The chain is deterministic and rule-based; it cannot perform the nuanced multi-step reasoning an LLM could (deferred per the project's LLM-hallucination guardrail).
- The "可能影响" hypotheses are coarse: a domain overlap is a coordinate containment observation, not a functional-effect prediction, and the hypothesis sentence makes this explicit but the gate is permissive.
- Wild-type validation against the UniProt sequence may flag legitimate isoform numbering differences as mismatches; the fact records this possibility rather than silently suppressing.
- The chain depends on all upstream layers being assembled first, so it is only available in the comprehensive and dedicated reasoning-chain endpoints.

### Neutral

- The chain introduces a new `ReasoningChain` model on `ComprehensiveAnalysis` and a new `/api/v1/reasoning/chain` endpoint, but does not change coverage scoring (the chain synthesizes already-scored layers).

## Alternatives Considered

- Use an LLM to generate the chain reasoning: rejected for this stage per the project's LLM-hallucination guardrail (AGENTS.md: "Do not present unsupported model output as biological fact"). A configurable LLM layer remains a later, separately-gated decision.
- Assert causal language ("导致功能丧失"): rejected; violates ADR-0012/0013 and AGENTS.md.
- Populate `CellularProcessEvidence.phenotype_hypotheses` from the chain: rejected; ADR-0012 keeps that field empty. Chain hypotheses live on the `ReasoningChain` and the separate `PhenotypeLiteratureEvidence`.
- Require a mutation for any chain: rejected; the plan's "功能三：蛋白到通路的证据链" does not require a mutation, so a no-mutation protein-function chain is provided.

## References

- ADR-0012: Normalize cellular processes without generating phenotype claims
- ADR-0013: Phenotype literature with conservative abstract-metadata claim-support classification
- Project plan: 创新点1 (从单点预测到多尺度推理), 功能二 (突变影响推理), 功能三 (蛋白到通路的证据链)

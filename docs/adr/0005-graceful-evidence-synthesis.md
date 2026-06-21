# ADR-0005: Use Graceful Degradation and Evidence Coverage Scoring

## Status

Accepted

## Context

A comprehensive analysis depends on UniProt, AlphaFold DB, STRING, and PubMed. Treating every upstream as mandatory would make the report fragile, while silently omitting failed layers would make it misleading. A single opaque confidence number could also be mistaken for biological correctness or pathogenicity.

## Decision

Require UniProt resolution as the identity and sequence anchor. Treat AlphaFold, STRING, and PubMed as optional evidence layers: catch their expected service/not-found errors, return available results, and expose warnings for missing layers. Generate one server-side report from the same structured response used by the Web UI.

Report an evidence coverage score with fixed, inspectable weights: sequence 15, UniProt annotation 20, structure 20, interactions 15, process/pathway terms 15, and literature retrieval 15. Award a component only when that layer contains usable data. Label totals below 50 as limited, 50-79 as partial, and 80-100 as comprehensive. Always state that coverage is not confidence, correctness, pathogenicity, or functional-effect probability.

## Consequences

### Positive

- One upstream outage no longer destroys all useful analysis.
- Every score is explainable from six visible components.
- The API, Web UI, and Markdown report share one structured result.
- Missing evidence is explicit rather than silently omitted.

### Negative

- The endpoint can take several seconds because it performs multiple live requests.
- Fixed weights reflect product completeness priorities, not learned evidence quality.
- A 100/100 report may still contain irrelevant literature or non-causal associations.

### Neutral

- UniProt remains mandatory because all downstream identifiers and sequence coordinates depend on it.

## Alternatives Considered

- Fail the entire request when any service fails: rejected as too fragile for exploratory research.
- Average upstream confidence metrics: rejected because pLDDT, STRING score, FDR, and literature retrieval are not commensurate.
- Use an LLM to assign confidence: rejected because it would be opaque and difficult to reproduce.

## References

- `docs/adr/0001-alphafold-on-demand-analysis.md`
- `docs/adr/0003-typed-in-memory-evidence-graph.md`
- `docs/adr/0004-pubmed-metadata-retrieval.md`

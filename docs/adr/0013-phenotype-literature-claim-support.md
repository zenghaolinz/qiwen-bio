# ADR-0013: Phenotype Literature with Conservative Abstract-Metadata Claim-Support Classification

## Status

Accepted

## Context

Stage 3E (ADR-0012) normalized cellular-process evidence across UniProt GO, direct KEGG, and seed-linked STRING enrichment, and deliberately kept `phenotype_hypotheses` empty because database association does not establish activity, direction, mechanism, causality, or phenotype. The project's stage 3F goal is to connect directly supported cellular processes to phenotype-oriented literature and to emit hypotheses only when each link has explicit evidence and uncertainty.

Stage 3B (ADR-0004) retrieves PubMed metadata only and explicitly deferred abstract-based claim classification until a claim-level validation design existed. Stage 3F supplies that design, but it must not overstate evidence: search retrieval is never support by itself, and abstract text is not full-text appraisal.

## Decision

Add a phenotype-literature layer that operates strictly at the abstract-metadata boundary:

- Build a `PhenotypeLiteratureEvidence` object separate from `CellularProcessEvidence`, so ADR-0012's empty `phenotype_hypotheses` invariant on the cellular-process layer is preserved. Hypotheses live only on the new object.
- Restrict hypothesis eligibility to processes with **direct** support (`database_annotation` from UniProt GO or `database_membership` from direct KEGG). Enrichment-only STRING processes are literature-context only.
- For each eligible process, run a phenotype-oriented PubMed query `("<GENE>"[Title/Abstract]) AND ("<process label>"[Title/Abstract])`, then fetch abstracts via efetch (`retmode=xml`) and parse structured `AbstractText` sections including their `Label` attribute.
- Classify each article as `supports`, `mentions`, or `no_abstract`. `supports` requires all of: the gene token present (word-boundary match so `TP53` does not match `TP53BP1`), at least one substantive process-label token present, and a pinned conservative outcome verb (e.g. `regulates`, `associated with`, `activates`, `suppresses`, `required for`) in the title or abstract.
- Emit a phenotype hypothesis only when a directly supported process has at least one `supports` article. Each hypothesis cites the direct support sources, the PMID(s), the matched verb(s), and an explicit uncertainty statement.
- Degrade gracefully on PubMed failure: skip the process link, emit no hypothesis, and let the main literature layer surface its own warning.

## Evidence

Live TP53 verification on 2026-06-22: the `p53 signaling pathway` query returned PMIDs 36859359 and 31365877. Both abstracts co-mention TP53 and the process label and use outcome verbs (`activate`, `suppress`, `regulates`, `associated with`, `correlation`), so both were classified `supports` and hypotheses were emitted. The efetch XML path was verified against the official NCBI endpoint.

## Consequences

### Positive

- Phenotype hypotheses are gated on explicit, auditable abstract-level evidence rather than database association alone.
- The cellular-process layer's empty-hypothesis invariant (ADR-0012) is preserved; provenance stays clean.
- Every hypothesis carries its uncertainty and the matched language, so a reader can trace why an article was labelled `supports`.
- Optional-layer degradation keeps the comprehensive report usable when PubMed is unavailable.

### Negative

- Abstract-metadata classification is not full-text appraisal and cannot detect negation, hedging, or contradictory results within an article.
- Substring verb matching is deliberately permissive (e.g. `deregulated` contains `regulate`), so the `supports` label is a coarse filter, not a graded evidence assessment.
- The pinned verb list is a behaviour-changing constant; extending it broadens what counts as `supports`.
- Additional efetch requests add latency and depend on NCBI availability.

### Neutral

- The classifier is deterministic and rule-based; no LLM is used for claim classification in this stage.

## Alternatives Considered

- Populate `CellularProcessEvidence.phenotype_hypotheses` directly: rejected because it violates ADR-0012 and conflates database association with literature support.
- Use an LLM to classify claim support: rejected for this stage because it introduces non-determinism and hallucination risk; the project guards against presenting unsupported model output as biological fact.
- Include enrichment-only STRING processes as hypothesis-eligible: rejected because neighbor-only enrichment must not become seed evidence (ADR-0012).
- Require full-text appraisal before any hypothesis: rejected as out of scope; full-text appraisal remains a later stage. The abstract-metadata boundary is the conservative floor.

## References

- ADR-0004: PubMed metadata retrieval
- ADR-0012: Normalize cellular processes without generating phenotype claims
- https://www.ncbi.nlm.nih.gov/books/NBK25501/
- https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi

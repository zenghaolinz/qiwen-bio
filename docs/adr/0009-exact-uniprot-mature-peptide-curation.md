# ADR-0009: Curate Only Exact UniProt Mature-Peptide Features

## Status

Accepted

## Context

The first AMP baseline labels complete UniProt records by keyword. Many antimicrobial records are precursor proteins, while the intended task is often mature-peptide classification. Training or evaluating on full precursors can teach signal peptides and processing context rather than mature antimicrobial sequence properties.

## Decision

Retrieve each frozen positive accession from the official UniProt JSON endpoint and accept only features with `type=Peptide`, exact start and end modifiers, valid 1-based closed coordinates, and canonical amino acids. Slice the mature sequence from the returned parent sequence, retain accession, description, coordinates, feature index, and parent-sequence digest, and audit every exclusion.

Deduplicate exact mature sequences across accessions. Keep sequence rows outside Git and commit an aggregate manifest bound to the frozen baseline digest. Keep three independent readiness gates: exact mature positives, defensible negatives, and an independent external benchmark. Product replacement requires all three; this stage only satisfies the first.

## Evidence

On 2026-06-22, all 50 frozen positive accessions returned successfully. Twenty entries (40%) contained exact `Peptide` features. The process extracted 26 instances and 26 unique mature sequences, with lengths from 9 to 56 residues and a mean of 31.12. Thirty entries had no exact peptide feature.

## Consequences

### Positive

- Mature sequences are derived from explicit database coordinates rather than length heuristics.
- Missing annotations remain visible and measurable.
- The audit is reproducible and bound to the original dataset SHA-256.
- Readiness cannot be inferred from positive coverage alone.

### Negative

- Only 40% of the current positives provide exact mature-chain annotations.
- UniProt feature annotation is not an independent external benchmark.
- Feature presence does not by itself establish assay conditions, potency, toxicity, or target organism.
- No defensible negative set is delivered by this decision.

## Alternatives Considered

- Strip signal/propeptide regions heuristically: rejected because processing boundaries are protein-specific.
- Treat the remaining C-terminal sequence as mature: rejected because it can include propeptides or unrelated chains.
- Download a third-party AMP database immediately: deferred until redistribution terms, versioning, label semantics, and negative construction are qualified.
- Call mature positives an external test set: rejected because they originate from the same accessions as the training baseline.

## References

- https://rest.uniprot.org/uniprotkb/P82270.json
- https://www.uniprot.org/help/license

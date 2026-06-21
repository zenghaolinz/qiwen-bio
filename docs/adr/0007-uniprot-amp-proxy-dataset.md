# ADR-0007: Start AMP Training Work with an Auditable UniProt Proxy Dataset

## Status

Accepted

## Context

A trained AMP classifier requires explicit label provenance and homology-aware evaluation. Random row splits can place close homologs in training and test sets, inflating metrics. A convenient query for sequences lacking an AMP annotation does not create experimentally confirmed negatives.

## Decision

Build the first small baseline from reviewed UniProtKB entries of length 5 to 200. Positives carry antimicrobial keyword `KW-0929`. Proxy negatives are reviewed entries without that keyword, and every manifest states that missing annotation is not proof of missing antimicrobial activity.

Normalize canonical amino-acid sequences, merge same-label duplicates, and exclude conflicting labels. Compute deterministic Needleman-Wunsch global identity and connect records at or above 80% identity. Assign entire connected components to a seeded 70/15/15 split. Keep sequence rows outside Git and commit a manifest containing queries, policies, license URL, counts, parameters, and a content digest.

## Consequences

### Positive

- Label policy, source, query, and processing parameters are auditable.
- Near-identical sequence clusters cannot cross evaluation splits.
- Conflict and duplicate handling are deterministic and tested.
- The small baseline can be reproduced with one CLI command.

### Negative

- Proxy negatives may contain unannotated AMPs and introduce label noise.
- UniProt keyword coverage and database contents can change over time.
- Pairwise dynamic-programming clustering is O(n^2) and unsuitable for large datasets.
- This baseline alone cannot justify replacing the demo predictor.

## Alternatives Considered

- Treat unannotated entries as confirmed negatives: rejected because annotation absence is not biological evidence.
- Use a random row split: rejected because homolog leakage would invalidate held-out metrics.
- Introduce MMseqs2 immediately: deferred until dataset scale requires an external dependency.
- Commit retrieved sequences: rejected for now; local rows are regenerated while the provenance manifest is versioned.

## References

- https://rest.uniprot.org/uniprotkb/search
- https://www.uniprot.org/help/license
- https://www.uniprot.org/help/keywords

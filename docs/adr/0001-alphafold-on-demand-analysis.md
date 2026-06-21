# ADR-0001: Analyze AlphaFold Structures On Demand

## Status

Accepted

## Context

The MVP needs pLDDT and mutation-site context while remaining lightweight on a local machine. AlphaFold DB file versions change, structures can be several megabytes, and the project does not yet have a persistence layer or background workers.

## Decision

Query the official AlphaFold prediction metadata endpoint for the current versioned PDB URL, download the PDB on demand, and parse one CA record per residue in a pure function. Keep the HTTP client behind FastAPI dependency injection. Do not persist structures in this stage.

## Consequences

### Positive

- No hard-coded AlphaFold model version.
- Low operational complexity and deterministic parser tests.
- Clear boundary for adding caching later.

### Negative

- Repeated requests download the same structure.
- Large proteins increase latency and upstream dependence.
- No offline mode yet.

### Neutral

- pLDDT is reported only as model confidence, not as functional or pathogenic evidence.

## Alternatives Considered

- Hard-code a `model_v4` or `model_v6` URL: rejected because releases change.
- Store all structures immediately: rejected because persistence and cache invalidation are premature.
- Parse mmCIF with a heavy structural biology dependency: deferred until richer structural features require it.

## References

- https://alphafold.ebi.ac.uk/api/prediction/P04637
- https://alphafold.ebi.ac.uk/entry/P04637


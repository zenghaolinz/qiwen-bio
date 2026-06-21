# ADR-0002: Use CA Geometry and a Dependency-Free Structure Viewer

## Status

Accepted

## Context

Stage 2B needs a contact map, mutation neighborhood, and interactive structure view. Adding a large browser molecular-viewer package or a CDN dependency would increase supply-chain and offline reliability risks. Returning every all-atom contact would also make responses excessively large.

## Decision

Parse and return one CA coordinate per standard residue. Define non-local contacts as CA distance at most 8 A with sequence separation greater than 2. Return at most 10,000 contact pairs while retaining the total and a truncation flag. Define the mutation neighborhood as all CA atoms within 8 A of the validated mutation residue. Render the backbone and contact map with native Canvas and pointer events.

## Consequences

### Positive

- No new runtime or CDN dependency.
- Small, inspectable visualization code and deterministic geometry tests.
- Contact and neighborhood definitions are explicit in the API.

### Negative

- CA contacts are a coarse proxy and omit side-chain orientation.
- The viewer lacks molecular surfaces, atoms, ligands, and publication-grade rendering.
- Pairwise contact calculation can become expensive for very large proteins.

### Neutral

- Geometry is presented as proximity only and is not treated as functional evidence.

## Alternatives Considered

- 3Dmol.js or NGL from a CDN: deferred to avoid a new remote dependency.
- Vendor a molecular viewer: rejected for the current MVP because of bundle size and maintenance cost.
- All-atom contact analysis: deferred until functional-site analysis requires it.

## References

- https://alphafold.ebi.ac.uk/entry/P04637
- `docs/adr/0001-alphafold-on-demand-analysis.md`

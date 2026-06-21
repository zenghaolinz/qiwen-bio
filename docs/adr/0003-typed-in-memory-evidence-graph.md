# ADR-0003: Normalize STRING Evidence into an In-Memory Typed Graph

## Status

Accepted

## Context

Stage 3A needs to connect a seed protein to interaction partners and biological processes/pathways while preserving source provenance. The current product analyzes one protein at a time and has no persistence or cross-user graph requirements. Introducing Neo4j or another graph database now would add deployment and schema-migration cost without improving this request-scoped workflow.

## Decision

Call STRING's official `network` and `enrichment` JSON endpoints. Normalize responses into Pydantic `GraphNode`, `GraphEdge`, and `EvidenceGraph` models. Use explicit node types (`protein`, `process`, `pathway`) and edge types (`interacts_with`, `annotated_to`). Preserve interaction scores, enrichment FDR, source names, and source URLs. Deduplicate undirected interactions and cap enrichment terms at the ten most significant supported records. Return the graph in memory and render it with a deterministic Canvas layout.

## Consequences

### Positive

- Every graph relationship remains traceable to STRING.
- Stable schemas support later PubMed and phenotype evidence sources.
- No new service or storage dependency.
- Deterministic layouts make visual regression review practical.

### Negative

- Graphs are rebuilt for every request and are not queryable across analyses.
- STRING availability and rate limits affect this feature.
- The current graph is intentionally small and request scoped.

### Neutral

- Pathway labels from STRING enrichment are not equivalent to a direct KEGG, Reactome, or WikiPathways integration.

## Alternatives Considered

- Neo4j: deferred until persistent multi-analysis graph queries are required.
- SQLite node/edge tables: useful for caching later, but unnecessary before persistence requirements are defined.
- Render STRING's image endpoint: rejected because it would hide normalized evidence and prevent local interaction with typed nodes.

## References

- https://string-db.org/help/api/
- https://string-db.org/api/json/network
- https://string-db.org/api/json/enrichment

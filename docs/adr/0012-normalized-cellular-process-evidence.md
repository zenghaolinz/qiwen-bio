# ADR-0012: Normalize Cellular Processes Without Generating Phenotype Claims

## Status

Accepted

## Context

The platform now has three overlapping process sources: UniProt GO biological-process annotations, direct KEGG pathway memberships, and STRING enrichment terms. Listing them independently creates duplicates and can make correlated evidence look like separate biological conclusions. STRING enrichment may also describe network neighbors rather than the seed protein.

## Decision

Build a deterministic, network-free cellular-process layer from already retrieved records. Merge nodes by stable GO or KEGG identifier while retaining source-specific supports:

- UniProt GO as `database_annotation`;
- KEGG as `database_membership`;
- STRING as `enrichment_statistic` with FDR.

Include a STRING term only when an `annotated_to` edge starts at a protein node whose label matches the graph seed. Prefer directly supported process labels when ordering PubMed context. Reuse the existing pathway coverage points instead of awarding a second score for derived normalized evidence.

Always return an empty `phenotype_hypotheses` list. The layer may organize process evidence but must not infer activity, directionality, mechanism, causality, or phenotype.

## Evidence

Live TP53 verification on 2026-06-22 produced 92 normalized process/pathway records. Supports were UniProt GO 67, KEGG 20, and STRING enrichment 7; overlapping stable IDs were merged. No phenotype hypotheses were generated.

## Consequences

### Positive

- Stable IDs remove duplicate process/pathway rows without hiding provenance.
- Multi-source support remains inspectable at record level.
- Neighbor-only STRING terms cannot become seed evidence.
- Literature queries prioritize direct annotations and memberships.
- The coverage score does not double-count a derived layer.

### Negative

- Database associations still do not measure cellular state or pathway activity.
- Different sources may use related but non-identical identifiers that cannot be merged automatically.
- The layer does not perform GO ontology ancestry or semantic similarity expansion.
- Phenotype inference remains a later, evidence-gated stage.

## Alternatives Considered

- Concatenate all labels: rejected because duplicates and provenance inflation remain.
- Merge by text similarity: rejected because similar labels are not stable biological identity.
- Include all STRING enrichment terms: rejected because some terms are connected only to network neighbors.
- Add new coverage points: rejected because normalized records derive from evidence already scored.
- Generate phenotype hypotheses from pathway names: rejected because association does not support causality.

# ADR-0011: Keep Direct KEGG Membership Separate from STRING Enrichment

## Status

Accepted

## Context

STRING enrichment can return KEGG-labelled terms, but those results are enrichment statistics over a protein set rather than direct KEGG records for the seed protein. Treating both as one source hides provenance and can make an enriched pathway appear to be a direct membership assertion.

KEGG is a copyrighted database product. Its official policy permits academic website/API use, but repository-scale redistribution is outside this project's needs.

## Decision

Use the official fixed-host KEGG REST API in three steps: convert a UniProt accession to KEGG gene identifiers, link those genes to pathways, then retrieve pathway flat-file metadata in batches of at most ten. Limit responses to 1-50 pathways and report both the total linked count and whether returned records are truncated.

Retrieve KEGG data on demand only. Do not bundle, persist, train on, or redistribute KEGG pathway content. Preserve KEGG URLs and copyright-policy URL in every response.

In synthesis, direct KEGG records are the primary pathway layer and PubMed context. If KEGG fails but STRING process/pathway enrichment exists, keep the pathway coverage points but label the detail `STRING enrichment fallback`. If neither exists, mark the layer unavailable.

## Evidence

Live verification on 2026-06-22 mapped TP53 UniProt `P04637` to KEGG gene `hsa:7157`. It returned 51 linked pathways; the default limit returned 20 with `truncated=true`. The result included `hsa04115`, p53 signaling pathway, classified under Cellular Processes / Cell growth and death.

## Consequences

### Positive

- Direct membership and enrichment evidence have distinct provenance.
- KEGG pathway names improve literature query context without relying on STRING labels.
- Response limits and truncation are explicit.
- The comprehensive report still degrades when KEGG is unavailable.

### Negative

- Membership does not indicate pathway activation, effect direction, or causal phenotype.
- KEGG availability and rate behavior are external dependencies.
- Only a limited number of pathways are expanded into metadata per request.
- KEGG content cannot be treated as a redistributable local dataset.

## Alternatives Considered

- Continue using STRING KEGG enrichment only: rejected because it does not provide direct KEGG provenance.
- Scrape KEGG HTML pages: rejected in favor of the documented REST interface.
- Cache or commit pathway records: rejected because on-demand retrieval is sufficient and respects the content boundary.
- Award separate coverage points to both KEGG and STRING: rejected because correlated pathway evidence should not inflate the 100-point coverage score.

## References

- https://www.kegg.jp/kegg/rest/keggapi.html
- https://www.kegg.jp/kegg/legal.html
- https://rest.kegg.jp/conv/genes/uniprot:P04637
- https://rest.kegg.jp/link/pathway/hsa:7157

# ADR-0004: Retrieve PubMed Metadata without Automated Claim Summaries

## Status

Accepted

## Context

Stage 3B needs literature provenance for proteins and graph terms. Titles and metadata can identify candidate reading, but automatically summarizing or classifying support from titles alone would overstate evidence. NCBI E-utilities also imposes usage guidance and rate limits, especially without an API key.

## Decision

Use the official PubMed ESearch JSON endpoint to obtain up to ten PMIDs and ESummary to retrieve bibliographic metadata. Build an explicit query from one protein and up to five sanitized context terms. Preserve the exact query in the response. Return PMID, title, up to ten authors, journal, publication date, DOI, and the official PubMed URL. Generate a Markdown citation section, but always state that retrieval does not validate a biological claim. Make two sequential requests per search and do not retrieve abstracts or full text in this stage.

## Consequences

### Positive

- Queries and citations are reproducible and auditable.
- No API key or paid literature service is required.
- The system does not invent summaries or support labels from metadata.

### Negative

- Relevance depends on PubMed indexing and query wording.
- Metadata alone cannot determine study quality or claim support.
- Live requests add latency and depend on NCBI availability.

### Neutral

- Context terms currently come from the most significant STRING enrichment results.

## Alternatives Considered

- Fetch abstracts and run LLM summaries: deferred until claim-level citation validation is designed.
- Use Europe PMC: potentially useful for open-access metadata, but PubMed is the requested primary literature index for this stage.
- Scrape PubMed HTML: rejected because E-utilities provides a supported structured API.

## References

- https://www.ncbi.nlm.nih.gov/books/NBK25501/
- https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi
- https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi

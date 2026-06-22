# ADR-0010: Retrieve InterPro and Pfam as Direct Coordinate Evidence

## Status

Accepted

## Context

The platform already connects sequence, AlphaFold geometry, interaction networks, pathways, and literature, but it lacks direct domain coordinates. UniProt feature text alone does not provide the combined InterPro/Pfam hierarchy and member-database locations needed to state whether a supplied mutation position lies inside an annotated domain.

## Decision

Query the official InterPro API separately for `interpro` and `pfam` entries by UniProt accession. Follow cursor pagination, cap traversal at 100 pages, reject pagination loops, retain source-specific accessions and URLs, and normalize all valid protein fragments. Merge duplicate entries and locations without collapsing InterPro and Pfam into one source.

Compute mutation overlap only as inclusive coordinate containment. Expose the result through a dedicated FastAPI endpoint and an optional comprehensive-report layer. Rebalance evidence coverage to sequence 10, annotation 15, domains 10, structure 20, interactions 15, pathways 15, and literature 15, preserving the 100-point total.

## Evidence

Live verification on 2026-06-22 for TP53 (`P04637`) returned 13 entries and 13 location fragments: nine InterPro and four Pfam. Position 175 overlapped five entries, including Pfam `PF00870` at residues 100-288.

## Consequences

### Positive

- Domain coordinates and hierarchy come from direct database records.
- InterPro and Pfam provenance remain distinguishable.
- Cursor pagination and malformed responses fail visibly.
- Pagination URLs are restricted to HTTPS on `www.ebi.ac.uk/interpro/api/`, preventing server-directed SSRF.
- Optional-service failure does not prevent the rest of the comprehensive report.

### Negative

- Broad families and superfamilies can overlap the same residue, so overlap counts are not independent evidence votes.
- Database coordinates do not establish mutation effect, binding disruption, conservation, or pathogenicity.
- InterPro availability and response time become another optional external dependency.
- The coverage score composition changes and must not be compared numerically with older six-layer reports.

## Alternatives Considered

- Use only UniProt domain text: rejected because it omits direct InterPro/Pfam pagination and normalized member coordinates.
- Query only InterPro integrated entries: rejected because Pfam-specific accessions, models, and scores would be lost.
- Turn every overlap into a functional-impact claim: rejected because coordinate containment alone cannot support that conclusion.
- Make domain retrieval mandatory: rejected because comprehensive analysis should degrade when optional databases fail.

## References

- https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/P04637/
- https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/P04637/

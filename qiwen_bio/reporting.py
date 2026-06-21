from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from qiwen_bio.models import AnalysisResponse
    from qiwen_bio.pubmed import LiteratureEvidence
    from qiwen_bio.synthesis import ComprehensiveAnalysis


def render_literature_section(evidence: "LiteratureEvidence") -> str:
    if evidence.articles:
        citations = "\n".join(
            (
                f"- {article.title.rstrip('.')} — [PMID {article.pmid}]({article.url}) — "
                f"{', '.join(article.authors[:3]) or 'Unknown authors'}; "
                f"*{article.journal}* ({article.published})."
            )
            for article in evidence.articles
        )
    else:
        citations = "- No PubMed records matched this query."
    return f"""## PubMed literature retrieval

**Query:** `{evidence.query}`

{citations}

> {evidence.disclaimer}
"""


def render_markdown_report(result: "AnalysisResponse") -> str:
    evidence_rows = "\n".join(
        f"| {item.stage} | {item.claim} | {item.source} | {item.confidence} |"
        for item in result.evidence_chain
    )
    limitations = "\n".join(f"- {item}" for item in result.prediction.limitations)
    return f"""# Qiwen Bio analysis: {result.name}

**Analysis ID:** `{result.analysis_id}`  
**Sequence length:** {result.features.length} aa  
**Mutation:** {result.mutation or "Not supplied"}

## Calculated properties

- Molecular weight: {result.features.molecular_weight_da:,.2f} Da
- Hydrophobic fraction: {result.features.hydrophobic_fraction:.1%}
- Charged fraction: {result.features.charged_fraction:.1%}
- Net-charge proxy: {result.features.net_charge_proxy:+d}

## Demo prediction

The transparent AMP baseline returned **{result.prediction.label}** with a score of **{result.prediction.score:.3f}**.

## Evidence chain

| Stage | Claim | Source | Confidence |
| --- | --- | --- | --- |
{evidence_rows}

## Limitations

{limitations}

## Suggested next validation

Replace the demo predictor with a classifier trained on a curated AMP dataset, split by sequence similarity, then report AUROC, AUPRC, calibration, and an external test result. Candidate activity should be confirmed with an appropriate antimicrobial assay.
"""


def render_comprehensive_report(result: "ComprehensiveAnalysis") -> str:
    coverage_rows = "\n".join(
        f"| {item.name} | {'yes' if item.available else 'no'} | {item.points}/{item.max_points} | {item.detail} |"
        for item in result.coverage.components
    )
    function_text = result.annotation.functions[0] if result.annotation.functions else "No function text available."
    sections = [
        f"""# Qiwen Bio comprehensive report: {result.annotation.protein_name}

- Analysis ID: `{result.analysis.analysis_id}`
- UniProt: [{result.annotation.accession}]({result.annotation.source_url})
- Organism: {result.annotation.organism}
- Mutation: {result.analysis.mutation or "Not supplied"}

## Evidence coverage: {result.coverage.score}/100 ({result.coverage.label})

> {result.coverage.interpretation}

| Layer | Available | Points | Detail |
| --- | --- | ---: | --- |
{coverage_rows}

## UniProt annotation

{function_text}

- Gene names: {', '.join(result.annotation.gene_names) or 'Not available'}
- GO terms: {len(result.annotation.go_terms)}

## Sequence evidence

- Length: {result.analysis.features.length} aa
- Molecular weight: {result.analysis.features.molecular_weight_da:,.2f} Da
- Hydrophobic fraction: {result.analysis.features.hydrophobic_fraction:.1%}
- Net-charge proxy: {result.analysis.features.net_charge_proxy:+d}
"""
    ]
    if result.structure:
        mutation_line = ""
        if result.structure.mutation_site:
            site = result.structure.mutation_site
            mutation_line = (
                f"- Mutation site {site.wild_type}{site.position}{site.mutant}: "
                f"pLDDT {site.plddt:.2f} ({site.confidence})\n"
                f"- CA neighbors within 8 A: {len(result.structure.mutation_neighborhood)}\n"
            )
        sections.append(
            f"""## Structure evidence

- AlphaFold mean pLDDT: {result.structure.mean_plddt:.2f}
- Residues represented: {result.structure.residue_count}
- Non-local CA contacts: {result.structure.contact_map.total_contacts}
{mutation_line}
> pLDDT and CA proximity do not predict pathogenicity, stability, or functional effect.
"""
        )
    if result.graph:
        interactions = sum(edge.type == "interacts_with" for edge in result.graph.edges)
        terms = [node for node in result.graph.nodes if node.type != "protein"]
        term_lines = "\n".join(
            f"- {node.label} (`{node.external_id}`), FDR {node.fdr:.2e}"
            for node in terms[:6]
            if node.fdr is not None
        )
        sections.append(
            f"""## Interaction and pathway evidence

- STRING protein nodes: {sum(node.type == 'protein' for node in result.graph.nodes)}
- STRING interaction edges: {interactions}
- Process/pathway terms: {len(terms)}

{term_lines}

> STRING interactions and enrichment are associative database evidence, not causal claims.
"""
        )
    if result.literature:
        sections.append(render_literature_section(result.literature))
    if result.warnings:
        warning_lines = "\n".join(f"- {warning}" for warning in result.warnings)
        sections.append(f"## Unavailable evidence layers\n\n{warning_lines}\n")
    sections.append(
        f"""## Demo predictor limitation

The AMP baseline returned **{result.analysis.prediction.label}** with score **{result.analysis.prediction.score:.3f}**. It is an untrained engineering baseline and must not be used as biological or clinical evidence.

## Interpretation boundary

This report organizes retrieved and calculated evidence. It does not establish a phenotype, mechanism, pathogenicity classification, or experimental outcome. Those require claim-level literature appraisal and experimental validation.
"""
    )
    return "\n".join(sections)

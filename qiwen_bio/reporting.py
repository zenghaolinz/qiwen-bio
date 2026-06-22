from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from qiwen_bio.models import AnalysisResponse
    from qiwen_bio.phenotype_literature import PhenotypeLiteratureEvidence
    from qiwen_bio.pubmed import LiteratureEvidence
    from qiwen_bio.reasoning_chain import ReasoningChain
    from qiwen_bio.structure_features import StructureFeatureSummary
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


def render_phenotype_literature_section(
    evidence: "PhenotypeLiteratureEvidence",
) -> str:
    if evidence.process_links:
        link_lines = []
        for link in evidence.process_links:
            phrases = (
                f" · verbs: {', '.join(link.matched_phrases)}"
                if link.matched_phrases
                else ""
            )
            link_lines.append(
                f"- [{link.support_level}] `{link.process_id}` {link.process_label} — "
                f"[PMID {link.pmid}](https://pubmed.ncbi.nlm.nih.gov/{link.pmid}/) "
                f"{link.title.rstrip('.')}{phrases} — {link.evidence_basis}"
            )
        links_text = "\n".join(link_lines)
    else:
        links_text = "- No phenotype literature links were classified for the directly supported processes."
    if evidence.hypotheses:
        hypothesis_lines = "\n".join(f"- {hypothesis}" for hypothesis in evidence.hypotheses)
    else:
        hypothesis_lines = "- No phenotype hypothesis was emitted: no directly supported process had an article classified as 'supports' at the abstract-metadata boundary."
    counts_text = ", ".join(
        f"{level}={count}" for level, count in evidence.counts_by_level.items()
    ) or "none"
    return f"""## Phenotype literature evidence

- Gene: {evidence.gene}
- Claim-support counts: {counts_text}

### Process-link classifications

{links_text}

### Phenotype hypotheses

{hypothesis_lines}

> {evidence.disclaimer}
> {evidence.boundary}
"""


def render_structure_summary_section(summary: "StructureFeatureSummary") -> str:
    if not summary.has_structure:
        return (
            "## Structure evidence summary\n\n"
            "- No AlphaFold structure is available for this entry.\n"
            "- Structure layer: missing.\n\n"
            "> No structural context is available; geometric interpretation is not possible.\n"
        )
    lines = [
        "## Structure evidence summary",
        "",
        f"- Structure source: AlphaFold predicted structure ({summary.source})",
        f"- Mean pLDDT: {summary.mean_plddt:.2f}" if summary.mean_plddt is not None else "- Mean pLDDT: unavailable",
    ]
    if summary.mutation is not None and summary.mutation_site_plddt is not None:
        lines.append(f"- Mutation site pLDDT: {summary.mutation_site_plddt:.2f}")
        lines.append(f"- Mutation confidence band: {summary.mutation_site_confidence_band}")
        if summary.low_confidence_region:
            lines.append("- Low-confidence predicted region: yes (structural interpretation limited)")
        else:
            lines.append("- Low-confidence predicted region: no")
    else:
        lines.append("- Mutation site pLDDT: not applicable (no mutation supplied)")
    if summary.domain_overlap is not None:
        lines.append(
            f"- Domain overlap: {'yes' if summary.domain_overlap else 'no'}"
            + (f" ({', '.join(summary.overlapping_domains)})" if summary.overlapping_domains else "")
        )
    if summary.contact_count_8a is not None:
        lines.append(f"- 8 A CA contacts: {summary.contact_count_8a}")
    if summary.neighbor_count_8a is not None:
        lines.append(f"- 8 A CA neighbors: {summary.neighbor_count_8a}")
    limits = "\n".join(f"  - {limit}" for limit in summary.interpretation_limits)
    return f"{chr(10).join(lines)}\n\nInterpretation limits:\n{limits}\n\n> AlphaFold pLDDT is local model confidence, not pathogenicity or functional-effect confidence. CA contacts are geometric proximity, not confirmed biochemical interactions.\n"


def render_reasoning_chain_section(chain: "ReasoningChain") -> str:
    step_blocks = []
    for step in chain.steps:
        status = "available" if step.available else "unavailable"
        facts = "\n".join(f"  - {fact}" for fact in step.evidence_facts) or "  - (no facts)"
        hypothesis_line = (
            f"  - 假设：{step.hypothesis}" if step.hypothesis else "  - (no hypothesis emitted)"
        )
        sources = ", ".join(step.evidence_sources) or "none"
        step_blocks.append(
            f"### {step.title} (`{step.step_id}`) — {status} · confidence: {step.confidence}\n\n"
            f"Facts:\n{facts}\n\n"
            f"{hypothesis_line}\n\n"
            f"Sources: {sources}\n\n"
            f"Uncertainty: {step.uncertainty}"
        )
    steps_text = "\n\n".join(step_blocks)
    missing = ", ".join(chain.missing_layers) or "none"
    return f"""## Reasoning chain

- Chain type: `{chain.chain_type}`
- Gene: {chain.gene}
- Mutation: {chain.mutation or "not supplied"}
- Missing layers: {missing}
- Summary: {chain.summary}

{steps_text}

> {chain.boundary}
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
    if result.reasoning_chain is not None and result.reasoning_chain.structure_summary is not None:
        sections.append(
            render_structure_summary_section(result.reasoning_chain.structure_summary)
        )
    if result.domains:
        domain_lines = []
        for entry in result.domains.entries[:12]:
            coordinates = ", ".join(
                f"{location.start}-{location.end}" for location in entry.locations
            )
            overlap = (
                f"; overlaps mutation position {result.domains.mutation_position}"
                if entry.overlaps_mutation
                else ""
            )
            domain_lines.append(
                f"- [{entry.accession}]({entry.source_url}) {entry.name} "
                f"({entry.source_database}; {coordinates}{overlap})"
            )
        sections.append(
            f"""## Domain evidence

- Direct InterPro/Pfam entries: {result.domains.entry_count}
- Annotated coordinate fragments: {result.domains.location_count}
- Entries overlapping the supplied mutation position: {len(result.domains.mutation_overlaps)}

{chr(10).join(domain_lines)}

> {result.domains.disclaimer}
"""
        )
    if result.kegg:
        pathway_lines = "\n".join(
            f"- [{pathway.pathway_id}]({pathway.source_url}) {pathway.name}"
            for pathway in result.kegg.pathways[:12]
        )
        sections.append(
            f"""## Direct KEGG pathway evidence

- KEGG gene IDs: {', '.join(result.kegg.gene_ids)}
- Returned pathway records: {result.kegg.pathway_count} of {result.kegg.linked_pathway_count} linked

{pathway_lines}

> {result.kegg.disclaimer}
> KEGG content is retrieved on demand for this report; it is not bundled or redistributed by Qiwen Bio.
"""
        )
    if result.cellular_processes.processes:
        process_lines = []
        for process in result.cellular_processes.processes[:12]:
            support_text = ", ".join(
                (
                    f"{support.source_name} ({support.evidence_type}"
                    f"; FDR {support.fdr:.2e})"
                    if support.fdr is not None
                    else f"{support.source_name} ({support.evidence_type})"
                )
                for support in process.supports
            )
            process_lines.append(
                f"- `{process.canonical_id}` {process.label}: {support_text}"
            )
        if result.phenotype_literature is not None and result.phenotype_literature.hypotheses:
            hypothesis_summary = (
                "Phenotype hypotheses were emitted for directly supported processes "
                "with abstract-level supporting literature; see the Phenotype literature "
                "evidence section."
            )
        else:
            hypothesis_summary = (
                "No phenotype hypothesis was emitted. Process association does not "
                "establish activity, direction, mechanism, causality, or phenotype."
            )
        sections.append(
            f"""## Cellular-process evidence

- Normalized process/pathway records: {len(result.cellular_processes.processes)}
- Supports by source: {', '.join(f'{name}={count}' for name, count in result.cellular_processes.counts_by_source.items())}

{chr(10).join(process_lines)}

{hypothesis_summary}

> {result.cellular_processes.interpretation_boundary}
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
    if result.phenotype_literature is not None:
        sections.append(render_phenotype_literature_section(result.phenotype_literature))
    if result.reasoning_chain is not None:
        sections.append(render_reasoning_chain_section(result.reasoning_chain))
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

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from qiwen_bio.models import AnalysisResponse
    from qiwen_bio.pubmed import LiteratureEvidence


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

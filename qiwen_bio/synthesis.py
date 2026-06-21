from typing import Literal

from pydantic import BaseModel

from qiwen_bio.alphafold import (
    AlphaFoldAnalysis,
    AlphaFoldNotFoundError,
    AlphaFoldServiceError,
    MutationMismatchError,
)
from qiwen_bio.models import AnalysisRequest, AnalysisResponse
from qiwen_bio.pubmed import LiteratureEvidence, PubMedServiceError
from qiwen_bio.stringdb import EvidenceGraph, StringNotFoundError, StringServiceError
from qiwen_bio.uniprot import UniProtAnnotation


class CoverageComponent(BaseModel):
    name: str
    available: bool
    points: int
    max_points: int
    detail: str


class EvidenceCoverage(BaseModel):
    score: int
    max_score: int = 100
    label: Literal["limited", "partial", "comprehensive"]
    components: list[CoverageComponent]
    interpretation: str = (
        "Coverage measures which evidence layers are present; it is not a confidence, "
        "pathogenicity, or correctness probability."
    )


class ComprehensiveAnalysis(BaseModel):
    annotation: UniProtAnnotation
    analysis: AnalysisResponse
    structure: AlphaFoldAnalysis | None
    graph: EvidenceGraph | None
    literature: LiteratureEvidence | None
    coverage: EvidenceCoverage
    warnings: list[str]
    report_markdown: str


def _component(name: str, available: bool, maximum: int, detail: str) -> CoverageComponent:
    return CoverageComponent(
        name=name,
        available=available,
        points=maximum if available else 0,
        max_points=maximum,
        detail=detail,
    )


def build_comprehensive_analysis(
    identifier: str,
    organism_id: int,
    mutation: str | None,
    pipeline,
    uniprot_client,
    alphafold_client,
    string_client,
    pubmed_client,
    string_limit: int = 8,
    required_score: int = 700,
    literature_limit: int = 5,
) -> ComprehensiveAnalysis:
    annotation = uniprot_client.resolve(identifier, organism_id)
    analysis = pipeline.analyze(
        AnalysisRequest(
            name=annotation.protein_name,
            sequence=annotation.sequence,
            mutation=mutation,
        )
    )
    warnings: list[str] = []

    structure = None
    if annotation.alphafold_url:
        try:
            structure = alphafold_client.analyze(annotation.accession, mutation)
        except (AlphaFoldNotFoundError, AlphaFoldServiceError, MutationMismatchError) as exc:
            warnings.append(f"AlphaFold: {exc}")
    else:
        warnings.append("AlphaFold: no model is linked from the UniProt record")

    graph = None
    try:
        graph = string_client.build_graph(
            identifier=annotation.gene_names[0] if annotation.gene_names else identifier,
            species=organism_id,
            limit=string_limit,
            required_score=required_score,
        )
    except (StringNotFoundError, StringServiceError) as exc:
        warnings.append(f"STRING: {exc}")

    context_terms = []
    if graph:
        context_terms = [node.label for node in graph.nodes if node.type != "protein"][:3]
    literature = None
    try:
        literature = pubmed_client.search(
            protein=annotation.gene_names[0] if annotation.gene_names else identifier,
            context_terms=context_terms,
            limit=literature_limit,
        )
    except PubMedServiceError as exc:
        warnings.append(f"PubMed: {exc}")

    interaction_available = bool(
        graph and any(edge.type == "interacts_with" for edge in graph.edges)
    )
    term_available = bool(graph and any(node.type != "protein" for node in graph.nodes))
    literature_available = bool(literature and literature.articles)
    components = [
        _component("sequence", True, 15, f"{analysis.features.length} canonical residues parsed"),
        _component("annotation", True, 20, f"Reviewed UniProt entry {annotation.accession}"),
        _component("structure", structure is not None, 20, "AlphaFold pLDDT and CA geometry"),
        _component("interactions", interaction_available, 15, "STRING scored interaction edges"),
        _component("pathways", term_available, 15, "STRING process/pathway enrichment terms"),
        _component("literature", literature_available, 15, "Context-bound PubMed records"),
    ]
    score = sum(item.points for item in components)
    label: Literal["limited", "partial", "comprehensive"]
    if score >= 80:
        label = "comprehensive"
    elif score >= 50:
        label = "partial"
    else:
        label = "limited"
    coverage = EvidenceCoverage(score=score, label=label, components=components)
    result = ComprehensiveAnalysis(
        annotation=annotation,
        analysis=analysis,
        structure=structure,
        graph=graph,
        literature=literature,
        coverage=coverage,
        warnings=warnings,
        report_markdown="",
    )
    from qiwen_bio.reporting import render_comprehensive_report

    result.report_markdown = render_comprehensive_report(result)
    return result


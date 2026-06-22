from typing import Literal

from pydantic import BaseModel

from qiwen_bio.alphafold import (
    AlphaFoldAnalysis,
    AlphaFoldNotFoundError,
    AlphaFoldServiceError,
    MutationMismatchError,
)
from qiwen_bio.models import AnalysisRequest, AnalysisResponse
from qiwen_bio.interpro import (
    DomainAnnotation,
    InterProNotFoundError,
    InterProServiceError,
)
from qiwen_bio.kegg import (
    KeggNotFoundError,
    KeggPathwayAnnotation,
    KeggServiceError,
)
from qiwen_bio.cellular_processes import (
    CellularProcessEvidence,
    build_cellular_process_evidence,
)
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
    domains: DomainAnnotation | None
    kegg: KeggPathwayAnnotation | None
    cellular_processes: CellularProcessEvidence
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
    interpro_client,
    kegg_client,
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

    mutation_position = None
    if mutation and mutation[1:-1].isdigit():
        mutation_position = int(mutation[1:-1])
    domains = None
    try:
        domains = interpro_client.fetch(
            annotation.accession, mutation_position=mutation_position
        )
    except (InterProNotFoundError, InterProServiceError) as exc:
        warnings.append(f"InterPro/Pfam: {exc}")

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

    kegg = None
    try:
        kegg = kegg_client.fetch(annotation.accession, limit=20)
    except (KeggNotFoundError, KeggServiceError) as exc:
        warnings.append(f"KEGG: {exc}")

    cellular_processes = build_cellular_process_evidence(annotation, kegg, graph)
    context_terms = cellular_processes.literature_context[:3]
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
    string_term_available = bool(
        graph and any(node.type != "protein" for node in graph.nodes)
    )
    direct_pathway_available = bool(kegg and kegg.pathways)
    pathway_available = direct_pathway_available or string_term_available
    pathway_detail = (
        f"{kegg.pathway_count} direct KEGG pathway records; "
        f"{len(cellular_processes.processes)} normalized process/pathway records"
        if direct_pathway_available
        else "STRING enrichment fallback (no direct KEGG records); "
        f"{len(cellular_processes.processes)} normalized process/pathway records"
        if string_term_available
        else "No direct KEGG or STRING pathway/process evidence"
    )
    literature_available = bool(literature and literature.articles)
    components = [
        _component("sequence", True, 10, f"{analysis.features.length} canonical residues parsed"),
        _component("annotation", True, 15, f"Reviewed UniProt entry {annotation.accession}"),
        _component(
            "domains",
            domains is not None,
            10,
            f"{domains.entry_count if domains else 0} direct InterPro/Pfam entries",
        ),
        _component("structure", structure is not None, 20, "AlphaFold pLDDT and CA geometry"),
        _component("interactions", interaction_available, 15, "STRING scored interaction edges"),
        _component("pathways", pathway_available, 15, pathway_detail),
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
        domains=domains,
        kegg=kegg,
        cellular_processes=cellular_processes,
        coverage=coverage,
        warnings=warnings,
        report_markdown="",
    )
    from qiwen_bio.reporting import render_comprehensive_report

    result.report_markdown = render_comprehensive_report(result)
    return result

from pathlib import Path

from fastapi import Depends, FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from qiwen_bio import __version__
from qiwen_bio.alphafold import (
    AlphaFoldAnalysis,
    AlphaFoldClient,
    AlphaFoldNotFoundError,
    AlphaFoldServiceError,
    MutationMismatchError,
)
from qiwen_bio.models import (
    AlphaFoldAnalysisRequest,
    AnalysisRequest,
    AnalysisResponse,
    ComprehensiveReportRequest,
    EvidenceItem,
    PubMedSearchRequest,
    StringGraphRequest,
    UniProtAnalysisRequest,
)
from qiwen_bio.pipeline import AnalysisPipeline
from qiwen_bio.pubmed import LiteratureEvidence, PubMedClient, PubMedServiceError
from qiwen_bio.reporting import render_literature_section, render_markdown_report
from qiwen_bio.stringdb import (
    EvidenceGraph,
    StringClient,
    StringNotFoundError,
    StringServiceError,
)
from qiwen_bio.synthesis import ComprehensiveAnalysis, build_comprehensive_analysis
from qiwen_bio.uniprot import (
    AmbiguousProteinError,
    ProteinNotFoundError,
    UniProtAnnotation,
    UniProtClient,
    UniProtServiceError,
)


STATIC_DIR = Path(__file__).parent / "static"
pipeline = AnalysisPipeline()
uniprot_client = UniProtClient()
alphafold_client = AlphaFoldClient()
string_client = StringClient()
pubmed_client = PubMedClient()
app = FastAPI(
    title="Qiwen Bio API",
    version=__version__,
    description="Evidence-oriented protein analysis MVP",
)
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")


@app.get("/", include_in_schema=False)
def index() -> FileResponse:
    return FileResponse(STATIC_DIR / "index.html")


@app.get("/api/health")
def health() -> dict[str, str]:
    return {"status": "ok", "version": __version__}


@app.post("/api/v1/analyze", response_model=AnalysisResponse)
def analyze(request: AnalysisRequest) -> AnalysisResponse:
    return pipeline.analyze(request)


def get_uniprot_client() -> UniProtClient:
    return uniprot_client


def get_alphafold_client() -> AlphaFoldClient:
    return alphafold_client


def get_string_client() -> StringClient:
    return string_client


def get_pubmed_client() -> PubMedClient:
    return pubmed_client


class UniProtAnalysisResponse(BaseModel):
    annotation: UniProtAnnotation
    analysis: AnalysisResponse


@app.post("/api/v1/analyze/uniprot", response_model=UniProtAnalysisResponse)
def analyze_uniprot(
    request: UniProtAnalysisRequest,
    client: UniProtClient = Depends(get_uniprot_client),
) -> UniProtAnalysisResponse:
    try:
        annotation = client.resolve(request.identifier, request.organism_id)
    except ProteinNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except AmbiguousProteinError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    except UniProtServiceError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc

    analysis = pipeline.analyze(
        AnalysisRequest(
            name=annotation.protein_name,
            sequence=annotation.sequence,
            mutation=request.mutation,
        )
    )
    analysis.evidence_chain.insert(
        0,
        EvidenceItem(
            stage="annotation",
            claim=(
                f"Resolved {request.identifier.upper()} to reviewed UniProt entry "
                f"{annotation.accession} ({annotation.protein_name})."
            ),
            evidence_type="database_record",
            source=annotation.source_url,
            confidence="high",
        ),
    )
    analysis.report_markdown = render_markdown_report(analysis)
    return UniProtAnalysisResponse(annotation=annotation, analysis=analysis)


class AlphaFoldAnalysisResponse(BaseModel):
    structure: AlphaFoldAnalysis
    interpretation: str


@app.post("/api/v1/structure/alphafold", response_model=AlphaFoldAnalysisResponse)
def analyze_alphafold(
    request: AlphaFoldAnalysisRequest,
    client: AlphaFoldClient = Depends(get_alphafold_client),
) -> AlphaFoldAnalysisResponse:
    try:
        structure = client.analyze(request.accession, request.mutation)
    except AlphaFoldNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except MutationMismatchError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except AlphaFoldServiceError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc
    return AlphaFoldAnalysisResponse(
        structure=structure,
        interpretation=(
            "pLDDT measures local model confidence; it is not a pathogenicity, "
            "stability, or functional-effect prediction. Contact pairs and mutation "
            "neighbors indicate CA geometric proximity only, not biochemical interaction."
        ),
    )


@app.post("/api/v1/graph/string", response_model=EvidenceGraph)
def build_string_graph(
    request: StringGraphRequest,
    client: StringClient = Depends(get_string_client),
) -> EvidenceGraph:
    try:
        return client.build_graph(
            identifier=request.identifier,
            species=request.organism_id,
            limit=request.limit,
            required_score=request.required_score,
        )
    except StringNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except StringServiceError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc


class PubMedSearchResponse(BaseModel):
    evidence: LiteratureEvidence
    markdown_section: str


@app.post("/api/v1/literature/pubmed", response_model=PubMedSearchResponse)
def search_pubmed(
    request: PubMedSearchRequest,
    client: PubMedClient = Depends(get_pubmed_client),
) -> PubMedSearchResponse:
    try:
        evidence = client.search(
            protein=request.protein,
            context_terms=request.context_terms,
            limit=request.limit,
        )
    except PubMedServiceError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc
    return PubMedSearchResponse(
        evidence=evidence,
        markdown_section=render_literature_section(evidence),
    )


@app.post("/api/v1/report/comprehensive", response_model=ComprehensiveAnalysis)
def comprehensive_report(
    request: ComprehensiveReportRequest,
    uniprot: UniProtClient = Depends(get_uniprot_client),
    alphafold: AlphaFoldClient = Depends(get_alphafold_client),
    string: StringClient = Depends(get_string_client),
    pubmed: PubMedClient = Depends(get_pubmed_client),
) -> ComprehensiveAnalysis:
    try:
        return build_comprehensive_analysis(
            identifier=request.identifier,
            organism_id=request.organism_id,
            mutation=request.mutation,
            pipeline=pipeline,
            uniprot_client=uniprot,
            alphafold_client=alphafold,
            string_client=string,
            pubmed_client=pubmed,
            string_limit=request.interaction_limit,
            required_score=request.required_score,
            literature_limit=request.literature_limit,
        )
    except ProteinNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except AmbiguousProteinError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    except UniProtServiceError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc

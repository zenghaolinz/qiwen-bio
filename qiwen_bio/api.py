from pathlib import Path

from fastapi import Depends, FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from qiwen_bio import __version__
from qiwen_bio.models import AnalysisRequest, AnalysisResponse, EvidenceItem, UniProtAnalysisRequest
from qiwen_bio.pipeline import AnalysisPipeline
from qiwen_bio.reporting import render_markdown_report
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

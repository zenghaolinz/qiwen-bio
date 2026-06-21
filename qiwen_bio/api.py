from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from qiwen_bio import __version__
from qiwen_bio.models import AnalysisRequest, AnalysisResponse
from qiwen_bio.pipeline import AnalysisPipeline


STATIC_DIR = Path(__file__).parent / "static"
pipeline = AnalysisPipeline()
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


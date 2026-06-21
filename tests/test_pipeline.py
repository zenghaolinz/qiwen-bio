import pytest
from pydantic import ValidationError

from qiwen_bio.models import AnalysisRequest
from qiwen_bio.pipeline import AnalysisPipeline


def test_sequence_is_normalized_and_analyzed() -> None:
    request = AnalysisRequest(name="Demo", sequence=" kwk lfkkl ")
    result = AnalysisPipeline().analyze(request)

    assert result.sequence == "KWKLFKKL"
    assert result.features.length == 8
    assert result.features.net_charge_proxy == 4
    assert result.prediction.model_version.endswith("-demo")
    assert "untrained heuristic" in result.report_markdown
    assert result.analysis_id == AnalysisPipeline().analyze(request).analysis_id


def test_invalid_amino_acid_is_rejected() -> None:
    with pytest.raises(ValidationError, match="unsupported amino-acid symbols: B, X"):
        AnalysisRequest(sequence="ACBX")


def test_mutation_is_preserved_as_unvalidated_evidence() -> None:
    result = AnalysisPipeline().analyze(
        AnalysisRequest(sequence="MKTIIALSYIFCLVFA", mutation="A12V")
    )

    mutation = next(item for item in result.evidence_chain if item.stage == "mutation")
    assert mutation.evidence_type == "user_input"
    assert mutation.confidence == "low"

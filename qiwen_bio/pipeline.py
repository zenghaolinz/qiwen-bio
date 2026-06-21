from hashlib import sha256

from qiwen_bio.features import extract_sequence_features
from qiwen_bio.models import AnalysisRequest, AnalysisResponse, EvidenceItem
from qiwen_bio.predictors import DemoAmpPredictor, ProteinPredictor
from qiwen_bio.reporting import render_markdown_report


class AnalysisPipeline:
    def __init__(self, predictor: ProteinPredictor | None = None) -> None:
        self.predictor = predictor or DemoAmpPredictor()

    def analyze(self, request: AnalysisRequest) -> AnalysisResponse:
        features = extract_sequence_features(request.sequence)
        prediction = self.predictor.predict(features)
        name = request.name or "Unnamed protein"
        analysis_id = sha256(
            f"{name}|{request.sequence}|{request.mutation or ''}".encode()
        ).hexdigest()[:12]
        evidence = [
            EvidenceItem(
                stage="sequence",
                claim=f"The input contains {features.length} canonical amino acids.",
                evidence_type="calculated",
                source="Qiwen Bio deterministic sequence parser",
                confidence="high",
            ),
            EvidenceItem(
                stage="properties",
                claim=(
                    f"Hydrophobic fraction is {features.hydrophobic_fraction:.1%}; "
                    f"net-charge proxy is {features.net_charge_proxy:+d}."
                ),
                evidence_type="calculated",
                source="Qiwen Bio residue-property calculator",
                confidence="high",
            ),
            EvidenceItem(
                stage="prediction",
                claim=(
                    f"Demo AMP baseline classified the sequence as {prediction.label} "
                    f"with score {prediction.score:.3f}."
                ),
                evidence_type="model_output",
                source=f"{prediction.model_name}@{prediction.model_version}",
                confidence="low",
            ),
        ]
        if request.mutation:
            evidence.insert(1, EvidenceItem(
                stage="mutation",
                claim=f"Mutation supplied by user: {request.mutation}.",
                evidence_type="user_input",
                source="User input; not structurally validated in MVP",
                confidence="low",
            ))

        response = AnalysisResponse(
            analysis_id=analysis_id,
            name=name,
            sequence=request.sequence,
            mutation=request.mutation,
            features=features,
            prediction=prediction,
            evidence_chain=evidence,
            report_markdown="",
        )
        response.report_markdown = render_markdown_report(response)
        return response


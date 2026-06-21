import math
from typing import Protocol

from qiwen_bio.models import Prediction, SequenceFeatures


class ProteinPredictor(Protocol):
    def predict(self, features: SequenceFeatures) -> Prediction: ...


class DemoAmpPredictor:
    """Transparent baseline that can later be replaced by an ESM classifier."""

    def predict(self, features: SequenceFeatures) -> Prediction:
        length_fit = max(0.0, 1.0 - abs(features.length - 25) / 35)
        charge_fit = min(max(features.net_charge_proxy, 0) / 6, 1.0)
        hydrophobic_fit = max(0.0, 1.0 - abs(features.hydrophobic_fraction - 0.45) / 0.45)
        raw = 2.2 * length_fit + 2.0 * charge_fit + 1.6 * hydrophobic_fit - 3.2
        score = round(1 / (1 + math.exp(-raw)), 4)
        return Prediction(
            task="antimicrobial_peptide_demo",
            label="candidate" if score >= 0.5 else "unlikely",
            score=score,
            model_name="transparent_amp_heuristic",
            model_version="0.1.0-demo",
            limitations=[
                "This is an untrained heuristic baseline, not a validated biological model.",
                "The result must not be used for clinical or experimental decisions.",
                "A trained classifier with homology-aware validation is required for research use.",
            ],
        )


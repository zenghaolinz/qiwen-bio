from typing import Literal

from pydantic import BaseModel, Field, field_validator


AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


class AnalysisRequest(BaseModel):
    sequence: str = Field(min_length=1, max_length=10_000)
    name: str | None = Field(default=None, max_length=120)
    mutation: str | None = Field(default=None, max_length=30)

    @field_validator("sequence")
    @classmethod
    def normalize_sequence(cls, value: str) -> str:
        normalized = "".join(value.split()).upper()
        invalid = sorted(set(normalized) - AMINO_ACIDS)
        if invalid:
            raise ValueError(f"unsupported amino-acid symbols: {', '.join(invalid)}")
        return normalized


class SequenceFeatures(BaseModel):
    length: int
    molecular_weight_da: float
    hydrophobic_fraction: float
    charged_fraction: float
    net_charge_proxy: int
    aromatic_fraction: float
    composition: dict[str, float]


class Prediction(BaseModel):
    task: Literal["antimicrobial_peptide_demo"]
    label: Literal["candidate", "unlikely"]
    score: float = Field(ge=0, le=1)
    model_name: str
    model_version: str
    limitations: list[str]


class EvidenceItem(BaseModel):
    stage: str
    claim: str
    evidence_type: Literal["calculated", "model_output", "user_input", "database_record"]
    source: str
    confidence: Literal["high", "medium", "low"]


class AnalysisResponse(BaseModel):
    analysis_id: str
    name: str
    sequence: str
    mutation: str | None
    features: SequenceFeatures
    prediction: Prediction
    evidence_chain: list[EvidenceItem]
    report_markdown: str


class UniProtAnalysisRequest(BaseModel):
    identifier: str = Field(min_length=1, max_length=40, pattern=r"^[A-Za-z0-9_.-]+$")
    organism_id: int = Field(default=9606, gt=0)
    mutation: str | None = Field(default=None, max_length=30)


class AlphaFoldAnalysisRequest(BaseModel):
    accession: str = Field(min_length=6, max_length=10, pattern=r"^[A-Za-z0-9]+$")
    mutation: str | None = Field(
        default=None,
        max_length=30,
        pattern=r"^[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]\d+[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]$",
    )


class StringGraphRequest(BaseModel):
    identifier: str = Field(min_length=1, max_length=40, pattern=r"^[A-Za-z0-9_.-]+$")
    organism_id: int = Field(default=9606, gt=0)
    limit: int = Field(default=10, ge=1, le=20)
    required_score: int = Field(default=700, ge=0, le=1000)


class PubMedSearchRequest(BaseModel):
    protein: str = Field(min_length=1, max_length=40, pattern=r"^[A-Za-z0-9_.-]+$")
    context_terms: list[str] = Field(default_factory=list, max_length=5)
    limit: int = Field(default=5, ge=1, le=10)


class ComprehensiveReportRequest(BaseModel):
    identifier: str = Field(min_length=1, max_length=40, pattern=r"^[A-Za-z0-9_.-]+$")
    organism_id: int = Field(default=9606, gt=0)
    mutation: str | None = Field(default=None, max_length=30)
    interaction_limit: int = Field(default=8, ge=1, le=20)
    required_score: int = Field(default=700, ge=0, le=1000)
    literature_limit: int = Field(default=5, ge=1, le=10)


class EmbeddingRequest(BaseModel):
    sequence: str = Field(min_length=1, max_length=1022)

    @field_validator("sequence")
    @classmethod
    def normalize_embedding_sequence(cls, value: str) -> str:
        normalized = "".join(value.split()).upper()
        invalid = sorted(set(normalized) - AMINO_ACIDS)
        if invalid:
            raise ValueError(f"unsupported amino-acid symbols: {', '.join(invalid)}")
        return normalized

import json
import math
from hashlib import sha256
from pathlib import Path
from typing import Protocol

from pydantic import BaseModel, ValidationError

from qiwen_bio.models import AMINO_ACIDS


ESM2_MODEL_ID = "facebook/esm2_t6_8M_UR50D"
ESM2_REVISION = "c731040fcd8d73dceaa04b0a8e6329b345b0f5df"


class ModelDependencyError(RuntimeError):
    pass


class ModelLoadError(RuntimeError):
    pass


class EmbeddingProvider(Protocol):
    model_id: str
    revision: str
    pooling: str
    dimension: int
    max_residues: int

    def embed(self, sequence: str) -> list[float]: ...


class ProteinEmbedding(BaseModel):
    model_id: str
    model_revision: str
    pooling: str
    sequence_sha256: str
    sequence_length: int
    dimension: int
    vector: list[float]
    cached: bool


class EmbeddingCache:
    def __init__(self, root: Path) -> None:
        self.root = root

    def key(self, provider: EmbeddingProvider, sequence: str) -> str:
        identity = f"{provider.model_id}|{provider.revision}|{provider.pooling}|{sequence}"
        return sha256(identity.encode()).hexdigest()

    def load(self, provider: EmbeddingProvider, sequence: str) -> ProteinEmbedding | None:
        path = self.root / f"{self.key(provider, sequence)}.json"
        if not path.exists():
            return None
        try:
            result = ProteinEmbedding.model_validate_json(path.read_text(encoding="utf-8"))
            expected_sha = sha256(sequence.encode()).hexdigest()
            if (
                result.model_id != provider.model_id
                or result.model_revision != provider.revision
                or result.pooling != provider.pooling
                or result.sequence_sha256 != expected_sha
                or result.dimension != provider.dimension
                or len(result.vector) != provider.dimension
            ):
                raise ValueError("cache metadata mismatch")
            result.cached = True
            return result
        except (OSError, ValueError, ValidationError):
            path.unlink(missing_ok=True)
            return None

    def store_for_sequence(
        self,
        provider: EmbeddingProvider,
        sequence: str,
        result: ProteinEmbedding,
    ) -> None:
        self.root.mkdir(parents=True, exist_ok=True)
        path = self.root / f"{self.key(provider, sequence)}.json"
        temporary = path.with_suffix(".tmp")
        temporary.write_text(result.model_dump_json(), encoding="utf-8")
        temporary.replace(path)


class EmbeddingService:
    def __init__(self, provider: EmbeddingProvider, cache: EmbeddingCache) -> None:
        self.provider = provider
        self.cache = cache

    def embed(self, sequence: str) -> ProteinEmbedding:
        normalized = "".join(sequence.split()).upper()
        invalid = sorted(set(normalized) - AMINO_ACIDS)
        if not normalized or invalid:
            raise ValueError("sequence must contain only canonical amino acids")
        if len(normalized) > self.provider.max_residues:
            raise ValueError(
                f"sequence length {len(normalized)} exceeds provider limit "
                f"{self.provider.max_residues}"
            )
        cached = self.cache.load(self.provider, normalized)
        if cached:
            return cached
        vector = self.provider.embed(normalized)
        if len(vector) != self.provider.dimension or not all(math.isfinite(item) for item in vector):
            raise ValueError("embedding provider returned an invalid vector")
        result = ProteinEmbedding(
            model_id=self.provider.model_id,
            model_revision=self.provider.revision,
            pooling=self.provider.pooling,
            sequence_sha256=sha256(normalized.encode()).hexdigest(),
            sequence_length=len(normalized),
            dimension=self.provider.dimension,
            vector=vector,
            cached=False,
        )
        self.cache.store_for_sequence(self.provider, normalized, result)
        return result


class Esm2EmbeddingProvider:
    model_id = ESM2_MODEL_ID
    revision = ESM2_REVISION
    pooling = "mean"
    dimension = 320
    max_residues = 1022

    def __init__(self, device: str = "auto") -> None:
        self.requested_device = device
        self._tokenizer = None
        self._model = None
        self._torch = None
        self.device = "unloaded"

    def _load(self) -> None:
        if self._model is not None:
            return
        try:
            import torch
            from transformers import AutoModel, AutoTokenizer
        except ImportError as exc:
            raise ModelDependencyError(
                'ESM dependencies are not installed; run: pip install -e ".[model]"'
            ) from exc
        self.device = (
            "cuda" if self.requested_device == "auto" and torch.cuda.is_available()
            else "cpu" if self.requested_device == "auto"
            else self.requested_device
        )
        try:
            self._tokenizer = AutoTokenizer.from_pretrained(self.model_id, revision=self.revision)
            self._model = AutoModel.from_pretrained(
                self.model_id,
                revision=self.revision,
                add_pooling_layer=False,
            )
        except (OSError, RuntimeError) as exc:
            raise ModelLoadError(f"Unable to load pinned ESM-2 model: {exc}") from exc
        self._model.to(self.device)
        self._model.eval()
        self._torch = torch

    def embed(self, sequence: str) -> list[float]:
        self._load()
        inputs = self._tokenizer(sequence, return_tensors="pt")
        inputs = {name: tensor.to(self.device) for name, tensor in inputs.items()}
        with self._torch.inference_mode():
            output = self._model(**inputs).last_hidden_state[0, 1 : len(sequence) + 1]
            pooled = output.mean(dim=0).float().cpu()
        return pooled.tolist()

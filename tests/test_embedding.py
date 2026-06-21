import json

from qiwen_bio.embedding import EmbeddingCache, EmbeddingService


class FakeProvider:
    model_id = "test/esm"
    revision = "revision-a"
    pooling = "mean"
    dimension = 3
    max_residues = 10

    def __init__(self) -> None:
        self.calls = 0

    def embed(self, sequence: str) -> list[float]:
        self.calls += 1
        return [float(len(sequence)), 2.0, 3.0]


def test_embedding_service_caches_by_sequence_and_model_version(tmp_path) -> None:
    provider = FakeProvider()
    service = EmbeddingService(provider, EmbeddingCache(tmp_path))

    first = service.embed("ACDE")
    second = service.embed("ACDE")

    assert provider.calls == 1
    assert first.cached is False
    assert second.cached is True
    assert second.vector == [4.0, 2.0, 3.0]
    assert second.dimension == 3
    assert second.model_revision == "revision-a"


def test_model_revision_change_uses_a_different_cache_entry(tmp_path) -> None:
    first_provider = FakeProvider()
    EmbeddingService(first_provider, EmbeddingCache(tmp_path)).embed("ACDE")
    second_provider = FakeProvider()
    second_provider.revision = "revision-b"

    result = EmbeddingService(second_provider, EmbeddingCache(tmp_path)).embed("ACDE")

    assert second_provider.calls == 1
    assert result.cached is False
    assert len(list(tmp_path.glob("*.json"))) == 2


def test_corrupt_cache_entry_is_recomputed(tmp_path) -> None:
    provider = FakeProvider()
    service = EmbeddingService(provider, EmbeddingCache(tmp_path))
    service.embed("ACDE")
    cache_file = next(tmp_path.glob("*.json"))
    cache_file.write_text(json.dumps({"vector": [999]}), encoding="utf-8")

    result = service.embed("ACDE")

    assert provider.calls == 2
    assert result.cached is False
    assert result.vector == [4.0, 2.0, 3.0]

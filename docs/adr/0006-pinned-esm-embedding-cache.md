# ADR-0006: Pin ESM-2 and Cache Embeddings by Content

## Status

Accepted

## Context

The project needs real protein language-model features before it can train a validated task classifier. Foundation-model outputs can drift when a repository's `main` revision changes, and repeated local inference is wasteful. PyTorch and Transformers are also much larger than the base Web dependencies.

## Decision

Define an `EmbeddingProvider` protocol and implement `Esm2EmbeddingProvider` with `facebook/esm2_t6_8M_UR50D` pinned to commit `c731040fcd8d73dceaa04b0a8e6329b345b0f5df`. Mean-pool residue hidden states while excluding special tokens, producing a 320-dimensional vector. Limit inputs to 1022 residues.

Load PyTorch, Transformers, and model weights lazily through the optional `model` dependency group. Select CUDA only when the installed PyTorch build reports it available; otherwise use CPU. Cache each result as validated JSON under `data/embeddings/`, keyed by model ID, revision, pooling method, and sequence. Use atomic file replacement and recompute corrupt or mismatched entries. Return dependency/loading failures as HTTP 503 and invalid input as 422.

## Consequences

### Positive

- Embeddings are reproducible across model-repository updates.
- Repeated inference becomes a fast local cache read.
- Base Web installation remains lightweight.
- Provider and cache can be tested without loading a neural model.

### Negative

- The JSON vector cache is larger and slower than a binary array format at scale.
- The 8M model is less expressive than larger ESM-2 variants.
- First use requires dependency and model downloads.
- CPU inference is slower when a CUDA-enabled PyTorch build is absent.

### Neutral

- Embeddings are features only; they do not replace a trained, evaluated classifier.

## Alternatives Considered

- Use an unpinned `main` revision: rejected because cached vectors could silently mix model versions.
- Install model dependencies in the base package: rejected because it would burden users who only need database analysis.
- Start with ESM-2 150M: deferred until batching, GPU installation, and memory profiling are established.
- Store vectors in NumPy/SQLite immediately: deferred until dataset-scale training defines access patterns.

## References

- https://huggingface.co/facebook/esm2_t6_8M_UR50D
- https://github.com/facebookresearch/esm

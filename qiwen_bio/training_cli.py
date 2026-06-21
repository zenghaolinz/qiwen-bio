import argparse
from pathlib import Path

from qiwen_bio.embedding import EmbeddingCache, EmbeddingService, Esm2EmbeddingProvider
from qiwen_bio.training import (
    load_embedding_dataset,
    precompute_embedding_dataset,
    render_model_card,
    train_calibrated_baseline,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare and train the offline AMP baseline.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    precompute = subparsers.add_parser("precompute", help="Generate pinned ESM-2 vectors.")
    precompute.add_argument(
        "--csv", type=Path, default=Path("data/datasets/uniprot_amp_samples.csv")
    )
    precompute.add_argument(
        "--manifest", type=Path, default=Path("data/manifests/uniprot_amp_baseline.json")
    )
    precompute.add_argument(
        "--output", type=Path, default=Path("data/training/amp_esm2_embeddings.json")
    )
    precompute.add_argument("--cache", type=Path, default=Path("data/embeddings"))
    precompute.add_argument("--device", choices=("auto", "cpu", "cuda"), default="auto")

    train = subparsers.add_parser("train", help="Fit, calibrate, evaluate, and export JSON.")
    train.add_argument(
        "--embeddings", type=Path, default=Path("data/training/amp_esm2_embeddings.json")
    )
    train.add_argument(
        "--model-output", type=Path, default=Path("models/amp_esm2_logistic.json")
    )
    train.add_argument(
        "--model-card",
        type=Path,
        default=Path("docs/model-cards/amp-esm2-logistic-v0.1.md"),
    )
    train.add_argument("--seed", type=int, default=42)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "precompute":
        service = EmbeddingService(
            Esm2EmbeddingProvider(device=args.device), EmbeddingCache(args.cache)
        )
        artifact = precompute_embedding_dataset(
            csv_path=args.csv,
            manifest_path=args.manifest,
            output_path=args.output,
            embedding_service=service,
        )
        print(
            f"Embedded {len(artifact.samples)} samples with {artifact.model_id} "
            f"at {artifact.model_revision}; output: {args.output}"
        )
        return

    artifact = load_embedding_dataset(args.embeddings)
    result = train_calibrated_baseline(
        artifact, output_path=args.model_output, seed=args.seed
    )
    args.model_card.parent.mkdir(parents=True, exist_ok=True)
    args.model_card.write_text(render_model_card(result.model), encoding="utf-8")
    print(
        f"Test ROC AUC={result.test_metrics.roc_auc:.4f}; "
        f"balanced accuracy={result.test_metrics.balanced_accuracy:.4f}; "
        f"model: {args.model_output}"
    )


if __name__ == "__main__":
    main()

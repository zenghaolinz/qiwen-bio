import argparse
from pathlib import Path

import httpx

from qiwen_bio.dataset import UniProtAmpDatasetClient, write_prepared_dataset


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prepare an auditable UniProt AMP proxy dataset with homology-aware splits."
    )
    parser.add_argument("--positive-limit", type=int, default=50)
    parser.add_argument("--negative-limit", type=int, default=50)
    parser.add_argument("--identity-threshold", type=float, default=0.8)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", type=Path, default=Path("data/datasets/uniprot_amp_samples.csv"))
    parser.add_argument("--manifest", type=Path, default=Path("data/manifests/uniprot_amp_baseline.json"))
    return parser


def main() -> None:
    args = build_parser().parse_args()
    with httpx.Client(timeout=60.0, follow_redirects=True) as http_client:
        prepared = UniProtAmpDatasetClient(http_client=http_client).prepare_dataset(
            positive_limit=args.positive_limit,
            negative_limit=args.negative_limit,
            identity_threshold=args.identity_threshold,
            seed=args.seed,
        )
    write_prepared_dataset(prepared, csv_path=args.output, manifest_path=args.manifest)
    print(
        f"Prepared {prepared.manifest.retained_samples} samples in "
        f"{prepared.manifest.cluster_count} clusters; manifest: {args.manifest}"
    )


if __name__ == "__main__":
    main()

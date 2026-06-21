import argparse
from pathlib import Path

import httpx

from qiwen_bio.mature_peptides import (
    UniProtMaturePeptideClient,
    curate_mature_positive_dataset,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Audit exact UniProt mature-peptide features in AMP positives."
    )
    parser.add_argument(
        "--baseline-csv",
        type=Path,
        default=Path("data/datasets/uniprot_amp_samples.csv"),
    )
    parser.add_argument(
        "--baseline-manifest",
        type=Path,
        default=Path("data/manifests/uniprot_amp_baseline.json"),
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("data/datasets/uniprot_mature_amp_sequences.csv"),
    )
    parser.add_argument(
        "--output-manifest",
        type=Path,
        default=Path("data/manifests/uniprot_mature_amp_audit.json"),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    with httpx.Client(timeout=30.0, follow_redirects=True) as http_client:
        result = curate_mature_positive_dataset(
            baseline_csv=args.baseline_csv,
            baseline_manifest=args.baseline_manifest,
            output_csv=args.output_csv,
            output_manifest=args.output_manifest,
            client=UniProtMaturePeptideClient(http_client=http_client),
        )
    manifest = result.manifest
    print(
        f"Exact mature peptides: {manifest.entries_with_exact_peptide}/"
        f"{manifest.fetch_successes} entries; unique sequences: "
        f"{manifest.unique_mature_sequences}; fetch errors: {manifest.fetch_errors}; "
        f"ready for model replacement: {manifest.ready_for_model_replacement}"
    )


if __name__ == "__main__":
    main()

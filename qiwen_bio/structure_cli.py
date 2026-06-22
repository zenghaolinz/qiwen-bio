"""Structure-feature export CLI (Stage 2, task 6).

Exports a flat structure-feature record for a future controlled
structure-enhanced prediction benchmark. This is NOT a predictor and does NOT
train a model. The record carries structure features only; it never carries
pathogenicity, functional-effect, or stability labels.

Usage:
    python -m qiwen_bio.structure_cli --accession P04637 --mutation R175H \
        --output data/exports/structure_features.jsonl

The export resolves a UniProt accession, fetches the AlphaFold structure and
InterPro domains, builds a StructureFeatureSummary, and writes one JSONL line.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from qiwen_bio.alphafold import AlphaFoldClient
from qiwen_bio.interpro import InterProClient
from qiwen_bio.structure_features import StructureFeatureSummary, build_structure_feature_summary
from qiwen_bio.uniprot import UniProtClient


def to_export_record(
    summary: StructureFeatureSummary,
    sequence_length: int,
    label: int | None = None,
) -> dict:
    """Flatten a StructureFeatureSummary into a benchmark-ready record.

    The record carries structure features only. ``label`` is an optional
    caller-supervised benchmark label (e.g. 0/1); it is never inferred from
    structure and carries no pathogenicity semantics.
    """
    return {
        "accession": summary.accession,
        "sequence_length": sequence_length,
        "mutation": summary.mutation,
        "mean_plddt": summary.mean_plddt,
        "mutation_site_plddt": summary.mutation_site_plddt,
        "mutation_confidence_band": summary.mutation_site_confidence_band,
        "contact_count_8a": summary.contact_count_8a,
        "neighbor_count_8a": summary.neighbor_count_8a,
        "domain_overlap": summary.domain_overlap,
        "overlapping_domains": summary.overlapping_domains,
        "label": label,
        "source": summary.structure_source_type,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Export a structure-feature record for a controlled structure-enhanced "
            "benchmark. Does not train or predict; carries structure features only."
        )
    )
    parser.add_argument("--accession", required=True, help="UniProt accession")
    parser.add_argument("--mutation", default=None, help="Mutation token (e.g. R175H)")
    parser.add_argument("--organism-id", type=int, default=9606)
    parser.add_argument("--label", type=int, default=None, help="Optional supervised benchmark label")
    parser.add_argument("--output", type=Path, default=Path("data/exports/structure_features.jsonl"))
    return parser


def main() -> None:
    args = build_parser().parse_args()
    uniprot = UniProtClient()
    annotation = uniprot.resolve(args.accession, args.organism_id)

    structure = None
    if annotation.alphafold_url:
        try:
            structure = AlphaFoldClient().analyze(annotation.accession, args.mutation)
        except Exception as exc:  # noqa: BLE001 - degrade gracefully
            print(f"warning: AlphaFold unavailable: {exc}", file=sys.stderr)

    domains = None
    try:
        from qiwen_bio.mutation import parse_mutation

        mutation_position = parse_mutation(args.mutation).position
        domains = InterProClient().fetch(annotation.accession, mutation_position=mutation_position)
    except Exception as exc:  # noqa: BLE001 - degrade gracefully
        print(f"warning: InterPro unavailable: {exc}", file=sys.stderr)

    summary = build_structure_feature_summary(
        annotation=annotation,
        structure=structure,
        domains=domains,
        mutation=args.mutation,
    )
    record = to_export_record(summary, sequence_length=len(annotation.sequence), label=args.label)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("a", encoding="utf-8") as f:
        f.write(json.dumps(record) + "\n")
    print(f"wrote 1 record to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()

from __future__ import annotations

import json
from collections import Counter
from csv import DictReader, DictWriter
from hashlib import sha256
from pathlib import Path
from typing import Any, Protocol

import httpx
from pydantic import BaseModel

from qiwen_bio.models import AMINO_ACIDS
from qiwen_bio.dataset import PreparedAmpSample, prepared_samples_sha256


class MaturePeptide(BaseModel):
    accession: str
    feature_index: int
    description: str
    start: int
    end: int
    sequence: str
    parent_sequence_sha256: str


class PeptideFeatureExclusion(BaseModel):
    accession: str
    feature_index: int
    description: str
    reason: str


class MaturePeptideEntry(BaseModel):
    accession: str
    parent_length: int
    peptides: list[MaturePeptide]
    exclusions: list[PeptideFeatureExclusion]


class MaturePeptideClient(Protocol):
    def fetch(self, accession: str) -> MaturePeptideEntry: ...


class CuratedMatureSequence(BaseModel):
    mature_id: str
    sequence: str
    source_accessions: list[str]
    descriptions: list[str]
    instance_count: int


class MaturePeptideFetchFailure(BaseModel):
    accession: str
    error_type: str
    message: str


class MaturePeptideAuditManifest(BaseModel):
    source_endpoint: str
    license_url: str
    baseline_dataset_sha256: str
    input_positive_records: int
    unique_positive_accessions: int
    fetch_successes: int
    fetch_errors: int
    fetch_failures: list[MaturePeptideFetchFailure]
    entries_with_exact_peptide: int
    entries_without_exact_peptide: int
    exact_feature_entry_coverage: float
    excluded_features: int
    exclusion_reasons: dict[str, int]
    extracted_instances: int
    unique_mature_sequences: int
    duplicate_instances: int
    min_mature_length: int | None
    max_mature_length: int | None
    mean_mature_length: float | None
    mature_sequence_sha256: str
    has_exact_mature_positives: bool
    has_defensible_negatives: bool = False
    has_independent_external_benchmark: bool = False
    ready_for_model_replacement: bool = False


class MaturePeptideCurationResult(BaseModel):
    sequences: list[CuratedMatureSequence]
    manifest: MaturePeptideAuditManifest


def parse_mature_peptides(
    payload: dict[str, Any], *, expected_accession: str
) -> MaturePeptideEntry:
    accession = payload.get("primaryAccession")
    if accession != expected_accession:
        raise ValueError(
            f"UniProt response accession {accession!r} does not match {expected_accession!r}"
        )
    sequence = payload.get("sequence", {}).get("value", "")
    if not sequence or set(sequence) - AMINO_ACIDS:
        raise ValueError("UniProt parent sequence is missing or non-canonical")

    peptides: list[MaturePeptide] = []
    exclusions: list[PeptideFeatureExclusion] = []
    for feature_index, feature in enumerate(payload.get("features", [])):
        if feature.get("type") != "Peptide":
            continue
        description = feature.get("description") or ""
        location = feature.get("location") or {}
        start_data = location.get("start") or {}
        end_data = location.get("end") or {}
        if start_data.get("modifier") != "EXACT" or end_data.get("modifier") != "EXACT":
            exclusions.append(
                PeptideFeatureExclusion(
                    accession=accession,
                    feature_index=feature_index,
                    description=description,
                    reason="non_exact_coordinate",
                )
            )
            continue
        start = start_data.get("value")
        end = end_data.get("value")
        if (
            not isinstance(start, int)
            or not isinstance(end, int)
            or start < 1
            or end < start
            or end > len(sequence)
        ):
            exclusions.append(
                PeptideFeatureExclusion(
                    accession=accession,
                    feature_index=feature_index,
                    description=description,
                    reason="out_of_range_coordinate",
                )
            )
            continue
        peptide_sequence = sequence[start - 1 : end]
        if set(peptide_sequence) - AMINO_ACIDS:
            exclusions.append(
                PeptideFeatureExclusion(
                    accession=accession,
                    feature_index=feature_index,
                    description=description,
                    reason="non_canonical_peptide",
                )
            )
            continue
        peptides.append(
            MaturePeptide(
                accession=accession,
                feature_index=feature_index,
                description=description,
                start=start,
                end=end,
                sequence=peptide_sequence,
                parent_sequence_sha256=sha256(sequence.encode()).hexdigest(),
            )
        )
    return MaturePeptideEntry(
        accession=accession,
        parent_length=len(sequence),
        peptides=peptides,
        exclusions=exclusions,
    )


class UniProtMaturePeptideClient:
    endpoint_template = "https://rest.uniprot.org/uniprotkb/{accession}.json"

    def __init__(self, http_client: httpx.Client | None = None) -> None:
        self.http_client = http_client or httpx.Client(timeout=30.0, follow_redirects=True)

    def fetch(self, accession: str) -> MaturePeptideEntry:
        response = self.http_client.get(self.endpoint_template.format(accession=accession))
        response.raise_for_status()
        return parse_mature_peptides(response.json(), expected_accession=accession)


def _read_baseline(csv_path: Path) -> list[PreparedAmpSample]:
    with csv_path.open(encoding="utf-8", newline="") as source:
        return [
            PreparedAmpSample(
                sample_id=row["sample_id"],
                sequence=row["sequence"],
                label=int(row["label"]),
                source_ids=row["source_ids"].split("|") if row["source_ids"] else [],
                sources=row["sources"].split("|") if row["sources"] else [],
                cluster_id=row["cluster_id"],
                split=row["split"],
            )
            for row in DictReader(source)
        ]


def curate_mature_positive_dataset(
    *,
    baseline_csv: Path,
    baseline_manifest: Path,
    output_csv: Path,
    output_manifest: Path,
    client: MaturePeptideClient,
) -> MaturePeptideCurationResult:
    samples = _read_baseline(baseline_csv)
    baseline_metadata = json.loads(baseline_manifest.read_text(encoding="utf-8"))
    baseline_sha = prepared_samples_sha256(samples)
    if baseline_sha != baseline_metadata.get("dataset_sha256"):
        raise ValueError("baseline CSV does not match manifest dataset SHA-256")
    positive_samples = [sample for sample in samples if sample.label == 1]
    accessions = list(
        dict.fromkeys(
            source_id
            for sample in positive_samples
            for source_id in sample.source_ids
        )
    )

    entries: list[MaturePeptideEntry] = []
    failures: list[MaturePeptideFetchFailure] = []
    for accession in accessions:
        try:
            entries.append(client.fetch(accession))
        except (httpx.HTTPError, ValueError, KeyError, TypeError) as exc:
            failures.append(
                MaturePeptideFetchFailure(
                    accession=accession,
                    error_type=type(exc).__name__,
                    message=str(exc),
                )
            )

    instances = [peptide for entry in entries for peptide in entry.peptides]
    by_sequence: dict[str, list[MaturePeptide]] = {}
    for peptide in instances:
        by_sequence.setdefault(peptide.sequence, []).append(peptide)
    curated = [
        CuratedMatureSequence(
            mature_id=f"AMP-MATURE-{index:04d}",
            sequence=sequence,
            source_accessions=list(dict.fromkeys(item.accession for item in peptides)),
            descriptions=list(
                dict.fromkeys(item.description for item in peptides if item.description)
            ),
            instance_count=len(peptides),
        )
        for index, (sequence, peptides) in enumerate(by_sequence.items(), start=1)
    ]
    digest_payload = [item.model_dump(mode="json") for item in curated]
    mature_sha = sha256(json.dumps(digest_payload, sort_keys=True).encode()).hexdigest()
    exclusions = [item for entry in entries for item in entry.exclusions]
    entries_with = sum(bool(entry.peptides) for entry in entries)
    entries_without = len(entries) - entries_with
    lengths = [len(item.sequence) for item in curated]
    reasons = Counter(item.reason for item in exclusions)
    if entries_without:
        reasons["no_exact_peptide_feature"] += entries_without
    manifest = MaturePeptideAuditManifest(
        source_endpoint=UniProtMaturePeptideClient.endpoint_template,
        license_url="https://www.uniprot.org/help/license",
        baseline_dataset_sha256=baseline_sha,
        input_positive_records=len(positive_samples),
        unique_positive_accessions=len(accessions),
        fetch_successes=len(entries),
        fetch_errors=len(failures),
        fetch_failures=failures,
        entries_with_exact_peptide=entries_with,
        entries_without_exact_peptide=entries_without,
        exact_feature_entry_coverage=entries_with / len(entries) if entries else 0.0,
        excluded_features=len(exclusions),
        exclusion_reasons=dict(sorted(reasons.items())),
        extracted_instances=len(instances),
        unique_mature_sequences=len(curated),
        duplicate_instances=len(instances) - len(curated),
        min_mature_length=min(lengths) if lengths else None,
        max_mature_length=max(lengths) if lengths else None,
        mean_mature_length=sum(lengths) / len(lengths) if lengths else None,
        mature_sequence_sha256=mature_sha,
        has_exact_mature_positives=bool(curated),
    )
    result = MaturePeptideCurationResult(sequences=curated, manifest=manifest)

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    csv_temporary = output_csv.with_suffix(output_csv.suffix + ".tmp")
    fieldnames = [
        "mature_id",
        "sequence",
        "source_accessions",
        "descriptions",
        "instance_count",
    ]
    with csv_temporary.open("w", encoding="utf-8", newline="") as output:
        writer = DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()
        for item in curated:
            row = item.model_dump()
            row["source_accessions"] = "|".join(item.source_accessions)
            row["descriptions"] = "|".join(item.descriptions)
            writer.writerow(row)
    csv_temporary.replace(output_csv)

    output_manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest_temporary = output_manifest.with_suffix(output_manifest.suffix + ".tmp")
    manifest_temporary.write_text(
        manifest.model_dump_json(indent=2) + "\n", encoding="utf-8"
    )
    manifest_temporary.replace(output_manifest)
    return result

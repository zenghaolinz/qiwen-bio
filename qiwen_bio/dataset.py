from __future__ import annotations

import json
import random
from csv import DictReader, DictWriter
from dataclasses import dataclass, field
from hashlib import sha256
from io import StringIO
from pathlib import Path
from typing import Literal

import httpx
from pydantic import BaseModel

from qiwen_bio.models import AMINO_ACIDS


Split = Literal["train", "validation", "test"]


@dataclass(frozen=True)
class RawAmpRecord:
    record_id: str
    sequence: str
    label: int
    source: str
    source_url: str | None = None


class AmpConflict(BaseModel):
    sequence: str
    source_ids: list[str]
    labels: list[int]


class PreparedAmpSample(BaseModel):
    sample_id: str
    sequence: str
    label: int
    source_ids: list[str]
    sources: list[str]
    cluster_id: str
    split: Split


class AmpDatasetManifest(BaseModel):
    input_records: int
    unique_sequences: int
    duplicate_records: int
    conflicting_sequences: int
    retained_samples: int
    identity_threshold: float
    seed: int
    cluster_count: int
    split_counts: dict[str, int]
    label_counts: dict[str, int]
    dataset_sha256: str
    positive_policy: str
    negative_policy: str
    source_url: str | None = None
    license_url: str | None = None
    positive_query: str | None = None
    negative_query: str | None = None


class PreparedAmpDataset(BaseModel):
    samples: list[PreparedAmpSample]
    conflicts: list[AmpConflict]
    manifest: AmpDatasetManifest


class UniProtAmpDatasetClient:
    endpoint = "https://rest.uniprot.org/uniprotkb/search"
    fields = "accession,id,sequence,length,keyword"

    @staticmethod
    def query_for_label(label: int) -> str:
        label_query = "(keyword:KW-0929)" if label == 1 else "NOT (keyword:KW-0929)"
        return f"(reviewed:true) AND {label_query} AND (length:[5 TO 200])"

    def __init__(self, http_client: httpx.Client | None = None) -> None:
        self.http_client = http_client or httpx.Client(timeout=30.0, follow_redirects=True)

    def fetch_records(self, *, label: int, limit: int, seed: int) -> list[RawAmpRecord]:
        if label not in (0, 1):
            raise ValueError("AMP labels must be 0 or 1")
        if limit < 1:
            raise ValueError("limit must be positive")
        query = self.query_for_label(label)
        response = self.http_client.get(
            self.endpoint,
            params={
                "query": query,
                "format": "tsv",
                "fields": self.fields,
                "size": min(max(limit * 5, limit), 500),
            },
        )
        response.raise_for_status()
        source = "UniProtKB reviewed AMP keyword" if label == 1 else "UniProtKB reviewed proxy-negative"
        records = [
            RawAmpRecord(
                record_id=row["Entry"],
                sequence=row["Sequence"],
                label=label,
                source=source,
                source_url=f"https://www.uniprot.org/uniprotkb/{row['Entry']}/entry",
            )
            for row in DictReader(StringIO(response.text), delimiter="\t")
            if row.get("Entry") and row.get("Sequence")
        ]
        if len(records) <= limit:
            return records
        return random.Random(seed).sample(records, limit)

    def prepare_dataset(
        self,
        *,
        positive_limit: int,
        negative_limit: int,
        seed: int,
        identity_threshold: float = 0.8,
    ) -> PreparedAmpDataset:
        positives = self.fetch_records(label=1, limit=positive_limit, seed=seed)
        negatives = self.fetch_records(label=0, limit=negative_limit, seed=seed + 1)
        prepared = prepare_amp_dataset(
            positives + negatives,
            identity_threshold=identity_threshold,
            seed=seed,
        )
        prepared.manifest.source_url = self.endpoint
        prepared.manifest.license_url = "https://www.uniprot.org/help/license"
        prepared.manifest.positive_query = self.query_for_label(1)
        prepared.manifest.negative_query = self.query_for_label(0)
        return prepared


def _normalize_sequence(sequence: str) -> str:
    normalized = "".join(sequence.split()).upper()
    invalid = sorted(set(normalized) - AMINO_ACIDS)
    if not normalized or invalid:
        raise ValueError("sequence must contain only canonical amino acids")
    return normalized


def global_sequence_identity(first: str, second: str) -> float:
    """Return matches/alignment length for a deterministic global alignment."""
    left = _normalize_sequence(first)
    right = _normalize_sequence(second)
    # Each cell stores alignment score, matches, and aligned length.
    previous = [(-2 * column, 0, column) for column in range(len(right) + 1)]
    for row, left_residue in enumerate(left, start=1):
        current = [(-2 * row, 0, row)]
        for column, right_residue in enumerate(right, start=1):
            diagonal = previous[column - 1]
            candidates = [
                (
                    diagonal[0] + (2 if left_residue == right_residue else -1),
                    diagonal[1] + int(left_residue == right_residue),
                    diagonal[2] + 1,
                ),
                (previous[column][0] - 2, previous[column][1], previous[column][2] + 1),
                (current[column - 1][0] - 2, current[column - 1][1], current[column - 1][2] + 1),
            ]
            current.append(max(candidates, key=lambda value: (value[0], value[1], -value[2])))
        previous = current
    _, matches, aligned_length = previous[-1]
    return round(matches / aligned_length, 10)


class _UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[max(left_root, right_root)] = min(left_root, right_root)


@dataclass
class _MergedRecord:
    sequence: str
    label: int
    source_ids: list[str] = field(default_factory=list)
    sources: list[str] = field(default_factory=list)


def _assign_cluster_splits(
    clusters: list[list[int]], labels: list[int], seed: int
) -> dict[int, Split]:
    rng = random.Random(seed)
    split_order: list[Split] = ["train", "validation", "test"]
    fractions = {"train": 0.7, "validation": 0.15, "test": 0.15}
    targets = {split: fraction * len(labels) for split, fraction in fractions.items()}
    totals_by_label = {label: labels.count(label) for label in (0, 1)}
    label_targets = {
        split: {label: fraction * totals_by_label[label] for label in (0, 1)}
        for split, fraction in fractions.items()
    }
    counts = {split: 0 for split in split_order}
    label_counts = {split: {0: 0, 1: 0} for split in split_order}
    assignments: dict[int, Split] = {}

    pure_clusters: dict[int, list[list[int]]] = {0: [], 1: []}
    mixed_clusters: list[list[int]] = []
    for cluster in clusters:
        cluster_labels = {labels[index] for index in cluster}
        if len(cluster_labels) == 1:
            pure_clusters[cluster_labels.pop()].append(cluster)
        else:
            mixed_clusters.append(cluster)
    ordered: list[list[int]] = []
    for label in (0, 1):
        rng.shuffle(pure_clusters[label])
        pure_clusters[label].sort(key=len, reverse=True)
        ordered.extend(pure_clusters[label])
    remaining_pure = {label: len(pure_clusters[label]) for label in (0, 1)}
    rng.shuffle(mixed_clusters)
    mixed_clusters.sort(key=len, reverse=True)
    ordered.extend(mixed_clusters)

    for cluster in ordered:
        additions = {label: sum(labels[index] == label for index in cluster) for label in (0, 1)}

        def allocation_cost(candidate: Split) -> tuple[float, float, int]:
            label_delta = sum(
                (label_counts[candidate][label] + additions[label] - label_targets[candidate][label]) ** 2
                - (label_counts[candidate][label] - label_targets[candidate][label]) ** 2
                for label in (0, 1)
            )
            total_delta = (
                (counts[candidate] + len(cluster) - targets[candidate]) ** 2
                - (counts[candidate] - targets[candidate]) ** 2
            )
            return label_delta, total_delta, split_order.index(candidate)

        present_labels = [label for label, addition in additions.items() if addition]
        candidates = split_order
        if len(present_labels) == 1:
            pure_label = present_labels[0]
            uncovered = [split for split in split_order if label_counts[split][pure_label] == 0]
            if remaining_pure[pure_label] <= len(uncovered):
                candidates = uncovered
            remaining_pure[pure_label] -= 1
        split = min(candidates, key=allocation_cost)
        counts[split] += len(cluster)
        for label in (0, 1):
            label_counts[split][label] += additions[label]
        for sample_index in cluster:
            assignments[sample_index] = split
    return assignments


def prepare_amp_dataset(
    records: list[RawAmpRecord], *, identity_threshold: float = 0.8, seed: int = 42
) -> PreparedAmpDataset:
    if not 0 < identity_threshold <= 1:
        raise ValueError("identity_threshold must be in (0, 1]")

    grouped: dict[str, list[RawAmpRecord]] = {}
    for record in records:
        if record.label not in (0, 1):
            raise ValueError("AMP labels must be 0 or 1")
        grouped.setdefault(_normalize_sequence(record.sequence), []).append(record)

    conflicts: list[AmpConflict] = []
    merged: list[_MergedRecord] = []
    duplicate_records = 0
    for sequence, sequence_records in grouped.items():
        labels = sorted({record.label for record in sequence_records})
        if len(labels) > 1:
            conflicts.append(
                AmpConflict(
                    sequence=sequence,
                    source_ids=[record.record_id for record in sequence_records],
                    labels=labels,
                )
            )
            continue
        duplicate_records += len(sequence_records) - 1
        merged.append(
            _MergedRecord(
                sequence=sequence,
                label=labels[0],
                source_ids=[record.record_id for record in sequence_records],
                sources=list(dict.fromkeys(record.source for record in sequence_records)),
            )
        )

    union_find = _UnionFind(len(merged))
    for left_index, left in enumerate(merged):
        for right_index in range(left_index + 1, len(merged)):
            right = merged[right_index]
            length_ratio = min(len(left.sequence), len(right.sequence)) / max(len(left.sequence), len(right.sequence))
            if length_ratio >= identity_threshold and global_sequence_identity(left.sequence, right.sequence) >= identity_threshold:
                union_find.union(left_index, right_index)

    cluster_members: dict[int, list[int]] = {}
    for index in range(len(merged)):
        cluster_members.setdefault(union_find.find(index), []).append(index)
    clusters = sorted(cluster_members.values(), key=lambda members: members[0])
    cluster_ids = {
        index: f"cluster-{cluster_number:04d}"
        for cluster_number, members in enumerate(clusters, start=1)
        for index in members
    }
    splits = _assign_cluster_splits(clusters, [record.label for record in merged], seed)
    samples = [
        PreparedAmpSample(
            sample_id=record.source_ids[0],
            sequence=record.sequence,
            label=record.label,
            source_ids=record.source_ids,
            sources=record.sources,
            cluster_id=cluster_ids[index],
            split=splits[index],
        )
        for index, record in enumerate(merged)
    ]
    digest_payload = [sample.model_dump(mode="json") for sample in samples]
    dataset_sha = sha256(json.dumps(digest_payload, sort_keys=True).encode()).hexdigest()
    split_counts = {split: sum(sample.split == split for sample in samples) for split in ("train", "validation", "test")}
    label_counts = {str(label): sum(sample.label == label for sample in samples) for label in (0, 1)}
    manifest = AmpDatasetManifest(
        input_records=len(records),
        unique_sequences=len(grouped),
        duplicate_records=duplicate_records,
        conflicting_sequences=len(conflicts),
        retained_samples=len(samples),
        identity_threshold=identity_threshold,
        seed=seed,
        cluster_count=len(clusters),
        split_counts=split_counts,
        label_counts=label_counts,
        dataset_sha256=dataset_sha,
        positive_policy="UniProtKB reviewed entries annotated with antimicrobial keyword KW-0929",
        negative_policy=(
            "proxy-negative: UniProtKB reviewed entries without KW-0929; absence of annotation "
            "does not prove absence of antimicrobial activity"
        ),
    )
    return PreparedAmpDataset(samples=samples, conflicts=conflicts, manifest=manifest)


def write_prepared_dataset(
    prepared: PreparedAmpDataset, *, csv_path: Path, manifest_path: Path
) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id",
        "sequence",
        "label",
        "source_ids",
        "sources",
        "cluster_id",
        "split",
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as output:
        writer = DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()
        for sample in prepared.samples:
            row = sample.model_dump()
            row["source_ids"] = "|".join(sample.source_ids)
            row["sources"] = "|".join(sample.sources)
            writer.writerow(row)
    audit = {
        **prepared.manifest.model_dump(mode="json"),
        "conflicts": [conflict.model_dump(mode="json") for conflict in prepared.conflicts],
    }
    manifest_path.write_text(
        json.dumps(audit, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

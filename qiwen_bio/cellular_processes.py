from __future__ import annotations

from collections import Counter
from typing import Literal

from pydantic import BaseModel, Field

from qiwen_bio.kegg import KeggPathwayAnnotation
from qiwen_bio.stringdb import EvidenceGraph
from qiwen_bio.uniprot import UniProtAnnotation


class ProcessSupport(BaseModel):
    source_name: Literal["UniProt GO", "KEGG", "STRING enrichment"]
    evidence_type: Literal[
        "database_annotation", "database_membership", "enrichment_statistic"
    ]
    source_url: str
    detail: str
    fdr: float | None = None


class CellularProcess(BaseModel):
    canonical_id: str
    label: str
    process_type: Literal["biological_process", "pathway"]
    supports: list[ProcessSupport]


class CellularProcessEvidence(BaseModel):
    protein_accession: str
    processes: list[CellularProcess]
    counts_by_source: dict[str, int]
    literature_context: list[str]
    phenotype_hypotheses: list[str] = Field(default_factory=list)
    interpretation_boundary: str = (
        "These records establish database annotation, membership, or enrichment only. "
        "A process record does not establish activity, effect direction, mechanism, causality, "
        "or a cellular or organismal phenotype."
    )


def build_cellular_process_evidence(
    annotation: UniProtAnnotation,
    kegg: KeggPathwayAnnotation | None,
    graph: EvidenceGraph | None,
) -> CellularProcessEvidence:
    merged: dict[str, CellularProcess] = {}

    def add_support(
        canonical_id: str,
        label: str,
        process_type: Literal["biological_process", "pathway"],
        support: ProcessSupport,
    ) -> None:
        process = merged.get(canonical_id)
        if process is None:
            process = CellularProcess(
                canonical_id=canonical_id,
                label=label,
                process_type=process_type,
                supports=[],
            )
            merged[canonical_id] = process
        identity = (support.source_name, support.evidence_type, support.source_url)
        if identity not in {
            (item.source_name, item.evidence_type, item.source_url)
            for item in process.supports
        }:
            process.supports.append(support)

    for term in annotation.go_terms:
        if term.aspect != "biological_process":
            continue
        add_support(
            term.id,
            term.name,
            "biological_process",
            ProcessSupport(
                source_name="UniProt GO",
                evidence_type="database_annotation",
                source_url=annotation.source_url,
                detail=f"UniProt annotates {annotation.accession} to {term.id}.",
            ),
        )

    if kegg:
        for pathway in kegg.pathways:
            add_support(
                pathway.pathway_id,
                pathway.name,
                "pathway",
                ProcessSupport(
                    source_name="KEGG",
                    evidence_type="database_membership",
                    source_url=pathway.source_url,
                    detail=(
                        f"KEGG links {', '.join(kegg.gene_ids)} to "
                        f"{pathway.pathway_id}."
                    ),
                ),
            )

    if graph:
        seed_node_ids = {
            node.id
            for node in graph.nodes
            if node.type == "protein" and node.label.upper() == graph.seed.upper()
        }
        seed_term_ids = {
            edge.target
            for edge in graph.edges
            if edge.type == "annotated_to" and edge.source in seed_node_ids
        }
        for node in graph.nodes:
            if node.id not in seed_term_ids or node.type not in ("process", "pathway"):
                continue
            add_support(
                node.external_id,
                node.label,
                "biological_process" if node.type == "process" else "pathway",
                ProcessSupport(
                    source_name="STRING enrichment",
                    evidence_type="enrichment_statistic",
                    source_url=node.source_url,
                    detail=(
                        f"STRING enrichment links seed {graph.seed} to "
                        f"{node.external_id}."
                    ),
                    fdr=node.fdr,
                ),
            )

    source_priority = {"UniProt GO": 0, "KEGG": 1, "STRING enrichment": 2}
    processes = sorted(merged.values(), key=lambda item: item.canonical_id.lower())
    for process in processes:
        process.supports.sort(key=lambda item: source_priority[item.source_name])
    counts = Counter(
        support.source_name
        for process in processes
        for support in process.supports
    )
    context_order = sorted(
        processes,
        key=lambda process: (
            not any(
                support.evidence_type != "enrichment_statistic"
                for support in process.supports
            ),
            process.canonical_id.lower(),
        ),
    )
    return CellularProcessEvidence(
        protein_accession=annotation.accession,
        processes=processes,
        counts_by_source=dict(sorted(counts.items())),
        literature_context=[process.label for process in context_order],
    )

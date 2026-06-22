"""Cross-scale reasoning chain (stage 3G).

This module assembles the already-retrieved evidence layers into a single
structured chain that connects a mutation (or, in its absence, the protein
itself) to structure, function, pathway, and phenotype evidence. It is a
**synthesis layer**: it performs no new network requests and introduces no new
biological claims of its own.

Biological-correctness guardrails (enforced here):

* The chain is a **pure evidence chain**. Each step states deterministic,
  source-bound facts first. A hypothesis sentence ("可能影响/可能相关/提示") is
  emitted only when explicit upstream evidence supports it, and every
  hypothesis carries an uncertainty clause. No step asserts causality,
  direction, mechanism, activity, or a phenotype outcome.
* For a mutation-impact chain, the wild-type residue is validated against the
  UniProt sequence. A mismatch suppresses downstream "可能影响" hypotheses:
  facts remain, but no functional-impact inference is drawn from a mutation
  that does not match the annotated sequence.
* Per ADR-0012/0013: process association and abstract-level literature do not
  establish activity, direction, mechanism, causality, or phenotype.
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel

from qiwen_bio.alphafold import AlphaFoldAnalysis
from qiwen_bio.cellular_processes import CellularProcessEvidence
from qiwen_bio.interpro import DomainAnnotation
from qiwen_bio.kegg import KeggPathwayAnnotation
from qiwen_bio.models import AnalysisResponse
from qiwen_bio.mutation import parse_mutation, wild_type_matches_sequence
from qiwen_bio.phenotype_literature import PhenotypeLiteratureEvidence
from qiwen_bio.stringdb import EvidenceGraph
from qiwen_bio.structure_features import (
    StructureFeatureSummary,
    build_structure_feature_summary,
)
from qiwen_bio.uniprot import UniProtAnnotation


ChainType = Literal["mutation_impact", "protein_function"]
StepConfidence = Literal["high", "medium", "low", "insufficient"]


class ChainStep(BaseModel):
    step_id: str
    title: str
    evidence_facts: list[str]
    hypothesis: str | None = None
    evidence_sources: list[str]
    confidence: StepConfidence
    uncertainty: str
    available: bool


class ReasoningChain(BaseModel):
    chain_type: ChainType
    gene: str
    mutation: str | None
    steps: list[ChainStep]
    structure_summary: StructureFeatureSummary | None = None
    missing_layers: list[str]
    summary: str
    boundary: str = (
        "This chain organizes retrieved and calculated evidence. It does not establish "
        "activity, effect direction, mechanism, causality, pathogenicity, or a cellular "
        "or organismal phenotype. Hypothesis sentences are explicitly labelled and require "
        "full-text appraisal and experimental validation before use."
    )


def _confidence_from_plddt(plddt: float | None) -> StepConfidence:
    """Map a pLDDT value to a reasoning-STEP confidence level.

    This is intentionally a DIFFERENT vocabulary from the AlphaFold pLDDT band
    (``very_high``/``confident``/``low``/``very_low`` in
    :func:`qiwen_bio.structure_features.plddt_confidence_band`). Here the
    output is the chain-step confidence ``high``/``medium``/``low``/
    ``insufficient``, which describes how much weight a reasoning step can
    place on the structural evidence — not the AlphaFold model band itself.
    The numeric cutoffs coincide but the semantics are distinct, so this is
    not a duplicate of the band function.
    """
    if plddt is None:
        return "insufficient"
    if plddt >= 90:
        return "high"
    if plddt >= 70:
        return "medium"
    if plddt >= 50:
        return "low"
    return "low"


def _direct_processes(cellular_processes: CellularProcessEvidence) -> list:
    """Processes with at least one non-enrichment support (direct annotation/membership)."""
    return [
        process
        for process in cellular_processes.processes
        if any(
            support.evidence_type in ("database_annotation", "database_membership")
            for support in process.supports
        )
    ]


def build_mutation_step(
    annotation: UniProtAnnotation,
    structure: AlphaFoldAnalysis | None,
    domains: DomainAnnotation | None,
    mutation: str | None,
) -> ChainStep:
    """Build the mutation step. Records a mismatch fact when the wild-type
    residue does not match the UniProt sequence, which downstream steps use to
    suppress functional-impact hypotheses."""
    parsed = parse_mutation(mutation)
    if parsed.status != "parsed" or parsed.wild_type is None or parsed.position is None or parsed.mutant is None:
        return ChainStep(
            step_id="mutation",
            title="Mutation input",
            evidence_facts=[f"Invalid or malformed mutation input: {mutation!r}."],
            evidence_sources=[],
            confidence="insufficient",
            uncertainty="No valid mutation could be parsed; downstream impact hypotheses are suppressed.",
            available=False,
        )

    wild_type, position, mutant = parsed.wild_type, parsed.position, parsed.mutant
    facts: list[str] = [f"User-supplied mutation: {wild_type}{position}{mutant}."]
    sources: list[str] = ["User input"]

    matches = wild_type_matches_sequence(annotation.sequence, wild_type, position)
    if matches:
        facts.append(
            f"Wild-type {wild_type} at position {position} matches the UniProt sequence."
        )
    else:
        actual = (
            annotation.sequence[position - 1]
            if 1 <= position <= len(annotation.sequence)
            else "absent"
        )
        facts.append(
            f"Wild-type mismatch: UniProt sequence has {actual} at position {position}, "
            f"not {wild_type}. The mutation may use different numbering or refer to "
            f"a different isoform; impact hypotheses are suppressed."
        )
    sources.append(annotation.source_url)

    domain_overlap_names: list[str] = []
    if domains and domains.mutation_overlaps:
        domain_overlap_names = [entry.accession for entry in domains.mutation_overlaps]
        facts.append(
            f"Mutation position overlaps {len(domains.mutation_overlaps)} annotated "
            f"domain entry/entries: {', '.join(domain_overlap_names)}."
        )
        sources.extend(entry.source_url for entry in domains.mutation_overlaps)
    elif domains:
        facts.append("Mutation position does not overlap any annotated InterPro/Pfam domain entry.")

    if structure and structure.mutation_site:
        site = structure.mutation_site
        facts.append(
            f"AlphaFold mutation site pLDDT: {site.plddt:.2f} ({site.confidence}); "
            f"{len(structure.mutation_neighborhood)} CA neighbours within 8 A."
        )
        sources.append(structure.structure_url)

    has_overlap = bool(domain_overlap_names)
    hypothesis: str | None = None
    if matches and has_overlap and structure and structure.mutation_site:
        hypothesis = (
            f"假设：该突变位点 {wild_type}{position}{mutant} 位于已注释结构域 "
            f"{', '.join(domain_overlap_names)} 内，因此可作为后续功能影响评估的候选位点。"
            f"该判断仅基于坐标重叠和预测结构上下文，不代表结构扰动、稳定性变化、功能改变或致病性结论；"
            f"需要文献证据或实验进一步验证。"
        )

    confidence: StepConfidence = "low"
    if not matches:
        confidence = "insufficient"
    elif has_overlap and structure and structure.mutation_site:
        confidence = _confidence_from_plddt(structure.mutation_site.plddt)

    return ChainStep(
        step_id="mutation",
        title="Mutation and structural context",
        evidence_facts=facts,
        hypothesis=hypothesis,
        evidence_sources=sources,
        confidence=confidence,
        uncertainty=(
            "Positional overlap is a coordinate containment observation, not a "
            "functional-effect or pathogenicity prediction."
        ),
        available=True,
    )


def build_structure_step(summary: StructureFeatureSummary) -> ChainStep:
    """Build the structure step from a unified StructureFeatureSummary.

    The step reports structural facts only; it never emits a functional-impact
    hypothesis (structure is geometric, not functional). When the mutation site
    lies in a low-confidence predicted region the uncertainty clause is
    strengthened to flag that structural interpretation is limited.
    """
    if not summary.has_structure:
        return ChainStep(
            step_id="structure",
            title="Predicted structure",
            evidence_facts=list(summary.evidence_facts),
            evidence_sources=[],
            confidence="insufficient",
            uncertainty="Structure layer unavailable; geometric context cannot be assessed.",
            available=False,
        )
    # Prefer the mutation-site pLDDT for step confidence when available,
    # otherwise fall back to the whole-model mean.
    plddt_for_confidence = summary.mutation_site_plddt or summary.mean_plddt
    confidence = _confidence_from_plddt(plddt_for_confidence)
    if summary.low_confidence_region:
        uncertainty = (
            "The mutation site lies in a low-confidence predicted region; structural "
            "interpretation is limited. pLDDT and CA proximity describe model confidence "
            "and geometry, not pathogenicity, stability, or functional effect."
        )
    else:
        uncertainty = (
            "pLDDT and CA proximity describe model confidence and geometry, not "
            "pathogenicity, stability, or functional effect."
        )
    return ChainStep(
        step_id="structure",
        title="Predicted structure",
        evidence_facts=list(summary.evidence_facts),
        hypothesis=None,
        evidence_sources=[summary.source_url] if summary.source_url else [summary.source],
        confidence=confidence,
        uncertainty=uncertainty,
        available=True,
    )


def build_function_step(
    annotation: UniProtAnnotation,
    domains: DomainAnnotation | None,
    mutation_step: ChainStep | None,
) -> ChainStep:
    facts: list[str] = []
    sources: list[str] = [annotation.source_url]
    if annotation.functions:
        facts.append(f"UniProt FUNCTION: {annotation.functions[0]}")
    else:
        facts.append("No UniProt FUNCTION annotation available.")
    bp_terms = [
        term for term in annotation.go_terms if term.aspect == "biological_process"
    ]
    facts.append(f"UniProt GO biological_process terms: {len(bp_terms)}.")
    if domains:
        names = ", ".join(entry.name for entry in domains.entries[:6])
        facts.append(f"InterPro/Pfam domain entries: {domains.entry_count} ({names}).")
        sources.extend(entry.source_url for entry in domains.entries[:6])
    else:
        facts.append("No InterPro/Pfam domain annotation available.")

    # A functional-impact hypothesis requires a matching mutation with a domain overlap.
    hypothesis: str | None = None
    if (
        mutation_step
        and mutation_step.available
        and mutation_step.hypothesis is not None
        and domains
        and domains.mutation_overlaps
    ):
        overlap_names = ", ".join(entry.accession for entry in domains.mutation_overlaps)
        hypothesis = (
            f"候选解释：该突变位点与功能注释结构域 {overlap_names} 存在坐标重叠，"
            f"因此可优先纳入功能影响评估。该判断不是功能改变结论，"
            f"需要文献证据或实验进一步验证。"
        )

    confidence: StepConfidence = "medium" if annotation.functions else "low"
    return ChainStep(
        step_id="function",
        title="Function and domain annotation",
        evidence_facts=facts,
        hypothesis=hypothesis,
        evidence_sources=sources,
        confidence=confidence,
        uncertainty=(
            "Database annotations describe known function; they do not predict how a "
            "specific mutation alters that function."
        ),
        available=True,
    )


def build_pathway_step(
    cellular_processes: CellularProcessEvidence,
    kegg: KeggPathwayAnnotation | None,
    graph: EvidenceGraph | None,
    function_step: ChainStep,
) -> ChainStep:
    facts: list[str] = []
    sources: list[str] = []
    direct = _direct_processes(cellular_processes)
    if direct:
        ids = ", ".join(process.canonical_id for process in direct[:8])
        facts.append(f"Directly supported processes/pathways: {len(direct)} ({ids}).")
        sources.extend(
            support.source_url
            for process in direct
            for support in process.supports
            if support.evidence_type in ("database_annotation", "database_membership")
        )
    else:
        facts.append("No directly supported process/pathway evidence (UniProt GO or direct KEGG).")

    if kegg:
        facts.append(
            f"Direct KEGG memberships: {kegg.pathway_count} pathway(s) for gene "
            f"{', '.join(kegg.gene_ids)}."
        )
        sources.append("https://www.kegg.jp/")

    if graph:
        interactions = sum(1 for edge in graph.edges if edge.type == "interacts_with")
        terms = [node for node in graph.nodes if node.type != "protein"]
        facts.append(f"STRING interaction edges: {interactions}; enrichment terms: {len(terms)}.")
        sources.append("https://string-db.org/")

    # A pathway-association hypothesis requires a functional-impact hypothesis upstream.
    hypothesis: str | None = None
    if function_step.hypothesis is not None and direct:
        hypothesis = (
            f"假设：如果该蛋白功能受影响，其参与的通路（如 {direct[0].label}）可能相关。"
            f"这是基于功能假设和数据库成员关系的二级假设，不是通路激活、方向性或表型因果结论；"
            f"需实验验证。"
        )

    confidence: StepConfidence = "medium" if direct else "low"
    return ChainStep(
        step_id="pathway",
        title="Pathway and interaction context",
        evidence_facts=facts,
        hypothesis=hypothesis,
        evidence_sources=sources,
        confidence=confidence,
        uncertainty=(
            "Pathway membership and interaction scores are database associations, not "
            "evidence of pathway activation, directionality, or phenotype causality."
        ),
        available=True,
    )


def build_phenotype_step(
    phenotype_literature: PhenotypeLiteratureEvidence | None,
    pathway_step: ChainStep,
) -> ChainStep:
    facts: list[str] = []
    sources: list[str] = []
    if phenotype_literature:
        counts = ", ".join(
            f"{level}={count}"
            for level, count in phenotype_literature.counts_by_level.items()
        ) or "none"
        facts.append(
            f"Phenotype literature claim-support: {counts}; "
            f"{len(phenotype_literature.process_links)} process-link classification(s)."
        )
        pmids = [link.pmid for link in phenotype_literature.process_links[:5]]
        if pmids:
            facts.append(f"Classified PMIDs: {', '.join(pmids)}.")
            sources.extend(
                f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" for pmid in pmids
            )
    else:
        facts.append("No phenotype literature layer available.")

    hypothesis: str | None = None
    if phenotype_literature and phenotype_literature.hypotheses:
        hypothesis = phenotype_literature.hypotheses[0]
    elif phenotype_literature:
        facts.append(
            "No phenotype hypothesis was emitted: no directly supported process had an "
            "article classified as 'supports' at the abstract-metadata boundary."
        )

    confidence: StepConfidence = "low"
    if phenotype_literature and phenotype_literature.counts_by_level.get("supports", 0) > 0:
        confidence = "medium"

    return ChainStep(
        step_id="phenotype",
        title="Phenotype literature evidence",
        evidence_facts=facts,
        hypothesis=hypothesis,
        evidence_sources=sources,
        confidence=confidence,
        uncertainty=(
            "Abstract-metadata classification is not full-text appraisal, claim-level "
            "grading, or causal inference."
        ),
        available=phenotype_literature is not None,
    )


def build_reasoning_chain(
    annotation: UniProtAnnotation,
    analysis: AnalysisResponse,
    structure: AlphaFoldAnalysis | None,
    domains: DomainAnnotation | None,
    kegg: KeggPathwayAnnotation | None,
    graph: EvidenceGraph | None,
    cellular_processes: CellularProcessEvidence,
    phenotype_literature: PhenotypeLiteratureEvidence | None,
    mutation: str | None,
) -> ReasoningChain:
    """Assemble the cross-scale reasoning chain from already-retrieved layers."""
    gene = annotation.gene_names[0] if annotation.gene_names else annotation.accession
    # chain_type reflects user intent: any supplied mutation string (even a
    # malformed one) triggers the mutation-impact chain type, so the chain
    # surfaces the invalid input rather than silently switching types.
    chain_type: ChainType = "mutation_impact" if mutation else "protein_function"

    steps: list[ChainStep] = []
    missing_layers: list[str] = []

    mutation_step: ChainStep | None = None
    if chain_type == "mutation_impact":
        mutation_step = build_mutation_step(annotation, structure, domains, mutation)
        steps.append(mutation_step)

    structure_summary = build_structure_feature_summary(
        annotation=annotation,
        structure=structure,
        domains=domains,
        mutation=mutation,
    )
    structure_step = build_structure_step(structure_summary)
    if not structure_step.available:
        missing_layers.append("structure")
    steps.append(structure_step)

    function_step = build_function_step(annotation, domains, mutation_step)
    steps.append(function_step)

    pathway_step = build_pathway_step(cellular_processes, kegg, graph, function_step)
    steps.append(pathway_step)
    if not cellular_processes.processes:
        missing_layers.append("cellular_processes")

    phenotype_step = build_phenotype_step(phenotype_literature, pathway_step)
    if not phenotype_step.available:
        missing_layers.append("phenotype_literature")
    steps.append(phenotype_step)

    if kegg is None:
        missing_layers.append("kegg")
    if graph is None:
        missing_layers.append("string_graph")

    if chain_type == "mutation_impact":
        summary = (
            f"Mutation-impact reasoning chain for {gene} {mutation}: connects the mutation "
            f"to structural context, function, pathway, and phenotype literature evidence."
        )
    else:
        summary = (
            f"Protein-function reasoning chain for {gene}: connects domain/function "
            f"annotation, pathway, and phenotype literature evidence (no mutation supplied)."
        )

    return ReasoningChain(
        chain_type=chain_type,
        gene=gene,
        mutation=mutation if chain_type == "mutation_impact" else None,
        steps=steps,
        structure_summary=structure_summary,
        missing_layers=missing_layers,
        summary=summary,
    )

"""Unified structure feature summary (Stage 2 convergence).

This module is the single source of truth for the pLDDT confidence-band
vocabulary and for a reusable :class:`StructureFeatureSummary` that aggregates
AlphaFold structure, CA geometry, and domain-overlap evidence into one object.

Biological-correctness guardrails (enforced here and in callers):

* The summary reports **structural facts only**. It never emits functional-
  effect, pathogenicity, stability, or phenotype statements.
* pLDDT is AlphaFold **local model confidence**, not pathogenicity or
  functional-effect confidence.
* CA contacts within 8 A are **geometric proximity**, not confirmed
  biochemical interactions.
* Domain overlap means the mutation position lies inside an annotated
  fragment; it is a coordinate-containment observation, not a
  functional-impact prediction.
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel

from qiwen_bio.alphafold import AlphaFoldAnalysis
from qiwen_bio.interpro import DomainAnnotation
from qiwen_bio.mutation import parse_mutation
from qiwen_bio.uniprot import UniProtAnnotation


# pLDDT band thresholds (AlphaFold convention). Single source of truth —
# alphafold.confidence_label delegates here, and reasoning_chain consumes the
# band string off the summary rather than re-deriving it.
PLDDT_BANDS = {"very_high": 90.0, "confident": 70.0, "low": 50.0}

PlddtBand = Literal["very_high", "confident", "low", "very_low", "unavailable"]
StructureSourceType = Literal["alphafold_predicted", "pdb_experimental", "none"]


def plddt_confidence_band(plddt: float | None) -> PlddtBand:
    """Map a pLDDT value to its AlphaFold confidence band.

    ``None`` maps to ``"unavailable"``. The numeric thresholds follow the
    AlphaFold convention (>=90 very_high, >=70 confident, >=50 low, else
    very_low). These bands describe local model confidence only.
    """
    if plddt is None:
        return "unavailable"
    if plddt >= PLDDT_BANDS["very_high"]:
        return "very_high"
    if plddt >= PLDDT_BANDS["confident"]:
        return "confident"
    if plddt >= PLDDT_BANDS["low"]:
        return "low"
    return "very_low"


class StructureFeatureSummary(BaseModel):
    """Aggregated, reusable structure-evidence summary.

    Built from an :class:`AlphaFoldAnalysis` and an optional
    :class:`DomainAnnotation` via :func:`build_structure_feature_summary`.
    The summary is consumed by the reasoning chain, the comprehensive report,
    and the structure-feature export so structure evidence is described
    consistently everywhere.
    """

    source: str
    accession: str
    mutation: str | None
    mutation_position: int | None
    has_structure: bool
    structure_source_type: StructureSourceType
    mean_plddt: float | None
    mutation_site_plddt: float | None
    mutation_site_confidence_band: PlddtBand
    low_confidence_region: bool
    contact_count_8a: int | None
    neighbor_count_8a: int | None
    nearest_neighbors: list[str]
    domain_overlap: bool | None
    overlapping_domains: list[str]
    interpretation_limits: list[str]
    evidence_facts: list[str]


_CANONICAL_LIMITS = [
    "AlphaFold is a predicted structure, not an experimental determination.",
    "pLDDT is local model confidence, not pathogenicity or functional-effect confidence.",
    "CA contacts within 8 A are geometric proximity, not confirmed biochemical interactions.",
]

_LOW_CONFIDENCE_LIMIT = (
    "The mutation site lies in a low-confidence predicted region; structural "
    "interpretation is limited."
)

_NO_STRUCTURE_LIMIT = "No AlphaFold structure is available for this entry."

_DOMAIN_OVERLAP_LIMIT = (
    "Domain overlap is a coordinate-containment observation, not a functional-impact "
    "or pathogenicity prediction."
)


def _format_neighbors(structure: AlphaFoldAnalysis, limit: int = 6) -> list[str]:
    neighbors = sorted(structure.mutation_neighborhood, key=lambda n: n.distance_angstrom)
    formatted = []
    for neighbor in neighbors[:limit]:
        formatted.append(
            f"{neighbor.amino_acid}{neighbor.position} ({neighbor.distance_angstrom:.1f} A)"
        )
    return formatted


def _domain_names(domains: DomainAnnotation | None) -> list[str]:
    if not domains or not domains.mutation_overlaps:
        return []
    return [
        f"{entry.accession} ({entry.name})" for entry in domains.mutation_overlaps
    ]


def build_structure_feature_summary(
    annotation: UniProtAnnotation,
    structure: AlphaFoldAnalysis | None,
    domains: DomainAnnotation | None,
    mutation: str | None,
) -> StructureFeatureSummary:
    """Build a unified structure-evidence summary from existing layer outputs.

    Reuses the already-parsed :class:`AlphaFoldAnalysis` and
    :class:`DomainAnnotation` fields; performs no new parsing and no new
    network requests. Wild-type validation is intentionally NOT repeated here
    (the reasoning chain owns that); the summary reports structural facts only.
    """
    parsed = parse_mutation(mutation)
    mutation_position = parsed.position if parsed.status == "parsed" else None

    if structure is None:
        return StructureFeatureSummary(
            source="AlphaFold DB",
            accession=annotation.accession,
            mutation=mutation,
            mutation_position=mutation_position,
            has_structure=False,
            structure_source_type="none",
            mean_plddt=None,
            mutation_site_plddt=None,
            mutation_site_confidence_band="unavailable",
            low_confidence_region=False,
            contact_count_8a=None,
            neighbor_count_8a=None,
            nearest_neighbors=[],
            domain_overlap=None,
            overlapping_domains=[],
            interpretation_limits=[_NO_STRUCTURE_LIMIT],
            evidence_facts=["No AlphaFold structure is available for this entry."],
        )

    site_plddt = structure.mutation_site.plddt if structure.mutation_site else None
    band = plddt_confidence_band(site_plddt)
    low_confidence = site_plddt is not None and site_plddt < PLDDT_BANDS["low"]
    neighbors = _format_neighbors(structure)
    overlap_names = _domain_names(domains)
    has_overlap = bool(overlap_names)

    limits = list(_CANONICAL_LIMITS)
    if low_confidence:
        limits.append(_LOW_CONFIDENCE_LIMIT)
    if has_overlap:
        limits.append(_DOMAIN_OVERLAP_LIMIT)

    facts: list[str] = [
        f"AlphaFold predicted structure: {structure.residue_count} residues, "
        f"mean pLDDT {structure.mean_plddt:.2f}.",
        f"Non-local CA contacts (8 A threshold): {structure.contact_map.total_contacts}.",
    ]
    if site_plddt is not None and structure.mutation_site:
        site = structure.mutation_site
        facts.append(
            f"Mutation site {site.wild_type}{site.position}{site.mutant} maps to a predicted "
            f"structural region with pLDDT {site.plddt:.2f} ({band})."
        )
    if neighbors:
        facts.append(
            f"CA geometric neighbors within 8 A: {len(structure.mutation_neighborhood)} "
            f"(nearest: {', '.join(neighbors[:3])})."
        )
    if has_overlap:
        facts.append(
            f"Mutation position overlaps annotated domain entry/entries: "
            f"{', '.join(overlap_names)}."
        )

    return StructureFeatureSummary(
        source="AlphaFold DB",
        accession=annotation.accession,
        mutation=mutation,
        mutation_position=mutation_position,
        has_structure=True,
        structure_source_type="alphafold_predicted",
        mean_plddt=structure.mean_plddt,
        mutation_site_plddt=site_plddt,
        mutation_site_confidence_band=band,
        low_confidence_region=low_confidence,
        contact_count_8a=structure.contact_map.total_contacts,
        neighbor_count_8a=len(structure.mutation_neighborhood),
        nearest_neighbors=neighbors,
        domain_overlap=has_overlap if domains is not None else None,
        overlapping_domains=overlap_names,
        interpretation_limits=limits,
        evidence_facts=facts,
    )

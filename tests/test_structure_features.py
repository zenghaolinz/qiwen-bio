"""Tests for the unified structure feature summary (Stage 2 convergence)."""

import pytest

from qiwen_bio.alphafold import (
    AlphaFoldAnalysis,
    ConfidenceDistribution,
    ContactMap,
    MutationNeighbor,
    MutationSiteAnalysis,
    confidence_label,
)
from qiwen_bio.interpro import (
    DomainAnnotation,
    DomainEntry,
    DomainLocation,
)
from qiwen_bio.structure_features import (
    PLDDT_BANDS,
    StructureFeatureSummary,
    build_structure_feature_summary,
    plddt_confidence_band,
)
from qiwen_bio.uniprot import UniProtAnnotation


# ---------------------------------------------------------------------------
# pLDDT band consolidation
# ---------------------------------------------------------------------------


def test_plddt_confidence_band_thresholds() -> None:
    assert plddt_confidence_band(95.0) == "very_high"
    assert plddt_confidence_band(90.0) == "very_high"
    assert plddt_confidence_band(89.9) == "confident"
    assert plddt_confidence_band(70.0) == "confident"
    assert plddt_confidence_band(69.9) == "low"
    assert plddt_confidence_band(50.0) == "low"
    assert plddt_confidence_band(49.9) == "very_low"
    assert plddt_confidence_band(0.0) == "very_low"


def test_plddt_confidence_band_none_is_unavailable() -> None:
    assert plddt_confidence_band(None) == "unavailable"


def test_plddt_bands_constant_has_four_thresholds() -> None:
    assert PLDDT_BANDS == {"very_high": 90.0, "confident": 70.0, "low": 50.0}


def test_alphafold_confidence_label_delegates_and_preserves_vocabulary() -> None:
    # Backward compat: confidence_label still returns the 4 AlphaFold bands.
    assert confidence_label(92.0) == "very_high"
    assert confidence_label(75.0) == "confident"
    assert confidence_label(55.0) == "low"
    assert confidence_label(30.0) == "very_low"


# ---------------------------------------------------------------------------
# Test fixtures
# ---------------------------------------------------------------------------


def _annotation(accession: str = "P04637", sequence: str = "M" * 200) -> UniProtAnnotation:
    return UniProtAnnotation(
        accession=accession,
        protein_name="test protein",
        gene_names=["TEST"],
        organism="Homo sapiens",
        taxonomy_id=9606,
        sequence=sequence,
        functions=["test function"],
        go_terms=[],
        alphafold_url="https://alphafold.ebi.ac.uk/entry/P04637",
        source_url="https://www.uniprot.org/uniprot/P04637",
    )


def _structure(
    mutation_site_plddt: float | None = 92.0,
    mutation: str = "R175H",
    mean_plddt: float = 78.2,
) -> AlphaFoldAnalysis:
    site = None
    if mutation_site_plddt is not None and mutation:
        site = MutationSiteAnalysis(
            wild_type=mutation[0],
            position=int(mutation[1:-1]),
            mutant=mutation[-1],
            plddt=mutation_site_plddt,
            confidence=confidence_label(mutation_site_plddt),
        )
    return AlphaFoldAnalysis(
        accession="P04637",
        entry_id="AF-P04637-F1",
        model_version=4,
        structure_url="https://alphafold.ebi.ac.uk/entry/P04637",
        residue_count=393,
        mean_plddt=mean_plddt,
        confidence_distribution=ConfidenceDistribution(
            very_high=0.4, confident=0.3, low=0.2, very_low=0.1
        ),
        mutation_site=site,
        mutation_neighborhood=[
            MutationNeighbor(position=174, amino_acid="R", distance_angstrom=1.5, plddt=90.0),
            MutationNeighbor(position=176, amino_acid="H", distance_angstrom=2.0, plddt=85.0),
        ],
        contact_map=ContactMap(
            threshold_angstrom=8.0, total_contacts=12, truncated=False, contacts=[]
        ),
        coordinates=[],
    )


def _domains(overlap: bool = True, mutation_position: int = 175) -> DomainAnnotation:
    entry = DomainEntry(
        accession="PF00870",
        name="P53 DNA-binding domain",
        source_database="pfam",
        entry_type="domain",
        integrated_accession="IPR011615",
        source_url="https://www.ebi.ac.uk/interpro/entry/pfam/PF00870/",
        locations=[DomainLocation(start=100, end=300, status="CONTINUOUS")],
        go_terms=[],
        overlaps_mutation=overlap,
    )
    return DomainAnnotation(
        protein_accession="P04637",
        protein_length=393,
        mutation_position=mutation_position,
        entries=[entry],
        mutation_overlaps=[entry] if overlap else [],
        entry_count=1,
        location_count=1,
        source_urls=["https://www.ebi.ac.uk/interpro/api/"],
    )


# ---------------------------------------------------------------------------
# Scenario 5.1: high pLDDT + domain overlap
# ---------------------------------------------------------------------------


def test_high_plddt_with_domain_overlap() -> None:
    summary = build_structure_feature_summary(
        annotation=_annotation(),
        structure=_structure(mutation_site_plddt=92.0, mutation="R175H"),
        domains=_domains(overlap=True, mutation_position=175),
        mutation="R175H",
    )
    assert summary.has_structure is True
    assert summary.structure_source_type == "alphafold_predicted"
    assert summary.mean_plddt == 78.2
    assert summary.mutation_site_plddt == 92.0
    assert summary.mutation_site_confidence_band == "very_high"
    assert summary.low_confidence_region is False
    assert summary.contact_count_8a == 12
    assert summary.neighbor_count_8a == 2
    assert summary.domain_overlap is True
    assert "PF00870" in summary.overlapping_domains[0]
    # Conservative wording: structural context, not functional disruption.
    # Check forbidden impact phrases (not bare substrings — "binding" appears
    # legitimately in the domain name "P53 DNA-binding domain").
    facts_joined = " ".join(summary.evidence_facts).lower()
    for forbidden in ("disrupts", "disrupt", "changes stability", "breaks binding",
                      "causes phenotype", "affects function"):
        assert forbidden not in facts_joined, f"forbidden phrase {forbidden!r} in facts"


# ---------------------------------------------------------------------------
# Scenario 5.2: low pLDDT mutation site
# ---------------------------------------------------------------------------


def test_low_plddt_mutation_site() -> None:
    summary = build_structure_feature_summary(
        annotation=_annotation(),
        structure=_structure(mutation_site_plddt=45.0, mutation="R175H"),
        domains=_domains(overlap=True, mutation_position=175),
        mutation="R175H",
    )
    assert summary.mutation_site_confidence_band == "very_low"
    assert summary.low_confidence_region is True
    limits_joined = " ".join(summary.interpretation_limits).lower()
    assert "low" in limits_joined or "limited" in limits_joined
    # No strong structural-impact wording.
    facts_joined = " ".join(summary.evidence_facts).lower()
    assert "disrupt" not in facts_joined
    assert "causes" not in facts_joined


# ---------------------------------------------------------------------------
# Scenario 5.3: no mutation
# ---------------------------------------------------------------------------


def test_no_mutation_describes_whole_structure_only() -> None:
    summary = build_structure_feature_summary(
        annotation=_annotation(),
        structure=_structure(mutation_site_plddt=None, mutation=None),
        domains=_domains(overlap=False, mutation_position=None),
        mutation=None,
    )
    assert summary.mutation is None
    assert summary.mutation_position is None
    assert summary.mutation_site_plddt is None
    assert summary.mutation_site_confidence_band == "unavailable"
    assert summary.domain_overlap is None or summary.domain_overlap is False
    # No mutation-specific statement in facts.
    facts_joined = " ".join(summary.evidence_facts).lower()
    assert "mutation site" not in facts_joined


# ---------------------------------------------------------------------------
# Scenario 5.4: wild-type mismatch — summary is structure-only, no impact
# ---------------------------------------------------------------------------


def test_wild_type_mismatch_summary_reports_structure_without_impact() -> None:
    # Sequence position 175 is 'A' (from "M"*200, all M), mutation says R175H.
    annotation = _annotation(sequence="M" * 200)
    summary = build_structure_feature_summary(
        annotation=annotation,
        structure=_structure(mutation_site_plddt=92.0, mutation="R175H"),
        domains=_domains(overlap=True, mutation_position=175),
        mutation="R175H",
    )
    # The summary reports structural facts; it does not emit impact statements.
    assert summary.has_structure is True
    assert summary.domain_overlap is True
    facts_joined = " ".join(summary.evidence_facts).lower()
    assert "disrupt" not in facts_joined
    assert "affects function" not in facts_joined
    assert "causes" not in facts_joined


# ---------------------------------------------------------------------------
# Scenario 5.5: no AlphaFold structure
# ---------------------------------------------------------------------------


def test_no_alphafold_structure() -> None:
    annotation = _annotation()
    annotation = annotation.model_copy(update={"alphafold_url": None})
    summary = build_structure_feature_summary(
        annotation=annotation,
        structure=None,
        domains=_domains(overlap=False, mutation_position=175),
        mutation="R175H",
    )
    assert summary.has_structure is False
    assert summary.structure_source_type == "none"
    assert summary.mean_plddt is None
    assert summary.mutation_site_plddt is None
    assert summary.mutation_site_confidence_band == "unavailable"
    assert summary.contact_count_8a is None
    assert summary.neighbor_count_8a is None
    assert any("no" in limit.lower() and "structure" in limit.lower()
               for limit in summary.interpretation_limits)


# ---------------------------------------------------------------------------
# Interpretation limits always present
# ---------------------------------------------------------------------------


def test_interpretation_limits_contain_canonical_disclaimers() -> None:
    summary = build_structure_feature_summary(
        annotation=_annotation(),
        structure=_structure(mutation_site_plddt=92.0, mutation="R175H"),
        domains=_domains(overlap=True, mutation_position=175),
        mutation="R175H",
    )
    limits_joined = " ".join(summary.interpretation_limits).lower()
    assert "predicted structure" in limits_joined
    assert "plddt" in limits_joined
    assert "pathogenicity" in limits_joined or "functional-effect" in limits_joined
    assert "geometric proximity" in limits_joined or "biochemical" in limits_joined


def test_nearest_neighbors_capped_and_formatted() -> None:
    summary = build_structure_feature_summary(
        annotation=_annotation(),
        structure=_structure(mutation_site_plddt=92.0, mutation="R175H"),
        domains=_domains(overlap=True, mutation_position=175),
        mutation="R175H",
    )
    assert len(summary.nearest_neighbors) <= 6
    # Format: single-letter AA + position + distance.
    assert summary.nearest_neighbors[0].startswith("R174") or "174" in summary.nearest_neighbors[0]

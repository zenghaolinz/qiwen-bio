import pytest

from qiwen_bio.mutation import (
    AMINO_ACIDS,
    MutationParseResult,
    parse_mutation,
    wild_type_matches_sequence,
)


def test_parse_standard_mutation() -> None:
    result = parse_mutation("R175H")
    assert result.status == "parsed"
    assert result.normalized == "R175H"
    assert result.wild_type == "R"
    assert result.position == 175
    assert result.mutant == "H"
    assert result.warnings == []


def test_parse_strips_whitespace_and_uppercases() -> None:
    result = parse_mutation("  r175h  ")
    assert result.status == "parsed"
    assert result.normalized == "R175H"
    assert result.wild_type == "R"


def test_parse_accepts_hgvs_short_prefix() -> None:
    result = parse_mutation("p.R175H")
    assert result.status == "parsed"
    assert result.normalized == "R175H"
    assert result.wild_type == "R"
    assert result.position == 175
    assert result.mutant == "H"


def test_parse_rejects_non_amino_acid_residue() -> None:
    # Z is not one of the 20 canonical amino acids.
    result = parse_mutation("Z175H")
    assert result.status == "invalid"
    assert result.normalized is None
    assert result.wild_type is None
    assert result.position is None
    assert result.mutant is None
    assert any("amino acid" in w.lower() or "format" in w.lower() for w in result.warnings)


def test_parse_rejects_missing_position() -> None:
    result = parse_mutation("RHH")
    assert result.status == "invalid"
    assert result.normalized is None
    assert result.warnings


def test_parse_rejects_empty_and_none() -> None:
    assert parse_mutation(None).status == "empty"
    assert parse_mutation("").status == "empty"
    assert parse_mutation("   ").status == "empty"


def test_parse_rejects_gene_prefixed_mutation_without_explicit_request() -> None:
    # "TP53 R175H" is an input-layer concern, not a single mutation token.
    result = parse_mutation("TP53 R175H")
    assert result.status == "invalid"
    assert result.normalized is None


def test_parse_preserves_raw_input() -> None:
    result = parse_mutation("p.r175h")
    assert result.raw == "p.r175h"
    assert result.normalized == "R175H"


def test_wild_type_matches_sequence_at_valid_position() -> None:
    assert wild_type_matches_sequence("MREPK", "R", 2) is True
    assert wild_type_matches_sequence("MREPK", "E", 3) is True


def test_wild_type_matches_sequence_detects_mismatch() -> None:
    # position 2 in "MEEPK" is E, not R
    assert wild_type_matches_sequence("MEEPK", "R", 2) is False


def test_wild_type_matches_sequence_out_of_range() -> None:
    assert wild_type_matches_sequence("ME", "R", 5) is False
    assert wild_type_matches_sequence("ME", "R", 0) is False


def test_amino_acids_constant_is_complete() -> None:
    assert AMINO_ACIDS == frozenset("ACDEFGHIKLMNPQRSTVWY")
    assert len(AMINO_ACIDS) == 20


def test_parse_result_is_pydantic_model() -> None:
    result = parse_mutation("R175H")
    assert isinstance(result, MutationParseResult)
    # Serialization for API responses
    dumped = result.model_dump()
    assert dumped["status"] == "parsed"
    assert dumped["wild_type"] == "R"


def test_invalid_result_has_consistent_null_fields() -> None:
    result = parse_mutation("garbage")
    assert result.status == "invalid"
    assert result.normalized is None
    assert result.wild_type is None
    assert result.position is None
    assert result.mutant is None
    assert len(result.warnings) > 0

"""Unified mutation parser for the Qiwen Bio reasoning pipeline.

A single canonical parser for protein mutation tokens such as ``R175H``.
All modules that need to interpret a mutation string (reasoning chain,
synthesis, API, AlphaFold) must go through :func:`parse_mutation` so the
format validation is consistent and the 20-canonical-amino-acid alphabet is
enforced everywhere.

Biological-correctness note
---------------------------
This parser only validates the *format* of a mutation token. It does **not**
validate that the wild-type residue matches a real protein sequence; that is
the caller's responsibility via :func:`wild_type_matches_sequence`. A
format-valid mutation whose wild type does not match the UniProt sequence
must suppress downstream impact hypotheses (see ADR-0014).
"""

from __future__ import annotations

import re
from typing import Literal

from pydantic import BaseModel


AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

# Strict 20-canonical-amino-acid alphabet. Deliberately rejects letters such
# as Z, B, X, J, U, O that are not one of the 20 standard residues, so a
# malformed or non-biological token never becomes a "parsed" mutation.
_AA = r"[ACDEFGHIKLMNPQRSTVWY]"
# Accept the bare form (R175H) and the HGVS short protein form (p.R175H).
# The optional "p." prefix is stripped on normalization.
_MUTATION_PATTERN = re.compile(rf"^(?:p\.)?({_AA})(\d+)({_AA})$", re.IGNORECASE)


MutationStatus = Literal["parsed", "invalid", "empty"]


class MutationParseResult(BaseModel):
    """Structured result of parsing a mutation token.

    ``status`` is one of:
    - ``parsed``: the token is a well-formed mutation; ``wild_type``,
      ``position``, ``mutant`` and ``normalized`` are populated.
    - ``invalid``: the token was supplied but does not match the expected
      format; all parsed fields are ``None`` and ``warnings`` explains why.
    - ``empty``: no mutation token was supplied (``None`` or whitespace).
    """

    raw: str | None
    normalized: str | None
    wild_type: str | None
    position: int | None
    mutant: str | None
    status: MutationStatus
    warnings: list[str]


def parse_mutation(mutation: str | None) -> MutationParseResult:
    """Parse a mutation token into a structured result.

    Accepts ``R175H`` and ``p.R175H`` (HGVS short protein form, case-
    insensitive, surrounding whitespace tolerated). Rejects non-amino-acid
    letters, missing positions, and multi-token inputs such as ``TP53 R175H``
    (which is an input-layer concern, not a single mutation token).
    """
    if mutation is None or not mutation.strip():
        return MutationParseResult(
            raw=mutation,
            normalized=None,
            wild_type=None,
            position=None,
            mutant=None,
            status="empty",
            warnings=[],
        )

    raw = mutation
    match = _MUTATION_PATTERN.fullmatch(mutation.strip())
    if not match:
        return MutationParseResult(
            raw=raw,
            normalized=None,
            wild_type=None,
            position=None,
            mutant=None,
            status="invalid",
            warnings=[
                "Could not parse mutation format. Expected a single token such "
                "as 'R175H' or 'p.R175H' using the 20 canonical amino acids."
            ],
        )

    wild_type, position_text, mutant = match.groups()
    wild_type = wild_type.upper()
    mutant = mutant.upper()
    position = int(position_text)
    return MutationParseResult(
        raw=raw,
        normalized=f"{wild_type}{position}{mutant}",
        wild_type=wild_type,
        position=position,
        mutant=mutant,
        status="parsed",
        warnings=[],
    )


def wild_type_matches_sequence(
    sequence: str, wild_type: str, position: int
) -> bool:
    """Return ``True`` if ``sequence[position-1] == wild_type``.

    Returns ``False`` for out-of-range positions. Positions are 1-indexed to
    match biological mutation notation (R175H means residue 175).
    """
    if position < 1 or position > len(sequence):
        return False
    return sequence[position - 1] == wild_type

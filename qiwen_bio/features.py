from collections import Counter

from qiwen_bio.models import AMINO_ACIDS, SequenceFeatures


RESIDUE_MASS = {
    "A": 89.09, "C": 121.16, "D": 133.10, "E": 147.13, "F": 165.19,
    "G": 75.07, "H": 155.16, "I": 131.17, "K": 146.19, "L": 131.17,
    "M": 149.21, "N": 132.12, "P": 115.13, "Q": 146.15, "R": 174.20,
    "S": 105.09, "T": 119.12, "V": 117.15, "W": 204.23, "Y": 181.19,
}
WATER_MASS = 18.015
HYDROPHOBIC = frozenset("AVILMFWY")
CHARGED = frozenset("DEKR")
AROMATIC = frozenset("FWY")


def extract_sequence_features(sequence: str) -> SequenceFeatures:
    counts = Counter(sequence)
    length = len(sequence)
    molecular_weight = sum(RESIDUE_MASS[aa] for aa in sequence)
    molecular_weight -= WATER_MASS * max(0, length - 1)

    return SequenceFeatures(
        length=length,
        molecular_weight_da=round(molecular_weight, 2),
        hydrophobic_fraction=round(sum(counts[aa] for aa in HYDROPHOBIC) / length, 4),
        charged_fraction=round(sum(counts[aa] for aa in CHARGED) / length, 4),
        net_charge_proxy=counts["K"] + counts["R"] - counts["D"] - counts["E"],
        aromatic_fraction=round(sum(counts[aa] for aa in AROMATIC) / length, 4),
        composition={aa: round(counts[aa] / length, 4) for aa in sorted(AMINO_ACIDS)},
    )


import re
from statistics import mean

import httpx
from pydantic import BaseModel


ALPHAFOLD_BASE_URL = "https://alphafold.ebi.ac.uk"
MUTATION_PATTERN = re.compile(r"^([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])$")
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


class AlphaFoldNotFoundError(LookupError):
    pass


class AlphaFoldServiceError(RuntimeError):
    pass


class MutationMismatchError(ValueError):
    pass


class ConfidenceDistribution(BaseModel):
    very_high: float
    confident: float
    low: float
    very_low: float


class MutationSiteAnalysis(BaseModel):
    wild_type: str
    position: int
    mutant: str
    plddt: float
    confidence: str


class AlphaFoldAnalysis(BaseModel):
    accession: str
    entry_id: str | None = None
    model_version: int | None = None
    structure_url: str
    residue_count: int
    mean_plddt: float
    confidence_distribution: ConfidenceDistribution
    mutation_site: MutationSiteAnalysis | None = None


def confidence_label(plddt: float) -> str:
    if plddt >= 90:
        return "very_high"
    if plddt >= 70:
        return "confident"
    if plddt >= 50:
        return "low"
    return "very_low"


def parse_alphafold_pdb(
    accession: str,
    pdb_text: str,
    structure_url: str,
    mutation: str | None = None,
) -> AlphaFoldAnalysis:
    residues: dict[int, tuple[str, float]] = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        try:
            residue_name = line[17:20].strip()
            position = int(line[22:26])
            plddt = float(line[60:66])
        except (ValueError, IndexError) as exc:
            raise AlphaFoldServiceError("AlphaFold PDB contains an invalid ATOM record") from exc
        if residue_name in THREE_TO_ONE:
            residues[position] = (THREE_TO_ONE[residue_name], plddt)

    if not residues:
        raise AlphaFoldServiceError("AlphaFold PDB contains no standard CA atoms")

    labels = [confidence_label(plddt) for _, plddt in residues.values()]
    total = len(labels)
    distribution = ConfidenceDistribution(
        very_high=round(labels.count("very_high") / total, 4),
        confident=round(labels.count("confident") / total, 4),
        low=round(labels.count("low") / total, 4),
        very_low=round(labels.count("very_low") / total, 4),
    )
    mutation_site = None
    if mutation:
        match = MUTATION_PATTERN.fullmatch(mutation.strip().upper())
        if not match:
            raise ValueError("mutation must use a format such as R175H")
        wild_type, position_text, mutant = match.groups()
        position = int(position_text)
        if position not in residues:
            raise MutationMismatchError(f"position {position} is absent from the AlphaFold model")
        observed, plddt = residues[position]
        if observed != wild_type:
            raise MutationMismatchError(
                f"expected {wild_type} at position {position}, found {observed}"
            )
        mutation_site = MutationSiteAnalysis(
            wild_type=wild_type,
            position=position,
            mutant=mutant,
            plddt=plddt,
            confidence=confidence_label(plddt),
        )

    return AlphaFoldAnalysis(
        accession=accession,
        structure_url=structure_url,
        residue_count=total,
        mean_plddt=round(mean(plddt for _, plddt in residues.values()), 2),
        confidence_distribution=distribution,
        mutation_site=mutation_site,
    )


class AlphaFoldClient:
    def __init__(
        self,
        base_url: str = ALPHAFOLD_BASE_URL,
        transport: httpx.BaseTransport | None = None,
        timeout: float = 30.0,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.transport = transport
        self.timeout = timeout

    def analyze(self, accession: str, mutation: str | None = None) -> AlphaFoldAnalysis:
        headers = {"Accept": "application/json", "User-Agent": "QiwenBio/0.3"}
        try:
            with httpx.Client(
                base_url=self.base_url,
                headers=headers,
                timeout=self.timeout,
                transport=self.transport,
                follow_redirects=True,
            ) as client:
                metadata_response = client.get(f"/api/prediction/{accession.upper()}")
                if metadata_response.status_code == 404:
                    raise AlphaFoldNotFoundError(f"No AlphaFold model found for {accession}")
                metadata_response.raise_for_status()
                records = metadata_response.json()
                if not records:
                    raise AlphaFoldNotFoundError(f"No AlphaFold model found for {accession}")
                metadata = records[0]
                structure_url = metadata["pdbUrl"]
                structure_response = client.get(
                    structure_url,
                    headers={"Accept": "chemical/x-pdb", "User-Agent": "QiwenBio/0.3"},
                )
                structure_response.raise_for_status()
        except (AlphaFoldNotFoundError, MutationMismatchError):
            raise
        except (httpx.HTTPError, KeyError, ValueError) as exc:
            raise AlphaFoldServiceError(f"AlphaFold request failed: {exc}") from exc

        analysis = parse_alphafold_pdb(
            accession=accession.upper(),
            pdb_text=structure_response.text,
            structure_url=structure_url,
            mutation=mutation,
        )
        analysis.entry_id = metadata.get("entryId")
        analysis.model_version = metadata.get("latestVersion")
        return analysis


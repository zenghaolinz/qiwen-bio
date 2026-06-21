import re
from typing import Any

import httpx
from pydantic import BaseModel


UNIPROT_BASE_URL = "https://rest.uniprot.org"
ACCESSION_PATTERN = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)
IDENTIFIER_PATTERN = re.compile(r"^[A-Za-z0-9_.-]{1,40}$")


class ProteinNotFoundError(LookupError):
    pass


class AmbiguousProteinError(LookupError):
    pass


class UniProtServiceError(RuntimeError):
    pass


class GOTerm(BaseModel):
    id: str
    name: str
    aspect: str


class UniProtAnnotation(BaseModel):
    accession: str
    protein_name: str
    gene_names: list[str]
    organism: str
    taxonomy_id: int
    sequence: str
    functions: list[str]
    go_terms: list[GOTerm]
    alphafold_url: str | None
    source_url: str


def _protein_name(record: dict[str, Any]) -> str:
    description = record.get("proteinDescription", {})
    recommended = description.get("recommendedName", {}).get("fullName", {}).get("value")
    if recommended:
        return recommended
    submissions = description.get("submissionNames") or []
    if submissions:
        return submissions[0].get("fullName", {}).get("value", "Uncharacterized protein")
    return "Uncharacterized protein"


def parse_uniprot_record(record: dict[str, Any]) -> UniProtAnnotation:
    accession = record["primaryAccession"]
    gene_names: list[str] = []
    for gene in record.get("genes", []):
        values = [gene.get("geneName", {}).get("value")]
        values.extend(item.get("value") for item in gene.get("synonyms", []))
        gene_names.extend(value for value in values if value and value not in gene_names)

    functions = [
        text["value"]
        for comment in record.get("comments", [])
        if comment.get("commentType") == "FUNCTION"
        for text in comment.get("texts", [])
        if text.get("value")
    ]
    aspect_names = {
        "F": "molecular_function",
        "P": "biological_process",
        "C": "cellular_component",
    }
    go_terms: list[GOTerm] = []
    has_alphafold = False
    for reference in record.get("uniProtKBCrossReferences", []):
        if reference.get("database") == "AlphaFoldDB":
            has_alphafold = True
        if reference.get("database") != "GO":
            continue
        term_value = next(
            (
                item.get("value", "")
                for item in reference.get("properties", [])
                if item.get("key") == "GoTerm"
            ),
            "",
        )
        prefix, _, name = term_value.partition(":")
        go_terms.append(
            GOTerm(
                id=reference["id"],
                name=name or term_value or "Unlabelled GO term",
                aspect=aspect_names.get(prefix, "unknown"),
            )
        )

    organism = record.get("organism", {})
    return UniProtAnnotation(
        accession=accession,
        protein_name=_protein_name(record),
        gene_names=gene_names,
        organism=organism.get("scientificName", "Unknown organism"),
        taxonomy_id=organism.get("taxonId", 0),
        sequence=record["sequence"]["value"],
        functions=functions,
        go_terms=go_terms,
        alphafold_url=(
            f"https://alphafold.ebi.ac.uk/entry/{accession}" if has_alphafold else None
        ),
        source_url=f"{UNIPROT_BASE_URL}/uniprotkb/{accession}",
    )


class UniProtClient:
    def __init__(
        self,
        base_url: str = UNIPROT_BASE_URL,
        transport: httpx.BaseTransport | None = None,
        timeout: float = 15.0,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.transport = transport
        self.timeout = timeout

    def resolve(self, identifier: str, organism_id: int = 9606) -> UniProtAnnotation:
        normalized = identifier.strip().upper()
        if not IDENTIFIER_PATTERN.fullmatch(normalized):
            raise ValueError("identifier contains unsupported characters")

        if ACCESSION_PATTERN.fullmatch(normalized):
            payload = self._get(f"/uniprotkb/{normalized}.json")
            return parse_uniprot_record(payload)

        query = f"(gene_exact:{normalized}) AND (organism_id:{organism_id}) AND (reviewed:true)"
        payload = self._get(
            "/uniprotkb/search",
            params={"query": query, "format": "json", "size": 2},
        )
        results = payload.get("results", [])
        if not results:
            raise ProteinNotFoundError(
                f"No reviewed UniProt entry found for {normalized} in taxonomy {organism_id}"
            )
        if len(results) > 1:
            raise AmbiguousProteinError(
                f"Multiple reviewed UniProt entries found for {normalized} in taxonomy {organism_id}"
            )
        return parse_uniprot_record(results[0])

    def _get(self, path: str, params: dict[str, Any] | None = None) -> dict[str, Any]:
        headers = {"Accept": "application/json", "User-Agent": "QiwenBio/0.2"}
        try:
            with httpx.Client(
                base_url=self.base_url,
                headers=headers,
                timeout=self.timeout,
                transport=self.transport,
            ) as client:
                response = client.get(path, params=params)
        except httpx.HTTPError as exc:
            raise UniProtServiceError(f"UniProt request failed: {exc}") from exc
        if response.status_code == 404:
            raise ProteinNotFoundError("UniProt entry was not found")
        try:
            response.raise_for_status()
            return response.json()
        except (httpx.HTTPError, ValueError) as exc:
            raise UniProtServiceError(
                f"UniProt returned an invalid response ({response.status_code})"
            ) from exc


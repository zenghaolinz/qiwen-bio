from __future__ import annotations

import re

import httpx
from pydantic import BaseModel


class KeggServiceError(RuntimeError):
    pass


class KeggNotFoundError(LookupError):
    pass


class KeggPathway(BaseModel):
    pathway_id: str
    name: str
    description: str
    classes: list[str]
    source_url: str


class KeggPathwayAnnotation(BaseModel):
    protein_accession: str
    gene_ids: list[str]
    pathways: list[KeggPathway]
    pathway_count: int
    linked_pathway_count: int
    truncated: bool
    query_urls: list[str]
    copyright_url: str = "https://www.kegg.jp/kegg/legal.html"
    disclaimer: str = (
        "KEGG pathway membership is a direct database association, not evidence that the "
        "protein activates the pathway or causes a cellular phenotype."
    )


def _parse_flat_file(text: str) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    for raw_record in text.split("///"):
        fields: dict[str, str] = {}
        current_key: str | None = None
        for line in raw_record.splitlines():
            if not line.strip():
                continue
            key = line[:12].strip()
            value = line[12:].strip()
            if key:
                current_key = key
                fields[key] = (
                    f"{fields[key]} {value}" if key in fields and value else value
                )
            elif current_key and value:
                fields[current_key] = f"{fields[current_key]} {value}".strip()
        if fields:
            records.append(fields)
    return records


class KeggClient:
    base_url = "https://rest.kegg.jp"

    def __init__(self, http_client: httpx.Client | None = None) -> None:
        self.http_client = http_client or httpx.Client(
            timeout=30.0,
            follow_redirects=False,
            headers={"User-Agent": "QiwenBio/0.1 (academic research prototype)"},
        )

    def _get(self, path: str) -> tuple[str, str]:
        url = f"{self.base_url}{path}"
        try:
            response = self.http_client.get(url)
            if response.status_code == 404:
                return "", url
            response.raise_for_status()
        except httpx.HTTPError as exc:
            raise KeggServiceError(f"KEGG REST request failed: {exc}") from exc
        return response.text, url

    def fetch(self, accession: str, limit: int = 20) -> KeggPathwayAnnotation:
        normalized_accession = accession.upper()
        if not re.fullmatch(r"[A-Z0-9]{6,10}", normalized_accession):
            raise ValueError("invalid UniProt accession")
        if not 1 <= limit <= 50:
            raise ValueError("KEGG pathway limit must be between 1 and 50")

        query_urls: list[str] = []
        conversion_text, conversion_url = self._get(
            f"/conv/genes/uniprot:{normalized_accession}"
        )
        query_urls.append(conversion_url)
        gene_ids: list[str] = []
        for line in conversion_text.splitlines():
            columns = line.split("\t")
            if len(columns) == 2 and columns[1] and ":" in columns[1]:
                gene_ids.append(columns[1])
        gene_ids = list(dict.fromkeys(gene_ids))
        if not gene_ids:
            raise KeggNotFoundError(
                f"No KEGG gene mapping found for UniProt {normalized_accession}"
            )

        pathway_ids: list[str] = []
        for gene_id in gene_ids[:20]:
            link_text, link_url = self._get(f"/link/pathway/{gene_id}")
            query_urls.append(link_url)
            for line in link_text.splitlines():
                columns = line.split("\t")
                if len(columns) != 2 or not columns[1].startswith("path:"):
                    continue
                pathway_ids.append(columns[1].removeprefix("path:"))
        all_pathway_ids = list(dict.fromkeys(pathway_ids))
        if not all_pathway_ids:
            raise KeggNotFoundError(
                f"No KEGG pathways found for UniProt {normalized_accession}"
            )
        pathway_ids = all_pathway_ids[:limit]

        records: dict[str, KeggPathway] = {}
        for offset in range(0, len(pathway_ids), 10):
            batch = pathway_ids[offset : offset + 10]
            flat_text, get_url = self._get(f"/get/{'+'.join(batch)}")
            query_urls.append(get_url)
            for fields in _parse_flat_file(flat_text):
                entry = fields.get("ENTRY", "").split()
                if not entry:
                    continue
                pathway_id = entry[0]
                if pathway_id not in pathway_ids:
                    continue
                records[pathway_id] = KeggPathway(
                    pathway_id=pathway_id,
                    name=fields.get("NAME", pathway_id),
                    description=fields.get("DESCRIPTION", ""),
                    classes=[
                        item.strip()
                        for item in fields.get("CLASS", "").split(";")
                        if item.strip()
                    ],
                    source_url=f"https://www.kegg.jp/entry/{pathway_id}",
                )
        missing = [pathway_id for pathway_id in pathway_ids if pathway_id not in records]
        if missing:
            raise KeggServiceError(
                f"KEGG REST omitted metadata for {len(missing)} linked pathways"
            )
        pathways = [records[pathway_id] for pathway_id in pathway_ids]
        return KeggPathwayAnnotation(
            protein_accession=normalized_accession,
            gene_ids=gene_ids,
            pathways=pathways,
            pathway_count=len(pathways),
            linked_pathway_count=len(all_pathway_ids),
            truncated=len(all_pathway_ids) > len(pathways),
            query_urls=query_urls,
        )

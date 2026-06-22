from __future__ import annotations

from typing import Any, Literal
from urllib.parse import urlsplit

import httpx
from pydantic import BaseModel


class InterProServiceError(RuntimeError):
    pass


class InterProNotFoundError(LookupError):
    pass


class DomainGoTerm(BaseModel):
    identifier: str
    name: str
    category: str


class DomainLocation(BaseModel):
    start: int
    end: int
    status: str
    model: str | None = None
    score: float | None = None

    def contains(self, position: int) -> bool:
        return self.start <= position <= self.end


class DomainEntry(BaseModel):
    accession: str
    name: str
    source_database: Literal["interpro", "pfam"]
    entry_type: str
    integrated_accession: str | None
    source_url: str
    locations: list[DomainLocation]
    go_terms: list[DomainGoTerm]
    overlaps_mutation: bool = False


class DomainAnnotation(BaseModel):
    protein_accession: str
    protein_length: int
    mutation_position: int | None
    entries: list[DomainEntry]
    mutation_overlaps: list[DomainEntry]
    entry_count: int
    location_count: int
    source_urls: list[str]
    disclaimer: str = (
        "Domain coordinates are database annotations. Positional overlap does not by itself "
        "establish functional impact, pathogenicity, or altered binding."
    )


class InterProClient:
    endpoint_template = (
        "https://www.ebi.ac.uk/interpro/api/entry/{source}/protein/uniprot/{accession}/"
    )
    entry_url_template = "https://www.ebi.ac.uk/interpro/entry/{source}/{entry}/"
    sources: tuple[Literal["interpro", "pfam"], ...] = ("interpro", "pfam")

    def __init__(self, http_client: httpx.Client | None = None) -> None:
        self.http_client = http_client or httpx.Client(
            timeout=30.0,
            follow_redirects=True,
            headers={"User-Agent": "QiwenBio/0.1 (evidence-oriented research prototype)"},
        )

    def _pages(self, initial_url: str) -> list[dict[str, Any]]:
        url: str | None = initial_url
        first = True
        visited: set[str] = set()
        results: list[dict[str, Any]] = []
        for _ in range(100):
            if url is None:
                return results
            parsed_url = urlsplit(url)
            if (
                parsed_url.scheme != "https"
                or parsed_url.hostname != "www.ebi.ac.uk"
                or parsed_url.port not in (None, 443)
                or not parsed_url.path.startswith("/interpro/api/")
            ):
                raise InterProServiceError("InterPro API returned an unsafe pagination URL")
            if url in visited:
                raise InterProServiceError("InterPro API pagination loop detected")
            visited.add(url)
            try:
                response = self.http_client.get(
                    url, params={"page_size": 200} if first else None
                )
                first = False
                if response.status_code == 404:
                    return results
                response.raise_for_status()
                payload = response.json()
            except (httpx.HTTPError, ValueError, TypeError) as exc:
                raise InterProServiceError(f"InterPro API request failed: {exc}") from exc
            page_results = payload.get("results")
            if not isinstance(page_results, list):
                raise InterProServiceError("InterPro API returned an invalid results payload")
            results.extend(item for item in page_results if isinstance(item, dict))
            next_url = payload.get("next")
            if next_url is not None and not isinstance(next_url, str):
                raise InterProServiceError("InterPro API returned an invalid next URL")
            url = next_url
        raise InterProServiceError("InterPro API exceeded the 100-page safety limit")

    def fetch(
        self, accession: str, mutation_position: int | None = None
    ) -> DomainAnnotation:
        normalized_accession = accession.upper()
        merged: dict[tuple[str, str], DomainEntry] = {}
        protein_lengths: list[int] = []
        source_urls: list[str] = []
        for source in self.sources:
            source_url = self.endpoint_template.format(
                source=source, accession=normalized_accession
            )
            source_urls.append(source_url)
            for result in self._pages(source_url):
                metadata = result.get("metadata") or {}
                entry_accession = metadata.get("accession")
                if not isinstance(entry_accession, str) or not entry_accession:
                    continue
                locations: list[DomainLocation] = []
                for protein in result.get("proteins") or []:
                    if str(protein.get("accession", "")).upper() != normalized_accession:
                        continue
                    protein_length = protein.get("protein_length")
                    if isinstance(protein_length, int) and protein_length > 0:
                        protein_lengths.append(protein_length)
                    for protein_location in protein.get("entry_protein_locations") or []:
                        model = protein_location.get("model")
                        score = protein_location.get("score")
                        for fragment in protein_location.get("fragments") or []:
                            start = fragment.get("start")
                            end = fragment.get("end")
                            if (
                                isinstance(start, int)
                                and isinstance(end, int)
                                and start >= 1
                                and end >= start
                            ):
                                locations.append(
                                    DomainLocation(
                                        start=start,
                                        end=end,
                                        status=str(fragment.get("dc-status") or "UNKNOWN"),
                                        model=str(model) if model is not None else None,
                                        score=float(score) if score is not None else None,
                                    )
                                )
                if not locations:
                    continue
                metadata_source = str(metadata.get("source_database") or source).lower()
                if metadata_source not in self.sources:
                    metadata_source = source
                key = (metadata_source, entry_accession)
                go_terms = [
                    DomainGoTerm(
                        identifier=str(term.get("identifier") or ""),
                        name=str(term.get("name") or ""),
                        category=str((term.get("category") or {}).get("name") or ""),
                    )
                    for term in metadata.get("go_terms") or []
                    if term.get("identifier")
                ]
                entry = merged.get(key)
                if entry is None:
                    entry = DomainEntry(
                        accession=entry_accession,
                        name=str(metadata.get("name") or entry_accession),
                        source_database=metadata_source,
                        entry_type=str(metadata.get("type") or "unknown"),
                        integrated_accession=(
                            str(metadata["integrated"])
                            if metadata.get("integrated") is not None
                            else None
                        ),
                        source_url=self.entry_url_template.format(
                            source=metadata_source, entry=entry_accession
                        ),
                        locations=[],
                        go_terms=go_terms,
                    )
                    merged[key] = entry
                existing = {
                    (item.start, item.end, item.status, item.model, item.score)
                    for item in entry.locations
                }
                for location in locations:
                    identity = (
                        location.start,
                        location.end,
                        location.status,
                        location.model,
                        location.score,
                    )
                    if identity not in existing:
                        entry.locations.append(location)
                        existing.add(identity)

        if not merged:
            raise InterProNotFoundError(
                f"No InterPro or Pfam entries found for {normalized_accession}"
            )
        entries = sorted(
            merged.values(), key=lambda item: (item.source_database, item.accession)
        )
        for entry in entries:
            entry.locations.sort(key=lambda item: (item.start, item.end))
            entry.overlaps_mutation = bool(
                mutation_position
                and any(location.contains(mutation_position) for location in entry.locations)
            )
        overlaps = [entry for entry in entries if entry.overlaps_mutation]
        return DomainAnnotation(
            protein_accession=normalized_accession,
            protein_length=max(protein_lengths, default=0),
            mutation_position=mutation_position,
            entries=entries,
            mutation_overlaps=overlaps,
            entry_count=len(entries),
            location_count=sum(len(entry.locations) for entry in entries),
            source_urls=source_urls,
        )

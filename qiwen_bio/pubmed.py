import re
from typing import Any

import httpx
from pydantic import BaseModel


EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


class PubMedServiceError(RuntimeError):
    pass


class PubMedArticle(BaseModel):
    pmid: str
    title: str
    authors: list[str]
    journal: str
    published: str
    doi: str | None
    url: str


class LiteratureEvidence(BaseModel):
    protein: str
    context_terms: list[str]
    query: str
    articles: list[PubMedArticle]
    disclaimer: str = (
        "Search retrieval does not by itself validate a biological claim; read and appraise "
        "each cited article before using it as evidence."
    )


def _clean_term(value: str, max_length: int) -> str:
    cleaned = re.sub(r'["\[\]\r\n]+', " ", value)
    return " ".join(cleaned.split())[:max_length]


def build_pubmed_query(protein: str, context_terms: list[str]) -> tuple[str, list[str]]:
    clean_protein = _clean_term(protein, 40).upper()
    clean_context = [
        term for term in (_clean_term(item, 120) for item in context_terms[:5]) if term
    ]
    protein_clause = f'"{clean_protein}"[Title/Abstract]'
    if not clean_context:
        return protein_clause, clean_context
    context_clause = " OR ".join(f'"{term}"[Title/Abstract]' for term in clean_context)
    return f"({protein_clause}) AND ({context_clause})", clean_context


def _parse_article(record: dict[str, Any]) -> PubMedArticle:
    pmid = str(record["uid"])
    doi = next(
        (
            item.get("value")
            for item in record.get("articleids", [])
            if item.get("idtype") == "doi"
        ),
        None,
    )
    return PubMedArticle(
        pmid=pmid,
        title=record.get("title") or "Untitled PubMed record",
        authors=[item["name"] for item in record.get("authors", [])[:10] if item.get("name")],
        journal=record.get("fulljournalname") or record.get("source") or "Unknown journal",
        published=record.get("pubdate") or "Unknown date",
        doi=doi,
        url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
    )


class PubMedClient:
    def __init__(
        self,
        base_url: str = EUTILS_BASE_URL,
        transport: httpx.BaseTransport | None = None,
        timeout: float = 20.0,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.transport = transport
        self.timeout = timeout

    def search(
        self,
        protein: str,
        context_terms: list[str] | None = None,
        limit: int = 5,
    ) -> LiteratureEvidence:
        query, clean_context = build_pubmed_query(protein, context_terms or [])
        search_payload = self._get(
            "/esearch.fcgi",
            {
                "db": "pubmed",
                "term": query,
                "retmode": "json",
                "retmax": limit,
                "sort": "relevance",
                "tool": "qiwen_bio",
            },
        )
        ids = search_payload.get("esearchresult", {}).get("idlist", [])
        if not ids:
            return LiteratureEvidence(
                protein=protein.upper(),
                context_terms=clean_context,
                query=query,
                articles=[],
            )

        summary_payload = self._get(
            "/esummary.fcgi",
            {
                "db": "pubmed",
                "id": ",".join(ids),
                "retmode": "json",
                "tool": "qiwen_bio",
            },
        )
        result = summary_payload.get("result", {})
        articles = [
            _parse_article(result[pmid])
            for pmid in result.get("uids", ids)
            if pmid in result
        ]
        return LiteratureEvidence(
            protein=protein.upper(),
            context_terms=clean_context,
            query=query,
            articles=articles,
        )

    def _get(self, path: str, params: dict[str, Any]) -> dict[str, Any]:
        try:
            with httpx.Client(
                base_url=self.base_url,
                headers={"Accept": "application/json", "User-Agent": "QiwenBio/0.5"},
                timeout=self.timeout,
                transport=self.transport,
                follow_redirects=True,
            ) as client:
                response = client.get(path, params=params)
            response.raise_for_status()
            payload = response.json()
            if not isinstance(payload, dict):
                raise ValueError("expected a JSON object")
            return payload
        except (httpx.HTTPError, ValueError) as exc:
            raise PubMedServiceError(f"PubMed request failed: {exc}") from exc


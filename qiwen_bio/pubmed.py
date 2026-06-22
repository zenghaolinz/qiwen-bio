import re
from typing import Any
from xml.etree import ElementTree

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
    abstract: str | None = None
    abstract_sections: list[tuple[str | None, str]] = []


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


def _parse_efetch_xml(
    text: str,
) -> dict[str, tuple[str | None, list[tuple[str | None, str]]]]:
    """Parse a PubMed efetch XML set into ``{pmid: (abstract, sections)}``.

    The abstract is the concatenation of every ``AbstractText`` child (in
    document order) joined by single spaces. Each section is a
    ``(label, text)`` tuple where ``label`` is the ``Label`` attribute or
    ``None``. Malformed XML or a missing Abstract element yields no entry; the
    caller treats a missing PMID as "no abstract available".
    """
    try:
        root = ElementTree.fromstring(text)
    except ElementTree.ParseError:
        return {}
    result: dict[str, tuple[str | None, list[tuple[str | None, str]]]] = {}
    for article in root.findall(".//PubmedArticle"):
        pmid_el = article.find(".//PMID")
        if pmid_el is None or not (pmid_el.text or "").strip():
            continue
        pmid = pmid_el.text.strip()
        sections: list[tuple[str | None, str]] = []
        for abstract_text in article.findall(".//Abstract/AbstractText"):
            label = abstract_text.get("Label")
            # ElementTree joins nested markup tails poorly, so reassemble all
            # text descendants to recover the full visible text.
            parts = [abstract_text.text or ""]
            for descendant in abstract_text.iter():
                if descendant is abstract_text:
                    continue
                if descendant.text:
                    parts.append(descendant.text)
                if descendant.tail:
                    parts.append(descendant.tail)
            text_value = " ".join(part for part in parts if part).strip()
            text_value = re.sub(r"\s+", " ", text_value)
            if text_value:
                sections.append((label, text_value))
        abstract = " ".join(section_text for _, section_text in sections).strip()
        abstract = re.sub(r"\s+", " ", abstract) or None
        result[pmid] = (abstract, sections)
    return result


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

    def search_with_abstracts(
        self,
        protein: str,
        context_terms: list[str] | None = None,
        limit: int = 5,
    ) -> LiteratureEvidence:
        """Run the same esearch+esummary path as ``search`` and enrich each
        article with its abstract via efetch (retmode=xml).

        Abstracts are parsed at the abstract-metadata boundary only. The
        structured ``AbstractText`` sections (with their ``Label`` attribute)
        are preserved so downstream layers can cite the exact excerpt. Records
        that have no abstract keep ``abstract=None`` and an empty section list.
        """
        evidence = self.search(protein, context_terms=context_terms, limit=limit)
        if not evidence.articles:
            return evidence
        pmids = [article.pmid for article in evidence.articles]
        abstracts = self._fetch_abstracts(pmids)
        enriched = [
            article.model_copy(
                update={
                    "abstract": abstracts[article.pmid][0],
                    "abstract_sections": abstracts[article.pmid][1],
                }
            )
            for article in evidence.articles
        ]
        return LiteratureEvidence(
            protein=evidence.protein,
            context_terms=evidence.context_terms,
            query=evidence.query,
            articles=enriched,
        )

    def _fetch_abstracts(self, pmids: list[str]) -> dict[str, tuple[str | None, list[tuple[str | None, str]]]]:
        if not pmids:
            return {}
        try:
            with httpx.Client(
                base_url=self.base_url,
                headers={"Accept": "application/xml", "User-Agent": "QiwenBio/0.5"},
                timeout=self.timeout,
                transport=self.transport,
                follow_redirects=True,
            ) as client:
                response = client.get(
                    "/efetch.fcgi",
                    params={
                        "db": "pubmed",
                        "id": ",".join(pmids),
                        "retmode": "xml",
                        "tool": "qiwen_bio",
                    },
                )
            response.raise_for_status()
        except httpx.HTTPError as exc:
            raise PubMedServiceError(f"PubMed efetch failed: {exc}") from exc
        return _parse_efetch_xml(response.text)

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


from typing import Any, Literal

import httpx
from pydantic import BaseModel


STRING_BASE_URL = "https://string-db.org"
NETWORK_SOURCE_URL = f"{STRING_BASE_URL}/api/json/network"
ENRICHMENT_SOURCE_URL = f"{STRING_BASE_URL}/api/json/enrichment"
SUPPORTED_TERM_TYPES = {
    "Process": "process",
    "KEGG": "pathway",
    "Reactome Pathways": "pathway",
    "WikiPathways": "pathway",
}


class StringNotFoundError(LookupError):
    pass


class StringServiceError(RuntimeError):
    pass


class GraphNode(BaseModel):
    id: str
    type: Literal["protein", "process", "pathway"]
    label: str
    external_id: str
    organism_id: int
    fdr: float | None = None
    source_url: str


class GraphEdge(BaseModel):
    id: str
    type: Literal["interacts_with", "annotated_to"]
    source: str
    target: str
    score: float | None = None
    fdr: float | None = None
    source_name: str
    source_url: str


class EvidenceGraph(BaseModel):
    seed: str
    organism_id: int
    nodes: list[GraphNode]
    edges: list[GraphEdge]


def build_evidence_graph(
    seed: str,
    organism_id: int,
    network_records: list[dict[str, Any]],
    enrichment_records: list[dict[str, Any]],
) -> EvidenceGraph:
    nodes: dict[str, GraphNode] = {}
    label_to_id: dict[str, str] = {}
    interaction_edges: dict[str, GraphEdge] = {}

    for record in network_records:
        endpoint_ids: list[str] = []
        for suffix in ("A", "B"):
            string_id = record[f"stringId_{suffix}"]
            label = record[f"preferredName_{suffix}"]
            node_id = f"string:{string_id}"
            endpoint_ids.append(node_id)
            nodes[node_id] = GraphNode(
                id=node_id,
                type="protein",
                label=label,
                external_id=string_id,
                organism_id=organism_id,
                source_url=f"{STRING_BASE_URL}/network/{string_id}",
            )
            label_to_id[label.upper()] = node_id

        source_id, target_id = sorted(endpoint_ids)
        edge_id = f"interaction:{source_id}|{target_id}"
        score = float(record.get("score", 0))
        previous = interaction_edges.get(edge_id)
        if previous is None or score > (previous.score or 0):
            interaction_edges[edge_id] = GraphEdge(
                id=edge_id,
                type="interacts_with",
                source=source_id,
                target=target_id,
                score=score,
                source_name="STRING",
                source_url=NETWORK_SOURCE_URL,
            )

    term_edges: dict[str, GraphEdge] = {}
    supported = [
        record for record in enrichment_records if record.get("category") in SUPPORTED_TERM_TYPES
    ]
    supported.sort(key=lambda item: float(item.get("fdr", 1)))
    for record in supported[:10]:
        category = record["category"]
        external_id = record["term"]
        node_id = f"term:{category}:{external_id}"
        fdr = float(record.get("fdr", 1))
        nodes[node_id] = GraphNode(
            id=node_id,
            type=SUPPORTED_TERM_TYPES[category],
            label=record.get("description") or external_id,
            external_id=external_id,
            organism_id=organism_id,
            fdr=fdr,
            source_url=ENRICHMENT_SOURCE_URL,
        )
        for label in record.get("preferredNames", []):
            protein_id = label_to_id.get(label.upper())
            if protein_id is None:
                continue
            edge_id = f"annotation:{protein_id}|{node_id}"
            term_edges[edge_id] = GraphEdge(
                id=edge_id,
                type="annotated_to",
                source=protein_id,
                target=node_id,
                fdr=fdr,
                source_name="STRING enrichment",
                source_url=ENRICHMENT_SOURCE_URL,
            )

    return EvidenceGraph(
        seed=seed.upper(),
        organism_id=organism_id,
        nodes=list(nodes.values()),
        edges=[*interaction_edges.values(), *term_edges.values()],
    )


class StringClient:
    def __init__(
        self,
        base_url: str = STRING_BASE_URL,
        transport: httpx.BaseTransport | None = None,
        timeout: float = 20.0,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.transport = transport
        self.timeout = timeout

    def build_graph(
        self,
        identifier: str,
        species: int = 9606,
        limit: int = 10,
        required_score: int = 700,
    ) -> EvidenceGraph:
        normalized = identifier.strip().upper()
        network = self._get(
            "/api/json/network",
            {
                "identifiers": normalized,
                "species": species,
                "required_score": required_score,
                "limit": limit,
                "caller_identity": "qiwen_bio",
            },
        )
        if not network:
            raise StringNotFoundError(f"STRING found no network for {normalized}")
        protein_names = sorted(
            {
                record[key]
                for record in network
                for key in ("preferredName_A", "preferredName_B")
            }
        )
        enrichment = self._get(
            "/api/json/enrichment",
            {
                "identifiers": "\r".join(protein_names),
                "species": species,
                "caller_identity": "qiwen_bio",
            },
        )
        return build_evidence_graph(normalized, species, network, enrichment)

    def _get(self, path: str, params: dict[str, Any]) -> list[dict[str, Any]]:
        try:
            with httpx.Client(
                base_url=self.base_url,
                headers={"Accept": "application/json", "User-Agent": "QiwenBio/0.4"},
                timeout=self.timeout,
                transport=self.transport,
                follow_redirects=True,
            ) as client:
                response = client.get(path, params=params)
            response.raise_for_status()
            payload = response.json()
            if not isinstance(payload, list):
                raise ValueError("expected a JSON list")
            return payload
        except (httpx.HTTPError, ValueError) as exc:
            raise StringServiceError(f"STRING request failed: {exc}") from exc


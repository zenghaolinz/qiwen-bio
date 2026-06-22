from qiwen_bio.cellular_processes import build_cellular_process_evidence
from qiwen_bio.kegg import KeggPathway, KeggPathwayAnnotation
from qiwen_bio.stringdb import build_evidence_graph
from qiwen_bio.uniprot import GOTerm, UniProtAnnotation
from tests.test_string_graph import NETWORK_RECORDS


def _annotation(go_terms: list[GOTerm]) -> UniProtAnnotation:
    return UniProtAnnotation(
        accession="P04637",
        protein_name="Cellular tumor antigen p53",
        gene_names=["TP53"],
        organism="Homo sapiens",
        taxonomy_id=9606,
        sequence="MEEPQSDPSV",
        functions=[],
        go_terms=go_terms,
        alphafold_url=None,
        source_url="https://rest.uniprot.org/uniprotkb/P04637",
    )


def _kegg() -> KeggPathwayAnnotation:
    return KeggPathwayAnnotation(
        protein_accession="P04637",
        gene_ids=["hsa:7157"],
        pathways=[
            KeggPathway(
                pathway_id="hsa04115",
                name="p53 signaling pathway - Homo sapiens (human)",
                description="p53 stress response.",
                classes=["Cellular Processes", "Cell growth and death"],
                source_url="https://www.kegg.jp/entry/hsa04115",
            )
        ],
        pathway_count=1,
        linked_pathway_count=1,
        truncated=False,
        query_urls=["https://rest.kegg.jp/link/pathway/hsa:7157"],
    )


def test_builder_merges_stable_ids_and_filters_string_to_seed_edges() -> None:
    enrichment = [
        {
            "category": "Process",
            "term": "GO:0072331",
            "description": "Signal transduction by p53 class mediator",
            "fdr": 1.2e-8,
            "preferredNames": ["TP53", "MDM2"],
        },
        {
            "category": "KEGG",
            "term": "hsa04115",
            "description": "p53 signaling pathway",
            "fdr": 2.4e-7,
            "preferredNames": ["TP53", "MDM2"],
        },
        {
            "category": "Process",
            "term": "GO:9999999",
            "description": "Neighbor-only process",
            "fdr": 1e-9,
            "preferredNames": ["MDM2"],
        },
    ]
    graph = build_evidence_graph("TP53", 9606, NETWORK_RECORDS, enrichment)
    annotation = _annotation(
        [
            GOTerm(
                id="GO:0072331",
                name="signal transduction by p53 class mediator",
                aspect="biological_process",
            ),
            GOTerm(id="GO:0003677", name="DNA binding", aspect="molecular_function"),
        ]
    )

    evidence = build_cellular_process_evidence(annotation, _kegg(), graph)

    assert [process.canonical_id for process in evidence.processes] == [
        "GO:0072331",
        "hsa04115",
    ]
    go_process = evidence.processes[0]
    assert [support.source_name for support in go_process.supports] == [
        "UniProt GO",
        "STRING enrichment",
    ]
    assert go_process.supports[1].fdr == 1.2e-8
    kegg_process = evidence.processes[1]
    assert [support.evidence_type for support in kegg_process.supports] == [
        "database_membership",
        "enrichment_statistic",
    ]
    assert "GO:9999999" not in {process.canonical_id for process in evidence.processes}
    assert evidence.counts_by_source == {
        "KEGG": 1,
        "STRING enrichment": 2,
        "UniProt GO": 1,
    }
    assert evidence.phenotype_hypotheses == []
    assert evidence.literature_context == [
        "signal transduction by p53 class mediator",
        "p53 signaling pathway - Homo sapiens (human)",
    ]
    assert "does not establish" in evidence.interpretation_boundary


def test_builder_returns_explicit_empty_layer_without_process_inputs() -> None:
    evidence = build_cellular_process_evidence(
        _annotation(
            [GOTerm(id="GO:0003677", name="DNA binding", aspect="molecular_function")]
        ),
        None,
        None,
    )

    assert evidence.processes == []
    assert evidence.counts_by_source == {}
    assert evidence.phenotype_hypotheses == []

from qiwen_bio.alphafold import parse_alphafold_pdb
from qiwen_bio.cellular_processes import (
    CellularProcess,
    CellularProcessEvidence,
    ProcessSupport,
    build_cellular_process_evidence,
)
from qiwen_bio.interpro import (
    DomainAnnotation,
    DomainEntry,
    DomainLocation,
)
from qiwen_bio.kegg import KeggPathway, KeggPathwayAnnotation
from qiwen_bio.models import AnalysisRequest
from qiwen_bio.phenotype_literature import (
    PhenotypeClaimLink,
    PhenotypeLiteratureEvidence,
)
from qiwen_bio.pipeline import AnalysisPipeline
from qiwen_bio.reasoning_chain import build_reasoning_chain
from qiwen_bio.stringdb import build_evidence_graph
from qiwen_bio.uniprot import parse_uniprot_record
from tests.test_alphafold import PDB_TEXT
from tests.test_string_graph import ENRICHMENT_RECORDS, NETWORK_RECORDS
from tests.test_uniprot import UNIPROT_RECORD


def _annotation():
    return parse_uniprot_record(UNIPROT_RECORD)


def _analysis(mutation=None):
    ann = _annotation()
    return AnalysisPipeline().analyze(
        AnalysisRequest(
            name=ann.protein_name,
            sequence=ann.sequence,
            mutation=mutation,
        )
    )


def _structure(mutation=None):
    return parse_alphafold_pdb(
        "P04637", PDB_TEXT, "https://example.test/model.pdb", mutation
    )


def _domains(mutation_position=None, overlap=True):
    entry = DomainEntry(
        accession="PF00870",
        name="P53 DNA-binding domain",
        source_database="pfam",
        entry_type="domain",
        integrated_accession="IPR011615",
        source_url="https://www.ebi.ac.uk/interpro/entry/pfam/PF00870/",
        locations=[DomainLocation(start=1, end=8, status="CONTINUOUS")],
        go_terms=[],
        overlaps_mutation=overlap and mutation_position == 2,
    )
    return DomainAnnotation(
        protein_accession="P04637",
        protein_length=10,
        mutation_position=mutation_position,
        entries=[entry],
        mutation_overlaps=[entry] if entry.overlaps_mutation else [],
        entry_count=1,
        location_count=1,
        source_urls=["https://www.ebi.ac.uk/interpro/api/"],
    )


def _kegg():
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


def _graph():
    return build_evidence_graph("TP53", 9606, NETWORK_RECORDS, ENRICHMENT_RECORDS)


def _cellular_processes(annotation, kegg, graph):
    return build_cellular_process_evidence(annotation, kegg, graph)


def _phenotype_literature():
    return PhenotypeLiteratureEvidence(
        gene="TP53",
        process_links=[
            PhenotypeClaimLink(
                process_id="hsa04115",
                process_label="p53 signaling pathway - Homo sapiens (human)",
                gene="TP53",
                support_level="supports",
                pmid="99999",
                title="TP53 regulates p53 signaling pathway.",
                evidence_basis="Direct support: KEGG. matched outcome verb(s): regulates",
                matched_phrases=["regulates"],
                uncertainty="Abstract-level; not causal.",
            )
        ],
        counts_by_level={"supports": 1},
        hypotheses=[
            "Because KEGG links TP53 to p53 signaling pathway (hsa04115), and "
            "PMID 99999 use 'regulates' in the abstract, TP53 may be associated "
            "with p53 signaling pathway. This is an abstract-level hypothesis."
        ],
        query_urls=['"TP53"[Title/Abstract]'],
        boundary="Gated on supports articles only.",
    )


def _full_bundle(mutation="R2H"):
    ann = _annotation()
    analysis = _analysis(mutation)
    structure = _structure(mutation) if mutation else None
    domains = _domains(mutation_position=2 if mutation else None)
    kegg = _kegg()
    graph = _graph()
    cellular = _cellular_processes(ann, kegg, graph)
    phenotype = _phenotype_literature()
    return ann, analysis, structure, domains, kegg, graph, cellular, phenotype


def test_mutation_chain_has_five_steps_with_gated_hypotheses() -> None:
    (ann, analysis, structure, domains, kegg, graph, cellular, phenotype) = _full_bundle("R2H")
    chain = build_reasoning_chain(
        annotation=ann,
        analysis=analysis,
        structure=structure,
        domains=domains,
        kegg=kegg,
        graph=graph,
        cellular_processes=cellular,
        phenotype_literature=phenotype,
        mutation="R2H",
    )

    assert chain.chain_type == "mutation_impact"
    assert chain.gene == "TP53"
    assert chain.mutation == "R2H"
    step_ids = [step.step_id for step in chain.steps]
    assert step_ids == ["mutation", "structure", "function", "pathway", "phenotype"]

    mutation_step = chain.steps[0]
    assert mutation_step.available is True
    assert any("R2H" in fact for fact in mutation_step.evidence_facts)
    # Wild-type R does not match UniProt sequence position 2 (E), so the
    # mismatch must be recorded and downstream "可能影响" hypotheses suppressed.
    assert any("mismatch" in fact.lower() or "does not match" in fact.lower() for fact in mutation_step.evidence_facts)
    assert mutation_step.hypothesis is None

    structure_step = chain.steps[1]
    assert structure_step.available is True
    assert any("pLDDT" in fact or "plddt" in fact.lower() for fact in structure_step.evidence_facts)

    function_step = chain.steps[2]
    assert function_step.available is True
    # No functional-impact hypothesis because the wild type did not match.
    assert function_step.hypothesis is None

    pathway_step = chain.steps[3]
    assert pathway_step.available is True
    assert any("hsa04115" in fact for fact in pathway_step.evidence_facts)

    phenotype_step = chain.steps[4]
    assert phenotype_step.available is True
    assert any("99999" in fact for fact in phenotype_step.evidence_facts)
    assert phenotype_step.hypothesis is not None
    assert "abstract-level hypothesis" in phenotype_step.hypothesis

    assert "not establish" in chain.boundary.lower() or "does not" in chain.boundary.lower()


def test_matching_wild_type_enables_function_and_pathway_hypotheses() -> None:
    # Build an annotation whose sequence has R at position 2 so R2H matches.
    ann = _annotation()
    ann_with_r = ann.model_copy(update={"sequence": "MR" + ann.sequence[2:]})
    analysis = _analysis("R2H")
    structure = _structure("R2H")
    domains = _domains(mutation_position=2, overlap=True)
    kegg = _kegg()
    graph = _graph()
    cellular = _cellular_processes(ann_with_r, kegg, graph)
    phenotype = _phenotype_literature()

    chain = build_reasoning_chain(
        annotation=ann_with_r,
        analysis=analysis,
        structure=structure,
        domains=domains,
        kegg=kegg,
        graph=graph,
        cellular_processes=cellular,
        phenotype_literature=phenotype,
        mutation="R2H",
    )

    mutation_step = chain.steps[0]
    assert mutation_step.hypothesis is not None
    # Conservative wording: candidate site for follow-up, not "可能影响".
    assert "候选" in mutation_step.hypothesis or "后续功能影响评估" in mutation_step.hypothesis
    assert "可能影响局部结构或功能" not in mutation_step.hypothesis
    assert "可能影响蛋白功能" not in mutation_step.hypothesis
    function_step = chain.steps[2]
    assert function_step.hypothesis is not None
    assert "候选" in function_step.hypothesis or "功能影响评估" in function_step.hypothesis
    assert "可能影响蛋白功能" not in function_step.hypothesis
    pathway_step = chain.steps[3]
    assert pathway_step.hypothesis is not None
    assert "可能" in pathway_step.hypothesis


def test_no_mutation_chain_is_protein_function_type_with_four_steps() -> None:
    (ann, analysis, _structure_unused, domains, kegg, graph, cellular, phenotype) = _full_bundle(None)
    chain = build_reasoning_chain(
        annotation=ann,
        analysis=analysis,
        structure=None,
        domains=domains,
        kegg=kegg,
        graph=graph,
        cellular_processes=cellular,
        phenotype_literature=phenotype,
        mutation=None,
    )

    assert chain.chain_type == "protein_function"
    assert chain.mutation is None
    step_ids = [step.step_id for step in chain.steps]
    assert step_ids == ["structure", "function", "pathway", "phenotype"]
    # No "可能影响" hypotheses when there is no mutation.
    for step in chain.steps:
        if step.step_id in ("function", "pathway"):
            assert step.hypothesis is None


def test_missing_layers_recorded_and_step_marked_unavailable() -> None:
    (ann, analysis, _structure_unused, domains, _kegg_unused, graph, cellular, phenotype) = _full_bundle("R2H")
    chain = build_reasoning_chain(
        annotation=ann,
        analysis=analysis,
        structure=None,
        domains=domains,
        kegg=None,
        graph=graph,
        cellular_processes=cellular,
        phenotype_literature=phenotype,
        mutation="R2H",
    )

    assert "structure" in chain.missing_layers
    assert "kegg" in chain.missing_layers
    structure_step = next(step for step in chain.steps if step.step_id == "structure")
    assert structure_step.available is False
    assert structure_step.confidence == "insufficient"


def test_malformed_mutation_marks_mutation_step_unavailable() -> None:
    (ann, analysis, structure, domains, kegg, graph, cellular, phenotype) = _full_bundle("R2H")
    chain = build_reasoning_chain(
        annotation=ann,
        analysis=analysis,
        structure=structure,
        domains=domains,
        kegg=kegg,
        graph=graph,
        cellular_processes=cellular,
        phenotype_literature=phenotype,
        mutation="NOTAMUTATION",
    )

    assert chain.chain_type == "mutation_impact"
    mutation_step = chain.steps[0]
    assert mutation_step.available is False
    assert mutation_step.confidence == "insufficient"
    assert any("invalid" in fact.lower() or "malformed" in fact.lower() for fact in mutation_step.evidence_facts)


def test_phenotype_step_handles_no_supports_articles() -> None:
    (ann, analysis, structure, domains, kegg, graph, cellular, _phenotype_unused) = _full_bundle("R2H")
    empty_phenotype = PhenotypeLiteratureEvidence(
        gene="TP53",
        process_links=[],
        counts_by_level={},
        hypotheses=[],
        query_urls=[],
        boundary="No supports articles.",
    )
    chain = build_reasoning_chain(
        annotation=ann,
        analysis=analysis,
        structure=structure,
        domains=domains,
        kegg=kegg,
        graph=graph,
        cellular_processes=cellular,
        phenotype_literature=empty_phenotype,
        mutation="R2H",
    )

    phenotype_step = next(step for step in chain.steps if step.step_id == "phenotype")
    assert phenotype_step.available is True
    assert phenotype_step.hypothesis is None
    assert any("no" in fact.lower() and "support" in fact.lower() for fact in phenotype_step.evidence_facts)

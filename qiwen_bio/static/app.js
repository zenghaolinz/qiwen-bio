const form = document.querySelector("#analysis-form");
const sequenceInput = document.querySelector("#sequence");
const lengthLabel = document.querySelector("#sequence-length");
const lookupButton = document.querySelector("#lookup-button");

function sequenceLength() {
  return sequenceInput.value.replace(/\s/g, "").length;
}
sequenceInput.addEventListener("input", () => { lengthLabel.textContent = `${sequenceLength()} aa`; });

form.addEventListener("submit", async (event) => {
  event.preventDefault();
  await requestAnalysis("/api/v1/analyze", {
    sequence: sequenceInput.value,
    name: document.querySelector("#name").value || null,
    mutation: document.querySelector("#mutation").value || null,
  });
});

lookupButton.addEventListener("click", async () => {
  const identifier = document.querySelector("#identifier").value.trim();
  if (!identifier) return showError("请输入基因名或 UniProt ID");
  await requestAnalysis("/api/v1/analyze/uniprot", {
    identifier,
    organism_id: 9606,
    mutation: document.querySelector("#mutation").value || null,
  });
});

async function requestAnalysis(url, payload) {
  const empty = document.querySelector("#empty-state");
  const loading = document.querySelector("#loading");
  const error = document.querySelector("#error");
  const results = document.querySelector("#results");
  empty.classList.add("hidden"); error.classList.add("hidden"); results.classList.add("hidden"); loading.classList.remove("hidden");

  try {
    const response = await fetch(url, {
      method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(payload),
    });
    const data = await response.json();
    if (!response.ok) throw new Error(data.detail?.[0]?.msg || "分析请求失败");
    const analysis = data.analysis || data;
    renderAnnotation(data.annotation || null);
    document.querySelector("#metric-length").textContent = analysis.features.length;
    document.querySelector("#metric-weight").textContent = analysis.features.molecular_weight_da.toLocaleString();
    document.querySelector("#metric-hydrophobic").textContent = `${(analysis.features.hydrophobic_fraction * 100).toFixed(1)}%`;
    document.querySelector("#metric-charge").textContent = `${analysis.features.net_charge_proxy >= 0 ? "+" : ""}${analysis.features.net_charge_proxy}`;
    document.querySelector("#prediction-label").textContent = analysis.prediction.label;
    document.querySelector("#prediction-score").textContent = analysis.prediction.score.toFixed(3);
    document.querySelector("#evidence-list").innerHTML = analysis.evidence_chain.map(item =>
      `<li><p>${escapeHtml(item.claim)}</p><small>${escapeHtml(item.source)} · ${item.confidence} confidence</small></li>`
    ).join("");
    document.querySelector("#report").textContent = analysis.report_markdown;
    results.classList.remove("hidden");
    if (data.annotation?.alphafold_url) {
      await renderStructure(data.annotation.accession, payload.mutation);
    } else {
      document.querySelector("#structure").classList.add("hidden");
    }
  } catch (err) {
    error.textContent = err.message; error.classList.remove("hidden");
  } finally { loading.classList.add("hidden"); }
}

async function renderStructure(accession, mutation) {
  const panel = document.querySelector("#structure");
  try {
    const response = await fetch("/api/v1/structure/alphafold", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ accession, mutation: mutation || null }),
    });
    const data = await response.json();
    if (!response.ok) throw new Error(typeof data.detail === "string" ? data.detail : "Structure analysis failed");
    const structure = data.structure;
    document.querySelector("#structure-plddt").textContent = structure.mean_plddt.toFixed(2);
    const distribution = structure.confidence_distribution;
    document.querySelector("#confidence-bar").innerHTML = ["very_high", "confident", "low", "very_low"]
      .map(key => `<span class="${key.replace("_", "-")}" style="width:${distribution[key] * 100}%"></span>`)
      .join("");
    const mutationLine = document.querySelector("#mutation-confidence");
    mutationLine.classList.toggle("hidden", !structure.mutation_site);
    if (structure.mutation_site) {
      const site = structure.mutation_site;
      mutationLine.textContent = `${site.wild_type}${site.position}${site.mutant}: pLDDT ${site.plddt.toFixed(1)} (${site.confidence})`;
    }
    document.querySelector("#structure-note").textContent = data.interpretation;
  } catch (err) {
    document.querySelector("#structure-plddt").textContent = "Unavailable";
    document.querySelector("#confidence-bar").innerHTML = "";
    document.querySelector("#structure-note").textContent = err.message;
  }
  panel.classList.remove("hidden");
}

function renderAnnotation(annotation) {
  const panel = document.querySelector("#annotation");
  if (!annotation) return panel.classList.add("hidden");
  document.querySelector("#annotation-name").textContent = annotation.protein_name;
  document.querySelector("#annotation-accession").textContent = annotation.accession;
  document.querySelector("#annotation-organism").textContent = annotation.organism;
  document.querySelector("#annotation-go").textContent = annotation.go_terms.length;
  document.querySelector("#annotation-function").textContent = annotation.functions[0] || "UniProt 暂无功能描述。";
  document.querySelector("#annotation-source").href = annotation.source_url;
  const alphafold = document.querySelector("#annotation-alphafold");
  alphafold.classList.toggle("hidden", !annotation.alphafold_url);
  if (annotation.alphafold_url) alphafold.href = annotation.alphafold_url;
  panel.classList.remove("hidden");
}

function showError(message) {
  const error = document.querySelector("#error");
  document.querySelector("#empty-state").classList.add("hidden");
  error.textContent = message; error.classList.remove("hidden");
}

function escapeHtml(value) {
  const node = document.createElement("span"); node.textContent = value; return node.innerHTML;
}

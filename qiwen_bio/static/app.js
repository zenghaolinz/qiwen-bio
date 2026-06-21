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
    if (data.annotation) {
      const enrichments = [renderEvidenceGraph(data.annotation.gene_names[0] || data.annotation.accession)];
      if (data.annotation.alphafold_url) enrichments.push(renderStructure(data.annotation.accession, payload.mutation));
      else document.querySelector("#structure").classList.add("hidden");
      await Promise.all(enrichments);
    } else {
      document.querySelector("#structure").classList.add("hidden");
      document.querySelector("#evidence-graph").classList.add("hidden");
    }
  } catch (err) {
    error.textContent = err.message; error.classList.remove("hidden");
  } finally { loading.classList.add("hidden"); }
}

async function renderStructure(accession, mutation) {
  const panel = document.querySelector("#structure");
  panel.classList.remove("hidden");
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
    renderBackbone(structure.coordinates, structure.mutation_site?.position || null);
    renderContactMap(structure.coordinates, structure.contact_map);
    renderNeighborhood(structure.mutation_neighborhood);
  } catch (err) {
    document.querySelector("#structure-plddt").textContent = "Unavailable";
    document.querySelector("#confidence-bar").innerHTML = "";
    document.querySelector("#structure-note").textContent = err.message;
  }
}

const viewerState = { coordinates: [], mutationPosition: null, angleX: -0.25, angleY: 0.45 };

function plddtColor(value) {
  if (value >= 90) return "#2f6de0";
  if (value >= 70) return "#65c7df";
  if (value >= 50) return "#f0d451";
  return "#ef762f";
}

function renderBackbone(coordinates, mutationPosition) {
  viewerState.coordinates = coordinates;
  viewerState.mutationPosition = mutationPosition;
  drawBackbone();
}

function drawBackbone() {
  const canvas = document.querySelector("#structure-viewer");
  const { ctx, width, height } = prepareCanvas(canvas);
  if (!viewerState.coordinates.length) return;
  const center = viewerState.coordinates.reduce((sum, point) => ({
    x: sum.x + point.x / viewerState.coordinates.length,
    y: sum.y + point.y / viewerState.coordinates.length,
    z: sum.z + point.z / viewerState.coordinates.length,
  }), { x: 0, y: 0, z: 0 });
  const sinY = Math.sin(viewerState.angleY), cosY = Math.cos(viewerState.angleY);
  const sinX = Math.sin(viewerState.angleX), cosX = Math.cos(viewerState.angleX);
  const rotated = viewerState.coordinates.map(point => {
    const x = point.x - center.x, y = point.y - center.y, z = point.z - center.z;
    const x1 = x * cosY + z * sinY, z1 = -x * sinY + z * cosY;
    return { ...point, rx: x1, ry: y * cosX - z1 * sinX, rz: y * sinX + z1 * cosX };
  });
  const extent = Math.max(...rotated.flatMap(point => [Math.abs(point.rx), Math.abs(point.ry)]), 1);
  const scale = Math.min(width, height) * 0.42 / extent;
  const projected = rotated.map(point => ({ ...point, sx: width / 2 + point.rx * scale, sy: height / 2 + point.ry * scale }));
  ctx.lineWidth = Math.max(1, width / 700); ctx.strokeStyle = "rgba(174, 204, 214, .38)"; ctx.beginPath();
  projected.forEach((point, index) => { if (index === 0) ctx.moveTo(point.sx, point.sy); else ctx.lineTo(point.sx, point.sy); }); ctx.stroke();
  [...projected].sort((a, b) => a.rz - b.rz).forEach(point => {
    ctx.beginPath(); ctx.fillStyle = point.position === viewerState.mutationPosition ? "#d7f45b" : plddtColor(point.plddt);
    ctx.arc(point.sx, point.sy, point.position === viewerState.mutationPosition ? 5 : 2.2, 0, Math.PI * 2); ctx.fill();
  });
}

function renderContactMap(coordinates, contactMap) {
  const canvas = document.querySelector("#contact-map");
  const { ctx, width, height } = prepareCanvas(canvas);
  const positions = coordinates.map(item => item.position);
  const min = Math.min(...positions), max = Math.max(...positions), span = Math.max(max - min, 1);
  ctx.strokeStyle = "rgba(174, 204, 214, .18)"; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(4, 4); ctx.lineTo(width - 4, height - 4); ctx.stroke();
  ctx.fillStyle = "rgba(101, 214, 239, .88)";
  contactMap.contacts.forEach(contact => {
    const xA = ((contact.residue_a - min) / span) * (width - 8) + 4;
    const yA = ((contact.residue_a - min) / span) * (height - 8) + 4;
    const xB = ((contact.residue_b - min) / span) * (width - 8) + 4;
    const yB = ((contact.residue_b - min) / span) * (height - 8) + 4;
    ctx.fillRect(xA, yB, 2.2, 2.2); ctx.fillRect(xB, yA, 2.2, 2.2);
  });
  document.querySelector("#contact-caption").textContent = `Contact map · ${contactMap.total_contacts.toLocaleString()} pairs${contactMap.truncated ? " · truncated" : ""}`;
}

function renderNeighborhood(neighbors) {
  const panel = document.querySelector("#mutation-neighborhood");
  panel.classList.toggle("hidden", !neighbors.length);
  document.querySelector("#neighbor-list").innerHTML = neighbors.map(item =>
    `<span class="neighbor-chip"><b>${item.amino_acid}${item.position}</b> · ${item.distance_angstrom.toFixed(2)} Å</span>`
  ).join("");
}

async function renderEvidenceGraph(identifier) {
  const panel = document.querySelector("#evidence-graph");
  panel.classList.remove("hidden");
  document.querySelector("#graph-title").textContent = identifier;
  document.querySelector("#graph-summary").textContent = "Loading STRING evidence...";
  try {
    const response = await fetch("/api/v1/graph/string", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ identifier, organism_id: 9606, limit: 8, required_score: 700 }),
    });
    const graph = await response.json();
    if (!response.ok) throw new Error(typeof graph.detail === "string" ? graph.detail : "STRING graph failed");
    const interactions = graph.edges.filter(edge => edge.type === "interacts_with").length;
    const terms = graph.nodes.filter(node => node.type !== "protein");
    document.querySelector("#graph-summary").textContent = `${graph.nodes.length} nodes · ${interactions} interactions · ${terms.length} terms`;
    drawEvidenceGraph(graph);
    document.querySelector("#graph-terms").innerHTML = terms.slice(0, 6).map(node =>
      `<div class="graph-term"><b>${escapeHtml(node.external_id)}</b>${escapeHtml(node.label)} · FDR ${formatScientific(node.fdr)}</div>`
    ).join("");
  } catch (err) {
    document.querySelector("#graph-summary").textContent = "Unavailable";
    document.querySelector("#graph-terms").innerHTML = `<div class="graph-term">${escapeHtml(err.message)}</div>`;
    prepareCanvas(document.querySelector("#graph-canvas"));
  }
}

function drawEvidenceGraph(graph) {
  const canvas = document.querySelector("#graph-canvas");
  const { ctx, width, height } = prepareCanvas(canvas);
  const proteins = graph.nodes.filter(node => node.type === "protein");
  const terms = graph.nodes.filter(node => node.type !== "protein");
  const seed = proteins.find(node => node.label.toUpperCase() === graph.seed.toUpperCase()) || proteins[0];
  const others = proteins.filter(node => node.id !== seed?.id);
  const positions = new Map();
  if (seed) positions.set(seed.id, { x: width / 2, y: height / 2 });
  placeRing(others, width / 2, height / 2, Math.min(width, height) * 0.23, positions, -Math.PI / 2);
  placeRing(terms, width / 2, height / 2, Math.min(width, height) * 0.43, positions, -Math.PI / 2);

  graph.edges.forEach(edge => {
    const source = positions.get(edge.source), target = positions.get(edge.target);
    if (!source || !target) return;
    ctx.beginPath(); ctx.moveTo(source.x, source.y); ctx.lineTo(target.x, target.y);
    ctx.strokeStyle = edge.type === "interacts_with" ? "rgba(45, 71, 62, .42)" : "rgba(179, 138, 29, .25)";
    ctx.lineWidth = edge.type === "interacts_with" ? 1.2 + (edge.score || 0) * 1.8 : 1;
    ctx.setLineDash(edge.type === "annotated_to" ? [3, 4] : []); ctx.stroke();
  });
  ctx.setLineDash([]);
  graph.nodes.forEach(node => {
    const point = positions.get(node.id); if (!point) return;
    const isSeed = node.id === seed?.id;
    ctx.beginPath();
    if (node.type === "pathway") ctx.rect(point.x - 6, point.y - 6, 12, 12);
    else if (node.type === "process") { ctx.moveTo(point.x, point.y - 7); ctx.lineTo(point.x + 7, point.y); ctx.lineTo(point.x, point.y + 7); ctx.lineTo(point.x - 7, point.y); ctx.closePath(); }
    else ctx.arc(point.x, point.y, isSeed ? 10 : 7, 0, Math.PI * 2);
    ctx.fillStyle = node.type === "protein" ? (isSeed ? "#d7f45b" : "#1d604b") : node.type === "pathway" ? "#d6a824" : "#3d88a8";
    ctx.fill(); ctx.strokeStyle = "rgba(20,33,29,.6)"; ctx.lineWidth = 1; ctx.stroke();
    ctx.fillStyle = "#26332f"; ctx.font = node.type === "protein" ? "bold 10px Candara" : "9px Candara";
    ctx.textAlign = "center"; ctx.fillText(shortLabel(node.label, node.type === "protein" ? 14 : 20), point.x, point.y + (node.type === "protein" ? 20 : 18));
  });
}

function placeRing(nodes, centerX, centerY, radius, positions, offset) {
  nodes.forEach((node, index) => {
    const angle = offset + (Math.PI * 2 * index) / Math.max(nodes.length, 1);
    positions.set(node.id, { x: centerX + Math.cos(angle) * radius, y: centerY + Math.sin(angle) * radius });
  });
}

function shortLabel(value, length) { return value.length > length ? `${value.slice(0, length - 1)}…` : value; }
function formatScientific(value) { return value == null ? "n/a" : Number(value).toExponential(2); }

function prepareCanvas(canvas) {
  const ratio = window.devicePixelRatio || 1, width = canvas.clientWidth, height = canvas.clientHeight;
  canvas.width = Math.round(width * ratio); canvas.height = Math.round(height * ratio);
  const ctx = canvas.getContext("2d"); ctx.scale(ratio, ratio); ctx.clearRect(0, 0, width, height);
  return { ctx, width, height };
}

const structureCanvas = document.querySelector("#structure-viewer");
let dragPoint = null;
structureCanvas.addEventListener("pointerdown", event => { dragPoint = { x: event.clientX, y: event.clientY }; structureCanvas.setPointerCapture(event.pointerId); });
structureCanvas.addEventListener("pointermove", event => {
  if (!dragPoint) return;
  viewerState.angleY += (event.clientX - dragPoint.x) * 0.01; viewerState.angleX += (event.clientY - dragPoint.y) * 0.01;
  dragPoint = { x: event.clientX, y: event.clientY }; drawBackbone();
});
structureCanvas.addEventListener("pointerup", () => { dragPoint = null; });
window.addEventListener("resize", () => { if (viewerState.coordinates.length) drawBackbone(); });

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

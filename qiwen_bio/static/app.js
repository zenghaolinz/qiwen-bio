const form = document.querySelector("#analysis-form");
const sequenceInput = document.querySelector("#sequence");
const lengthLabel = document.querySelector("#sequence-length");

function sequenceLength() {
  return sequenceInput.value.replace(/\s/g, "").length;
}
sequenceInput.addEventListener("input", () => { lengthLabel.textContent = `${sequenceLength()} aa`; });

form.addEventListener("submit", async (event) => {
  event.preventDefault();
  const empty = document.querySelector("#empty-state");
  const loading = document.querySelector("#loading");
  const error = document.querySelector("#error");
  const results = document.querySelector("#results");
  empty.classList.add("hidden"); error.classList.add("hidden"); results.classList.add("hidden"); loading.classList.remove("hidden");

  const payload = {
    sequence: sequenceInput.value,
    name: document.querySelector("#name").value || null,
    mutation: document.querySelector("#mutation").value || null,
  };
  try {
    const response = await fetch("/api/v1/analyze", {
      method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(payload),
    });
    const data = await response.json();
    if (!response.ok) throw new Error(data.detail?.[0]?.msg || "分析请求失败");
    document.querySelector("#metric-length").textContent = data.features.length;
    document.querySelector("#metric-weight").textContent = data.features.molecular_weight_da.toLocaleString();
    document.querySelector("#metric-hydrophobic").textContent = `${(data.features.hydrophobic_fraction * 100).toFixed(1)}%`;
    document.querySelector("#metric-charge").textContent = `${data.features.net_charge_proxy >= 0 ? "+" : ""}${data.features.net_charge_proxy}`;
    document.querySelector("#prediction-label").textContent = data.prediction.label;
    document.querySelector("#prediction-score").textContent = data.prediction.score.toFixed(3);
    document.querySelector("#evidence-list").innerHTML = data.evidence_chain.map(item =>
      `<li><p>${escapeHtml(item.claim)}</p><small>${escapeHtml(item.source)} · ${item.confidence} confidence</small></li>`
    ).join("");
    document.querySelector("#report").textContent = data.report_markdown;
    results.classList.remove("hidden");
  } catch (err) {
    error.textContent = err.message; error.classList.remove("hidden");
  } finally { loading.classList.add("hidden"); }
});

function escapeHtml(value) {
  const node = document.createElement("span"); node.textContent = value; return node.innerHTML;
}


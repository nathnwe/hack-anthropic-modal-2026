import {
  type GeneRecord,
  type Mode,
  modeLabels,
  isMode,
  formatNumber as fmt,
  regionLabel,
  isPlaceholder,
  contactsFor,
  contactStrength,
  sortElements,
  selectedElementFor,
  candidatesFor,
  evidenceLink,
  escapeHTML as esc,
} from "../lib/records";

const get = <T extends HTMLElement = HTMLElement>(id: string) =>
  document.getElementById(id) as T;
const app = get("results-app");
const manifest: { symbol: string; file: string; direction: string }[] =
  JSON.parse(app.dataset.manifest || "[]");
const base = app.dataset.base || "";
let record: GeneRecord;
let mode: Mode;
let activeElement: string | null = null;
let genomeViewer: import("../lib/genome-viewer").GenomeViewer | null = null;
let viewerLoading = false;
let viewerFailed = false;

function syncGenomeViewer() {
  const selected =
    record.locus.cres.find((c) => c.id === activeElement) ?? null;
  if (genomeViewer) {
    genomeViewer.update(record, selected);
    return;
  }
  if (viewerLoading || viewerFailed) return;
  viewerLoading = true;
  void import("../lib/genome-viewer")
    .then(({ GenomeViewer }) => {
      genomeViewer = new GenomeViewer();
      genomeViewer.update(
        record,
        record.locus.cres.find((c) => c.id === activeElement) ?? null,
      );
    })
    .catch(() => {
      viewerFailed = true;
      get("dna-loading").hidden = true;
      get("dna-fallback").hidden = false;
      get("dna-fallback-message").textContent =
        "This browser could not start the 3D renderer. The contact map remains available.";
      get("genome-actions").hidden = true;
      const switchView = (threeD: boolean) => {
        get("dna-panel").hidden = !threeD;
        get("contact-panel").hidden = threeD;
        get("view-3d").setAttribute("aria-pressed", String(threeD));
        get("view-contacts").setAttribute("aria-pressed", String(!threeD));
        drawMap();
      };
      get("view-3d").addEventListener("click", () => switchView(true));
      get("view-contacts").addEventListener("click", () => switchView(false));
      get("dna-fallback-map").addEventListener("click", () =>
        switchView(false),
      );
    })
    .finally(() => {
      viewerLoading = false;
    });
}
window.addEventListener("pagehide", (event) => {
  if (!event.persisted) genomeViewer?.dispose();
});
if (import.meta.hot) import.meta.hot.dispose(() => genomeViewer?.dispose());

function evidenceMarkup() {
  return (
    record.gene.evidence
      .map((source) => {
        const link = evidenceLink(source);
        return `<li>${isPlaceholder(source) ? '<span class="tag placeholder-tag">Placeholder source</span> ' : ""}${link ? `<a href="${esc(link)}" target="_blank" rel="noopener noreferrer">${esc(source)} ↗</a>` : esc(source)}</li>`;
      })
      .join("") || "<li>No sources supplied.</li>"
  );
}

function renderElements() {
  const criterion = get<HTMLSelectElement>("rank-select").value;
  get("rank-explanation").textContent =
    criterion === "re2g"
      ? "Ordered by supplied rE2G score. This is not enhancer potency or a probability that stapling will work."
      : "Ordered by supplied promoter-contact strength. Missing contact evidence sorts last. These values do not predict expression change.";
  const list = sortElements(record, criterion);
  activeElement = selectedElementFor(record, criterion)?.id ?? null;
  get("element-list").innerHTML =
    list
      .map((cre, i) => {
        const strength = contactStrength(record, cre);
        const matches = contactsFor(record, cre);
        return `<details class="element" data-element="${esc(cre.id)}" data-selected="${activeElement === cre.id}">
      <summary><span class="element-index mono">${String(i + 1).padStart(2, "0")}</span><span class="element-name"><strong>${esc(cre.id)}</strong><span class="small muted">${esc(cre.type)} · ${fmt(cre.start)}–${fmt(cre.end)}</span></span><span class="element-metric"><span class="small muted">${criterion === "re2g" ? "rE2G score" : "Contact strength"}</span><span class="mono">${criterion === "re2g" ? esc(cre.score_rE2G) : strength === null ? "Not recorded" : esc(strength)}</span></span><span class="expand-icon" aria-hidden="true">+</span></summary>
      <div class="element-detail"><div><p class="eyebrow">Recorded criteria ${record.illustrative ? "· placeholder" : ""}</p><dl class="criteria"><div><dt>rE2G score</dt><dd>${esc(cre.score_rE2G)} <span class="small muted">raw model score</span></dd></div><div><dt>Promoter contact</dt><dd>${strength === null ? "Not recorded" : esc(strength) + " · supplied strength"}</dd></div><div><dt>Matched contacts</dt><dd>${matches.length}</dd></div><div><dt>Potency / achievability</dt><dd>Not supplied</dd></div><div><dt>Element composite / tier</dt><dd>Not supplied (staples only)</dd></div></dl></div><div><p class="eyebrow">Gene-level literature</p><ul class="evidence-list">${evidenceMarkup()}</ul><p class="small muted">Context for the gene; not element-specific validation.</p></div></div>
    </details>`;
      })
      .join("") ||
    '<p class="empty-inline">No regulatory elements are supplied for this locus.</p>';
  get("element-list")
    .querySelectorAll<HTMLDetailsElement>("details")
    .forEach((detail) => {
      detail.querySelector("summary")!.addEventListener("click", () => {
        activeElement = detail.dataset.element!;
        get("element-list")
          .querySelectorAll<HTMLDetailsElement>("details")
          .forEach((other) => {
            other.dataset.selected = String(other === detail);
            if (other !== detail) other.open = false;
          });
        drawMap();
      });
    });
}

function drawMap() {
  syncGenomeViewer();
  const { tad, cres, contacts } = record.locus;
  const svg = document.getElementById("contact-map")!;
  const selected = cres.find((c) => c.id === activeElement);
  const selectedContacts = selected ? contactsFor(record, selected) : [];
  get("viewer-element").textContent = selected?.id ?? "No element supplied";
  get("viewer-context").textContent = selected
    ? `${selected.type} · ${regionLabel({ chrom: tad.chrom, start: selected.start, end: selected.end })}`
    : "This record has no regulatory elements to display.";
  svg.setAttribute(
    "aria-label",
    `${record.gene.symbol} contact schematic${selected ? ` for ${selected.id}` : ""}${record.illustrative ? " — placeholder data" : ""}`,
  );
  const start = tad.start,
    end = tad.end;
  if (end <= start) {
    svg.innerHTML =
      '<text x="40" y="150" fill="currentColor">The supplied TAD interval cannot be plotted.</text>';
    get("map-note").textContent =
      "Invalid or empty TAD interval in the record.";
    return;
  }
  const width = Math.max(280, svg.parentElement!.clientWidth || 640);
  const plotRight = width - 32;
  svg.setAttribute("viewBox", `0 0 ${width} 305`);
  const x = (n: number) => 32 + ((n - start) / (end - start)) * (width - 64);
  const inside = (n: number) => n >= start && n <= end;
  const shown = contacts.filter((c) => inside(c.a) && inside(c.b));
  const plotCREs = cres.filter((c) => c.end >= start && c.start <= end);
  const promoterInside = inside(record.gene.tss);
  let markup = `<title>${esc(record.gene.symbol)} contact map${selected ? ` — ${esc(selected.id)}` : ""}${record.illustrative ? " — placeholder coordinates" : ""}</title><desc>Linear positions across the supplied TAD. Square: gene promoter. Circles: regulatory elements. Arcs: supplied contacts, with strength in each arc title. Arc emphasis identifies the selected element, not likelihood.</desc>`;
  markup += `<line x1="32" y1="209" x2="${plotRight}" y2="209" stroke="currentColor" opacity=".28"/>`;
  const intervals = width < 480 ? 2 : 4;
  for (let i = 0; i <= intervals; i++) {
    const p = start + ((end - start) * i) / intervals;
    markup += `<line x1="${x(p)}" y1="203" x2="${x(p)}" y2="215" stroke="currentColor" opacity=".25"/><text x="${x(p)}" y="286" text-anchor="${i === 0 ? "start" : i === intervals ? "end" : "middle"}" class="axis-label">${fmt(Math.round(p))}</text>`;
  }
  shown.forEach((c) => {
    const a = x(c.a),
      b = x(c.b),
      peak = 209 - Math.min(190, Math.max(45, Math.abs(b - a) * 0.48));
    const highlighted = selectedContacts.includes(c);
    markup += `<path d="M${a} 209 C${a} ${peak} ${b} ${peak} ${b} 209" fill="none" stroke="#A54C38" stroke-width="${highlighted ? 4 : 1.8}" opacity="${selected && !highlighted ? 0.12 : 0.9}"><title>Contact: ${fmt(c.a)} ↔ ${fmt(c.b)}; strength ${esc(c.strength)}${record.illustrative ? " (placeholder)" : ""}</title></path>`;
  });
  plotCREs.forEach((cre) => {
    const midpoint = Math.max(start, Math.min(end, (cre.start + cre.end) / 2));
    const pos = x(midpoint);
    markup += `<g><title>${esc(cre.id)} · ${esc(cre.type)} · ${fmt(cre.start)}–${fmt(cre.end)}</title><circle cx="${pos}" cy="209" r="${activeElement === cre.id ? 7 : 4}" fill="#A54C38" opacity="${activeElement === cre.id ? 1 : 0.3}"/>${activeElement === cre.id ? `<text x="${pos}" y="256" text-anchor="${pos > width - 75 ? "end" : pos < 75 ? "start" : "middle"}" class="cre-label" font-weight="600">${esc(cre.id)}</text>` : ""}</g>`;
  });
  if (promoterInside)
    markup += `<rect x="${x(record.gene.tss) - 7}" y="202" width="14" height="14" fill="#2A6FA8"/><text x="${x(record.gene.tss)}" y="239" text-anchor="middle" class="promoter-label">${esc(record.gene.symbol)} promoter</text>`;
  if (!shown.length)
    markup += `<text x="${width / 2}" y="100" text-anchor="middle" class="axis-label">No in-domain contacts supplied</text>`;
  svg.innerHTML = markup;
  const omitted = contacts.length - shown.length;
  get("map-note").textContent =
    `${record.illustrative ? "Illustrative coordinates. " : ""}Schematic contacts, not a predicted 3D structure. ${omitted ? `${omitted} out-of-domain contact(s) omitted. ` : ""}${!promoterInside ? "The supplied promoter is outside this TAD. " : ""}${selected && !selectedContacts.length ? "No promoter contact is recorded for this element." : ""}`;
}

function renderStaples() {
  const staples = candidatesFor(record, mode);
  get("staple-count").textContent =
    `${staples.length} matching ${staples.length === 1 ? "design" : "designs"}`;
  if (!staples.length) {
    get("staples-content").innerHTML =
      `<div class="empty-inline"><h3>No ${modeLabels[mode].toLowerCase()} staple is supplied.</h3><p>The record has no candidate for this direction. Recorded contacts remain available to explore.</p></div>`;
    return;
  }
  get("staples-content").innerHTML =
    `<div class="table-scroll" role="region" aria-label="Candidate staple anchors and risks, scroll horizontally on small screens" tabindex="0"><table class="staple-table"><caption class="sr-only">Candidate staples for ${esc(record.gene.symbol)} — ${modeLabels[mode]}${record.illustrative ? " (placeholder data)" : ""}</caption><thead><tr><th scope="col">Candidate</th><th scope="col">Anchor pair</th><th scope="col">Mode / change</th><th scope="col">Tier</th><th scope="col">Risk flags</th></tr></thead><tbody>${staples
      .map((staple) => {
        const ranking = record.ranking.find((r) => r.staple_id === staple.id);
        const risks = record.risks
          .filter((r) => r.staple_id === staple.id)
          .flatMap((r) => r.flags);
        return `<tr><th scope="row" class="mono">${esc(staple.id)}</th><td class="anchor-pair"><span><i>A</i>${esc(regionLabel(staple.anchor_a))}</span><span><i>B</i>${esc(regionLabel(staple.anchor_b))}</span></td><td>${modeLabels[staple.mode]}<small>${esc(staple.predicted_delta.replaceAll("_", " "))}</small></td><td>${ranking ? `<span class="tier tier-${ranking.tier}">${ranking.tier}</span>` : "Not supplied"}</td><td>${risks.length ? risks.map((f) => `<div class="risk-flag"><strong>${esc(f.severity)} · ${esc(f.type)}</strong><span>${esc(f.detail)}${f.gene ? ` · ${esc(f.gene)}` : ""}</span></div>`).join("") : "No flags supplied; safety not established."}</td></tr>`;
      })
      .join("")}</tbody></table></div><div class="design-notes">${staples
      .map((staple) => {
        const r = record.ranking.find((r) => r.staple_id === staple.id);
        return `<div class="design-note"><div><p class="eyebrow">${esc(staple.id)} / Rationale</p>${r ? `<ul>${r.rationale.map((s) => `<li>${esc(s)}</li>`).join("")}</ul>` : "<p>No rationale supplied.</p>"}</div><div><p class="eyebrow">Recommended method ${record.illustrative ? "· placeholder" : ""}</p><p>${r ? esc(r.method_recommendation) : "No method supplied."}</p></div></div>`;
      })
      .join("")}</div>`;
}

function updateMode() {
  get("intent-title").textContent =
    `${mode === "up" ? "↑" : mode === "down" ? "↓" : "⊘"} ${modeLabels[mode]}`;
  get<HTMLSelectElement>("intent-select").value = mode;
  renderStaples();
  get("results-status").textContent =
    `${record.gene.symbol}, ${modeLabels[mode]}. ${candidatesFor(record, mode).length} matching candidates.${record.illustrative ? " Placeholder data." : ""}`;
}

function renderRecord() {
  get("gene-title").textContent = record.gene.symbol;
  document.title = `${record.gene.symbol} — Results — Rewire Bio`;
  get("disease-line").textContent = record.gene.disease;
  get("record-notice").innerHTML = record.illustrative
    ? '<span class="notice-icon" aria-hidden="true">◇</span><div><strong>Placeholder record — not a biological finding.</strong><p>Disease annotations, cell type, coordinates, scores, tiers and recommendations are illustrative. In silico and experimentally unvalidated.</p></div>'
    : `<span class="notice-icon" aria-hidden="true">◇</span><div><strong>In silico record — experimentally unvalidated.</strong><p>${record.illustrative === undefined ? "Illustrative status is unspecified. Review source provenance before interpreting this record." : "Review the sources and limitations before interpreting any proposed intervention."}</p></div>`;
  get("tad-label").textContent = regionLabel(record.locus.tad);
  get("map-label").textContent = record.illustrative
    ? "Placeholder schematic"
    : "Contact schematic";
  get("locus-facts").innerHTML =
    `<div><dt>Cell context</dt><dd>${isPlaceholder(record.locus.cell_type) ? "Not specified (placeholder)" : esc(record.locus.cell_type)}</dd></div><div><dt>Regulatory elements</dt><dd>${record.locus.cres.length}</dd></div><div><dt>Recorded contacts</dt><dd>${record.locus.contacts.length}</dd></div><div><dt>Data status</dt><dd>${record.illustrative ? "Placeholder" : record.illustrative === false ? "In silico" : "Provenance unspecified"}</dd></div>`;
  renderElements();
  drawMap();
  updateMode();
  get("loading-state").hidden = true;
  get("record-view").hidden = false;
}

async function loadRecord() {
  const params = new URLSearchParams(location.search);
  const symbol = (params.get("gene") || manifest[0]?.symbol || "")
    .trim()
    .toUpperCase();
  const entry = manifest.find((item) => item.symbol === symbol);
  try {
    if (!entry)
      throw new Error(
        `No record for “${symbol}” in this build. Available records: ${manifest.map((r) => r.symbol).join(", ")}.`,
      );
    const requested = params.get("intent");
    if (requested !== null && !isMode(requested))
      throw new Error(
        "The requested regulation direction is not supported. Use up, down or off.",
      );
    mode = requested || (entry.direction === "too_much" ? "down" : "up");
    const response = await fetch(
      `${base}/data/${encodeURIComponent(entry.file)}`,
    );
    if (!response.ok)
      throw new Error(
        `The ${symbol} record could not be loaded (HTTP ${response.status}).`,
      );
    record = await response.json();
    if (
      !record.gene ||
      record.gene.symbol !== entry.symbol ||
      !record.locus ||
      !Array.isArray(record.staples) ||
      !Array.isArray(record.ranking)
    )
      throw new Error(
        "The returned record is incomplete or does not match this gene.",
      );
    get<HTMLAnchorElement>("json-link").href =
      `${base}/data/${encodeURIComponent(entry.file)}`;
    renderRecord();
  } catch (error) {
    get("loading-state").hidden = true;
    get("record-view").hidden = true;
    get("error-state").hidden = false;
    get("error-message").textContent =
      error instanceof Error
        ? error.message
        : "The record could not be loaded. Try opening it from the data API.";
  }
}
get<HTMLSelectElement>("rank-select").addEventListener("change", () => {
  renderElements();
  drawMap();
});
get<HTMLSelectElement>("intent-select").addEventListener("change", () => {
  const value = get<HTMLSelectElement>("intent-select").value;
  if (!isMode(value)) return;
  mode = value;
  const url = new URL(location.href);
  url.searchParams.set("gene", record.gene.symbol);
  url.searchParams.set("intent", mode);
  history.replaceState(null, "", url);
  updateMode();
});
let lastMapWidth = 0;
new ResizeObserver(([entry]) => {
  const width = Math.round(entry.contentRect.width);
  if (record && width > 0 && width !== lastMapWidth) {
    lastMapWidth = width;
    drawMap();
  }
}).observe(document.querySelector(".map-scroll")!);
void loadRecord();

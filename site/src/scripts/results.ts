import { viewerDialog } from "../lib/viewer-dialog";
import { rankingColors } from "../lib/element-colors";
import { liveGeneRecord } from "../lib/live-genes";
import { k562Record, K562_VIEWER_IDS, type K562Dataset } from "../lib/k562";
import {
  elementLabel,
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
const manifest: { symbol: string; file: string; direction: string | null }[] =
  JSON.parse(app.dataset.manifest || "[]");
const base = app.dataset.base || "";
const disposeDialog = viewerDialog();
let record: GeneRecord;
let mode: Mode;
let activeElement: string | null = null;
let genomeViewer: import("../lib/genome-viewer").GenomeViewer | null = null;
let viewerLoading = false;
let viewerFailed = false;

let browserViewer: import("../lib/genome-browser").GenomeBrowser | null = null;
let browserLoading = false;
let currentView: "3d" | "contacts" | "browser" = "3d";
const selected = () =>
  record.locus.cres.find((c) => c.id === activeElement) ?? null;
const colors = () =>
  rankingColors(record, get<HTMLSelectElement>("rank-select").value);
const selectedColor = () =>
  selected() ? colors().color(selected()!) : "#ad5641";
function syncBrowserViewer() {
  if (browserViewer) {
    browserViewer.update(record, selected(), selectedColor());
    return;
  }
  if (currentView !== "browser" || browserLoading) return;
  browserLoading = true;
  get("igv-loading").hidden = false;
  void import("../lib/genome-browser")
    .then(({ GenomeBrowser }) => {
      browserViewer = new GenomeBrowser();
      browserViewer.update(record, selected(), selectedColor());
      browserViewer.setActive(currentView === "browser");
    })
    .catch(() => {
      get("igv-error").hidden = false;
      get("igv-retry").hidden = true;
      get("igv-error-message").textContent =
        "The genome browser could not be loaded. Reload to retry.";
    })
    .finally(() => {
      browserLoading = false;
      if (!browserViewer) get("igv-loading").hidden = true;
    });
}
function setView(view: typeof currentView) {
  currentView = view;
  for (const [key, panel] of [
    ["3d", "dna-panel"],
    ["contacts", "contact-panel"],
    ["browser", "browser-panel"],
  ] as const) {
    get(panel).hidden = key !== view;
    get(`view-${key}`).setAttribute("aria-pressed", String(key === view));
  }
  for (const id of ["dna-reset", "dna-zoom-in", "dna-zoom-out"])
    get(id).hidden = view !== "3d";
  genomeViewer?.setActive(view === "3d");
  browserViewer?.setActive(view === "browser");
  if (record) syncBrowserViewer();
}
for (const view of ["3d", "contacts", "browser"] as const)
  get(`view-${view}`).addEventListener("click", () => setView(view));
get("dna-fallback-map").addEventListener("click", () => setView("contacts"));
function syncGenomeViewer() {
  get("genome-viewer").style.setProperty(
    "--selected-element-color",
    selectedColor(),
  );
  syncBrowserViewer();
  if (genomeViewer) {
    genomeViewer.update(record, selected(), selectedColor());
    return;
  }
  if (viewerLoading || viewerFailed || !selected()) return;
  viewerLoading = true;
  void import("../lib/genome-viewer")
    .then(({ GenomeViewer }) => {
      genomeViewer = new GenomeViewer();
      genomeViewer.update(record, selected(), selectedColor());
      genomeViewer.setActive(currentView === "3d");
    })
    .catch(() => {
      viewerFailed = true;
      get("dna-loading").hidden = true;
      get("dna-fallback").hidden = false;
      get("dna-fallback-message").textContent =
        "3D rendering is unavailable. The contact map and genome browser remain available.";
      for (const id of ["dna-reset", "dna-zoom-in", "dna-zoom-out"])
        get<HTMLButtonElement>(id).disabled = true;
    })
    .finally(() => {
      viewerLoading = false;
    });
}
const dispose = () => {
  disposeDialog();
  genomeViewer?.dispose();
  browserViewer?.dispose();
};
window.addEventListener("pagehide", (event) => {
  if (!event.persisted) dispose();
});
if (import.meta.hot) import.meta.hot.dispose(dispose);

function evidenceMarkup() {
  return (
    record.gene.evidence
      .filter(
        (source) =>
          !record.caseStudy || !source.startsWith("Reference transcript"),
      )
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
    criterion === "distance"
      ? "Nearest reference TSS first. Colour uses the relative distance range: nearer is darker."
      : criterion === "re2g"
        ? "Highest supplied enhancer–gene score first. Colour uses the relative score range, not a probability of engineering success."
        : "Highest supplied promoter-contact strength first. Missing values remain unscored.";
  const list = sortElements(record, criterion),
    scale = colors();
  activeElement = selectedElementFor(record, criterion)?.id ?? null;
  const metricLabel =
    criterion === "distance"
      ? "TSS distance"
      : criterion === "re2g"
        ? "rE2G score"
        : "Contact strength";
  const display = (v: number | null) =>
    v === null
      ? "Not supplied"
      : criterion === "distance"
        ? `${fmt(v)} bp`
        : Number(v.toFixed(3)).toString();
  get("score-colour-key").hidden = scale.min === null;
  get("score-colour-key").innerHTML =
    `<span>${display(criterion === "distance" ? scale.max : scale.min)}</span><i aria-hidden="true"></i><span>${display(criterion === "distance" ? scale.min : scale.max)}</span><small>${esc(metricLabel)} · relative to these hits</small>`;
  get("element-list").innerHTML =
    list
      .map((cre, i) => {
        const strength = contactStrength(record, cre),
          distance = Math.abs((cre.start + cre.end) / 2 - record.gene.tss),
          name = elementLabel(record, cre);
        const metric =
          criterion === "distance"
            ? `${fmt(distance)} bp`
            : criterion === "re2g"
              ? cre.score_rE2G.toFixed(3)
              : strength === null
                ? "Not recorded"
                : String(strength);
        return `<article class="element" data-element="${esc(cre.id)}" data-selected="${activeElement === cre.id}" style="--hit-color:${scale.color(cre)}"><div class="element-row"><button type="button" class="element-select" aria-pressed="${activeElement === cre.id}" aria-label="View ${name}"><span class="element-index mono"><i class="hit-swatch" aria-hidden="true"></i>${String(i + 1).padStart(2, "0")}</span><span class="element-name"><strong>${name}</strong><span class="small muted">${fmt(cre.start)}–${fmt(cre.end)}</span></span><span class="element-metric"><span class="small muted">${metricLabel}</span><span class="mono">${metric}</span></span></button><button type="button" class="element-expand" aria-expanded="false" aria-controls="element-detail-${i}" aria-label="Details for ${name}">+</button></div><div class="element-detail" id="element-detail-${i}" hidden><div><p class="eyebrow">Recorded criteria ${record.illustrative ? "· placeholder" : ""}</p><dl class="criteria"><div><dt>rE2G score</dt><dd>${esc(cre.score_rE2G)} <span class="small muted">raw score</span></dd></div><div><dt>Promoter contact</dt><dd>${strength === null ? "Not recorded" : esc(strength)}</dd></div><div><dt>Distance to reference TSS</dt><dd>${fmt(distance)} bp · midpoint</dd></div><div><dt>Potency / achievability</dt><dd>Not supplied</dd></div></dl></div><div><p class="eyebrow">Gene-level literature</p><ul class="evidence-list">${evidenceMarkup()}</ul></div></div></article>`;
      })
      .join("") ||
    '<p class="empty-inline">No regulatory predictions have been supplied for this gene yet. Explore its reference annotation in the genome browser.</p>';
  get("element-list")
    .querySelectorAll<HTMLElement>("[data-element]")
    .forEach((row) => {
      row.querySelector(".element-select")!.addEventListener("click", () => {
        activeElement = row.dataset.element!;
        get("element-list")
          .querySelectorAll<HTMLElement>("[data-element]")
          .forEach((other) => {
            const active = other === row;
            other.dataset.selected = String(active);
            other
              .querySelector(".element-select")!
              .setAttribute("aria-pressed", String(active));
          });
        drawMap();
      });
      row
        .querySelector(".element-expand")!
        .addEventListener("click", (event) => {
          const button = event.currentTarget as HTMLButtonElement,
            details = row.querySelector<HTMLElement>(".element-detail")!;
          const open = details.hidden;
          details.hidden = !open;
          button.setAttribute("aria-expanded", String(open));
          button.textContent = open ? "−" : "+";
        });
    });
}

function drawMap() {
  syncGenomeViewer();
  const { tad, cres, contacts } = record.locus;
  const svg = document.getElementById("contact-map")!;
  const selected = cres.find((c) => c.id === activeElement);
  const selectedContacts = selected ? contactsFor(record, selected) : [];
  get("viewer-selection-kind").textContent = selected
    ? "Selected element"
    : "Target gene";
  get("viewer-element").textContent = selected
    ? elementLabel(record, selected)
    : record.gene.symbol;
  get("viewer-context").textContent = selected
    ? regionLabel({
        chrom: tad.chrom,
        start: selected.start,
        end: selected.end,
      })
    : record.reference
      ? regionLabel(record.reference)
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
  let markup = `<title>${esc(record.gene.symbol)} contact map${selected ? ` — ${esc(elementLabel(record, selected))}` : ""}${record.illustrative ? " — placeholder coordinates" : ""}</title><desc>Linear positions across ${record.caseStudy || record.reference ? "the displayed genomic span (not a TAD)" : "the supplied TAD"}. Square: gene promoter. Circles: regulatory elements. Arcs: supplied contacts, with strength in each arc title. Arc emphasis identifies the selected element, not likelihood.</desc>`;
  markup += `<line x1="32" y1="209" x2="${plotRight}" y2="209" stroke="currentColor" opacity=".28"/>`;
  const intervals = width < 480 ? 2 : 4;
  for (let i = 0; i <= intervals; i++) {
    const p = start + ((end - start) * i) / intervals;
    markup += `<line x1="${x(p)}" y1="203" x2="${x(p)}" y2="215" stroke="currentColor" opacity=".25"/><text x="${x(p)}" y="286" text-anchor="${i === 0 ? "start" : i === intervals ? "end" : "middle"}" class="axis-label">${fmt(Math.round(p))}</text>`;
  }
  if (record.caseStudy)
    plotCREs.forEach((cre) => {
      const a = x(record.gene.tss),
        b = x((cre.start + cre.end) / 2);
      const peak = 209 - Math.min(190, Math.max(45, Math.abs(b - a) * 0.48));
      markup += `<path d="M${a} 209 C${a} ${peak} ${b} ${peak} ${b} 209" fill="none" stroke="${colors().color(cre)}" stroke-dasharray="4 4" stroke-width="${cre.id === activeElement ? 3 : 1}" opacity="${cre.id === activeElement ? 0.95 : 0.13}"><title>${esc(elementLabel(record, cre))} · predicted regulatory link, not a measured contact</title></path>`;
    });
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
    markup += `<g><title>${esc(elementLabel(record, cre))} · ${esc(cre.type)} · ${fmt(cre.start)}–${fmt(cre.end)}</title><circle cx="${pos}" cy="209" r="${activeElement === cre.id ? 7 : 4}" fill="${colors().color(cre)}" stroke="#756a50" stroke-width="0.6" opacity="${activeElement === cre.id ? 1 : 0.8}"/>${activeElement === cre.id ? `<text x="${pos}" y="256" text-anchor="${pos > width - 75 ? "end" : pos < 75 ? "start" : "middle"}" class="cre-label" font-weight="600">${esc(elementLabel(record, cre))}</text>` : ""}</g>`;
  });
  if (promoterInside)
    markup += `<rect x="${x(record.gene.tss) - 7}" y="202" width="14" height="14" fill="#2A6FA8"/><text x="${x(record.gene.tss)}" y="239" text-anchor="middle" class="promoter-label">${esc(record.gene.symbol)} promoter</text>`;
  if (!shown.length && !record.caseStudy)
    markup += `<text x="${width / 2}" y="100" text-anchor="middle" class="axis-label">No in-domain contacts supplied</text>`;
  svg.innerHTML = markup;
  const omitted = contacts.length - shown.length;
  get("map-note").textContent = record.reference
    ? "GRCh38 reference coordinates. No regulatory links supplied; this span is not an identified TAD."
    : record.caseStudy
      ? "GRCh38, BED coordinates. Dashed arcs show predicted enhancer–gene links, not measured contacts. The displayed span is not an identified TAD."
      : `${record.illustrative ? "Illustrative coordinates. " : ""}Schematic contacts, not a predicted 3D structure. ${omitted ? `${omitted} out-of-domain contact(s) omitted. ` : ""}${!promoterInside ? "The supplied promoter is outside this TAD. " : ""}${selected && !selectedContacts.length ? "No promoter contact is recorded for this element." : ""}`;
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
  get("disease-line").textContent = record.reference
    ? record.reference.name
    : record.caseStudy
      ? ""
      : record.gene.disease;
  get("record-notice").innerHTML = record.reference
    ? '<span class="notice-icon">◇</span><div><strong>Live reference annotation.</strong><p>Regulatory predictions have not been supplied for this gene.</p></div>'
    : record.illustrative
      ? '<span class="notice-icon">◇</span><div><strong>Illustrative sample — not a biological finding.</strong></div>'
      : '<span class="notice-icon">◇</span><div><strong>Sample predictions · experimentally unvalidated.</strong></div>';
  get("map-label").textContent = record.illustrative
    ? "Placeholder schematic"
    : record.caseStudy
      ? "Predicted links"
      : "Reference annotation";
  const rank = get<HTMLSelectElement>("rank-select");
  rank.querySelector<HTMLOptionElement>('[value="contact"]')!.disabled =
    !record.locus.contacts.length;
  if (record.caseStudy || record.reference) rank.value = "re2g";
  rank.disabled = !record.locus.cres.length;
  get("elements-count").textContent =
    `${record.locus.cres.length} regulatory elements`;
  get("loading-state").hidden = true;
  get("record-view").hidden = false;
  renderElements();
  drawMap();
  updateMode();
  get<HTMLButtonElement>("view-3d").disabled = !record.locus.cres.length;
  if (record.reference) setView("browser");
}

async function loadRecord() {
  const params = new URLSearchParams(location.search);
  const symbol = (params.get("gene") || manifest[0]?.symbol || "")
    .trim()
    .toUpperCase();
  const entry = manifest.find((item) => item.symbol === symbol);
  try {
    const requested = params.get("intent");
    if (requested !== null && !isMode(requested))
      throw new Error("The requested regulation direction is not supported.");
    mode = requested || (entry?.direction === "too_much" ? "down" : "up");
    if (params.get("source") === "reference" || !entry) {
      record = await liveGeneRecord(symbol);
      get("json-link").hidden = true;
    } else {
      const file =
        symbol === "MYC"
          ? "case-studies/myc-k562.json"
          : encodeURIComponent(entry.file);
      const response = await fetch(`${base}/data/${file}`);
      if (!response.ok)
        throw new Error(
          `The sample could not be loaded (HTTP ${response.status}).`,
        );
      const data = await response.json();
      record =
        symbol === "MYC"
          ? k562Record(data as K562Dataset, K562_VIEWER_IDS)
          : data;
      get<HTMLAnchorElement>("json-link").href = `${base}/data/${file}`;
    }
    if (
      !record.gene ||
      (!record.reference && record.gene.symbol !== entry?.symbol) ||
      !record.locus ||
      !Array.isArray(record.staples) ||
      !Array.isArray(record.ranking)
    )
      throw new Error(
        "The returned record is incomplete or does not match this gene.",
      );
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

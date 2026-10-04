import {
  escapeHTML as esc,
  regionLabel,
  type GeneRecord,
  type Mode,
} from "../lib/records";
import {
  assertPipelineRecord,
  rankedEnhancers,
  designElement,
  designColor,
  riskLabel,
  responseLabel,
  confidenceReasons,
  percent,
  type PipelineAnalysis,
  type PipelineDesign,
} from "../lib/pipeline-results";
import { strategies } from "../data/strategies";

const el = <T extends HTMLElement = HTMLElement>(id: string) =>
  document.getElementById(id) as T;
const format = (n: number) => n.toLocaleString("en-GB");
const sourceLink = (url: string, label: string) =>
  `<a href="${esc(url)}" target="_blank" rel="noopener noreferrer">${esc(label)} ↗</a>`;
// A small literature-based set of direct looping approaches, not a per-locus ranking.
const methodOptions = strategies.filter((s) =>
  ["CLOuD9", "LADL", "BPCL"].includes(s.name),
);

export class PipelineResults {
  private record: GeneRecord;
  private analysis: PipelineAnalysis;
  private selected: PipelineDesign | null = null;
  private mode: Mode;
  private showDesign = false;
  private view: "3d" | "contacts" | "browser" = "3d";
  private dna: import("../lib/genome-viewer").GenomeViewer | null = null;
  private browser: import("../lib/genome-browser").GenomeBrowser | null = null;
  private dnaLoading = false;
  private browserLoading = false;
  private disposed = false;
  private abort = new AbortController();

  constructor(record: GeneRecord, mode: Mode, availableModes: Mode[] = record.analysis!.supported_modes) {
    assertPipelineRecord(record);
    this.record = record;
    this.analysis = record.analysis!;
    this.mode = mode;
    el("results-app").classList.add("pipeline-results");
    el("gene-title").textContent = record.gene.symbol;
    document.title = `${record.gene.symbol} · K562 — Rewire Bio`;
    el("disease-line").textContent = record.gene.disease;
    el("record-notice").innerHTML =
      '<span class="status-dot"></span><strong>K562 · GRCh38 · 5 kb</strong><span class="analysis-status">In silico · unvalidated</span>';
    document.querySelector("#elements h2")!.textContent = "Regulatory elements";
    el("elements-count").textContent =
      `${this.rows().length} candidate region${this.rows().length === 1 ? "" : "s"}`;
    el("rank-select").closest("label")!.hidden = true;
    el("score-colour-key").hidden = true;
    el("view-contacts").textContent = "Genomic map";
    el("dna-fallback-map").textContent = "Open genomic map →";
    el("viewer-selection-kind").hidden = true;
    el("viewer-context").hidden = true;
    el("design-nav").textContent = "Stapling strategy";
    el("staples").hidden = true;
    const panel = document.createElement("section");
    panel.id = "selected-design";
    panel.className = "associated-design";
    panel.setAttribute("aria-labelledby", "design-heading");
    el("element-list").after(panel);
    const dropped = this.analysis.source.dropped;
    document.querySelector("#elements .methodology")!.innerHTML =
      `<summary>About these results <span aria-hidden="true">+</span></summary><p>Each enhancer is associated with one exported anchor pair. Regions are ordered by that pair’s reference-setting ΔT/T: simulated regulatory-input change, not expression change or an intrinsic enhancer score.</p><p>The shortlist excludes ${dropped.blocked_by_risk || 0} risk-gated and ${dropped.below_min_effect || 0} below-threshold designs. Retained alternatives and the original scores are available in View JSON.</p><p>This is a K562 demonstration, not validation in disease-relevant tissue. The 2 Mb interval is an analysis window, not a called TAD. Only ClinGen dosage sensitivity was screened; COSMIC and DepMap were not assessed.</p><p>${sourceLink("https://github.com/nathnwe/hack-anthropic-modal-2026/blob/main/docs/method.md", "Method and assumptions")}</p>`;
    el("design-visibility").hidden = false;
    el("show-design").addEventListener(
      "click",
      () => {
        this.showDesign = !this.showDesign;
        el("show-design").setAttribute("aria-checked", String(this.showDesign));
        el("show-design-state").textContent = this.showDesign ? "Yes" : "No";
        el("genome-viewer").dataset.designVisible = String(this.showDesign);
        if (this.selected) {
          this.drawMap(this.selected);
          void this.syncViewers();
        }
      },
      { signal: this.abort.signal },
    );
    el("genome-viewer").dataset.designVisible = "false";
    el("intent-select").innerHTML =
      ([['up', 'Upregulate'], ['down', 'Downregulate']] as const).map(([value,label]) =>
        `<option value="${value}">${label}${availableModes.includes(value) ? '' : ' — unavailable'}</option>`).join('');
    el("intent-select").addEventListener(
      "change",
      () => {
        this.mode = el<HTMLSelectElement>("intent-select").value as Mode;
        const url = new URL(location.href);
        url.searchParams.set("intent", this.mode);
        if (availableModes.includes(this.mode) && !this.analysis.supported_modes.includes(this.mode as "up" | "down")) {
          url.searchParams.set("source", "analysis");
          url.searchParams.delete("design");
          location.assign(url);
          return;
        }
        history.replaceState(null, "", url);
        this.updateMode();
      },
      { signal: this.abort.signal },
    );
    el("igv-retry").addEventListener(
      "click",
      () => {
        if (!this.browser) void this.syncViewers();
      },
      { signal: this.abort.signal },
    );
    el("loading-state").hidden = true;
    el("record-view").hidden = false;
    this.updateMode();
  }
  private rows() {
    return rankedEnhancers(this.analysis);
  }
  private name(d: PipelineDesign) {
    return `Enhancer ${String(this.rows().findIndex((row) => row.id === d.id) + 1).padStart(2, "0")}`;
  }
  private color(d: PipelineDesign) {
    return designColor(this.rows(), d, "pipeline");
  }
  private updateMode() {
    el("intent-title").textContent =
      this.mode === "up" ? "↑ Upregulate" : "↓ Downregulate";
    el<HTMLSelectElement>("intent-select").value = this.mode;
    const unsupported = !this.analysis.supported_modes.includes(this.mode as "up" | "down");
    el("landscape").hidden =
      el("selected-design").hidden =
      el("design-nav").hidden =
        unsupported;
    if (unsupported) {
      this.dna?.setActive(false);
      this.browser?.setActive(false);
      el("element-list").innerHTML =
        `<div class="empty-inline"><h3>${this.mode === "down" ? "Downregulation" : "Upregulation"} is not modelled in this analysis.</h3><p>These K562 results evaluate ${this.analysis.supported_modes[0] === "down" ? "decreased" : "increased"} regulatory input.</p><button class="text-link" id="switch-up" type="button">View ${this.analysis.supported_modes[0] === "down" ? "downregulation" : "upregulation"} results →</button></div>`;
      el("switch-up").addEventListener(
        "click",
        () => {
          el<HTMLSelectElement>("intent-select").value = this.analysis.supported_modes[0];
          el("intent-select").dispatchEvent(new Event("change"));
        },
        { signal: this.abort.signal },
      );
      return;
    }
    this.renderList();
    this.dna?.setActive(this.view === "3d");
    this.browser?.setActive(this.view === "browser");
  }
  private renderList() {
    const rows = this.rows();
    el("element-list").innerHTML =
      rows
        .map(
          (
            d,
            i,
          ) => `<article class="element pipeline-enhancer" data-design="${esc(d.id)}" style="--hit-color:${this.color(d)}">
      <div class="element-row"><button type="button" class="element-select" data-select-design="${esc(d.id)}" aria-pressed="false" aria-label="View ${this.name(d)}">
      <span class="element-index mono"><i class="hit-swatch" aria-hidden="true"></i>${String(i + 1).padStart(2, "0")}</span>
      <span class="element-name"><strong>${this.name(d)}</strong><span class="small muted">${esc(regionLabel(d.element))}</span></span>
      </button>
      <button type="button" class="element-expand" data-expand-design="${esc(d.id)}" aria-label="Details for ${this.name(d)}" aria-expanded="false" aria-controls="detail-${esc(d.id)}">+</button></div>
      <div id="detail-${esc(d.id)}" class="element-detail enhancer-detail" hidden><div class="reference-score"><span>Reference-setting ΔT/T</span><strong>${percent(d.delta_target)}</strong></div><p>Modelled input, not expression. p99: ${d.delta_at_p99 === null ? "not supplied" : percent(d.delta_at_p99, 2)}.</p><p class="response-conclusion">${esc(responseLabel(d))}.</p></div></article>`,
        )
        .join("") ||
      '<p class="empty-inline">No regulatory region has a design that passes the current analysis thresholds.</p>';
    el("element-list")
      .querySelectorAll<HTMLButtonElement>("[data-select-design]")
      .forEach((button) =>
        button.addEventListener(
          "click",
          () => this.select(button.dataset.selectDesign!),
          { signal: this.abort.signal },
        ),
      );
    el("element-list")
      .querySelectorAll<HTMLButtonElement>("[data-expand-design]")
      .forEach((button) =>
        button.addEventListener(
          "click",
          () => {
            const detail = el(button.getAttribute("aria-controls")!);
            detail.hidden = !detail.hidden;
            button.setAttribute("aria-expanded", String(!detail.hidden));
            button.textContent = detail.hidden ? "+" : "−";
          },
          { signal: this.abort.signal },
        ),
      );
    const requested =
      this.selected?.id || new URLSearchParams(location.search).get("design");
    const first = rows.find((d) => d.id === requested) || rows[0];
    if (first) this.select(first.id);
    else {
      this.selected = null;
      el("selected-design").hidden =
        el("landscape").hidden =
        el("design-nav").hidden =
          true;
    }
  }
  private select(id: string) {
    const d = this.rows().find((d) => d.id === id);
    if (!d) return;
    this.selected = d;
    document
      .querySelectorAll<HTMLElement>(".pipeline-enhancer[data-design]")
      .forEach((row) => {
        const active = row.dataset.design === id;
        row.dataset.selected = String(active);
        row
          .querySelector("[data-select-design]")!
          .setAttribute("aria-pressed", String(active));
      });
    el("results-app").dataset.design = id;
    const url = new URL(location.href);
    url.searchParams.set("design", id);
    history.replaceState(null, "", url);
    el("viewer-selection-kind").textContent = "Selected regulatory region";
    el("viewer-element").textContent =
      `${this.name(d)} → ${this.record.gene.symbol}`;
    el("map-label").textContent = "Schematic";
    el("viewer-context").textContent = `${regionLabel(d.element)} · GRCh38`;
    el("genome-viewer").style.setProperty(
      "--selected-element-color",
      this.color(d),
    );
    this.renderDesign(d);
    this.drawMap(d);
    void this.syncViewers();
    el("results-status").textContent =
      `${this.name(d)} selected. Associated anchor pair and viewer updated. ${responseLabel(d)}.`;
  }
  private renderDesign(d: PipelineDesign) {
    const confidence =
      d.confidence === "CONFIDENT" ? "Stable direction" : "Tentative";
    el("selected-design").dataset.design = d.id;
    el("selected-design").innerHTML =
      `<div class="section-heading"><h2 id="design-heading">Stapling strategy</h2><span class="tag">${confidence}</span></div>
      <div class="staple-regions"><h3>Staple target regions</h3>
      <dl class="anchor-pair"><div><dt><span>A</span> Near enhancer</dt><dd>${esc(regionLabel(d.anchor_a))}</dd></div><div><dt><span>B</span> ${d.mode === "down" ? "Decoy region" : `Near ${esc(this.record.gene.symbol)}`}</dt><dd>${esc(regionLabel(d.anchor_b))}</dd></div></dl>
      <p class="anchor-distance">${format(d.separation_kb)} kb apart</p>
      <details class="design-checks"><summary>Model checks <span aria-hidden="true">+</span></summary><p>5 kb candidate regions; guide sequences and binding sites have not been designed.${d.mode === "down" ? " Decoy-mediated lowering depends on assumed contact conservation; it is not experimentally validated." : ""}</p><p>${esc(responseLabel(d))}.</p><ul>${confidenceReasons(
        this.analysis,
        d,
      )
        .map((r) => `<li>${esc(r)}</li>`)
        .join(
          "",
        )}</ul><p><strong>${esc(riskLabel(d.risk))}.</strong> ${esc(d.risk.reason)}.</p>${d.unflagged_collateral_over_tau.length ? `<p>Other modelled changes: ${d.unflagged_collateral_over_tau.map((x) => `${esc(x.gene)} (${percent(x.delta)})`).join(", ")}. These genes are not flagged by the available ClinGen list.</p>` : ""}<p>ClinGen-only screening; COSMIC and DepMap were not assessed. An empty screen does not establish safety.</p>${d.sweep ? `<table class="sweep-table"><caption>Regulatory-input change by contact-strength setting</caption><thead><tr><th>Setting</th><th>ΔT/T</th></tr></thead><tbody>${d.sweep.map((v, i) => `<tr><th>p${this.analysis.sweep_percentiles[i]}</th><td>${percent(v, 2)}</td></tr>`).join("")}<tr><th>Reference</th><td>${percent(d.delta_target)}</td></tr></tbody></table><p>The reference uses the strongest natural contact at this separation, above p99. These settings are not experimental doses or probabilities.</p>` : "<p>No strength sweep was retained.</p>"}</details></div>
      <div class="suggested-strategies"><h3>Suggested approaches</h3>${methodOptions.map((s) => `<details class="method-option"><summary><strong>${esc(s.name)}</strong><span>${esc(s.mechanism)}</span><b aria-hidden="true">+</b></summary><div><p>${esc(s.attraction)}</p><p>Suitability for this locus has not been assessed.</p><p><strong>Evidence:</strong> ${esc(s.evidence)}.</p><p><strong>Limitation:</strong> ${esc(s.issue)}</p>${sourceLink(s.source, s.citation)}</div></details>`).join("")}</div>`;
  }
  setView(view: "3d" | "contacts" | "browser") {
    this.view = view;
    el("map-label").hidden = view === "browser";
    for (const [key, panel] of [
      ["3d", "dna-panel"],
      ["contacts", "contact-panel"],
      ["browser", "browser-panel"],
    ]) {
      el(panel).hidden = key !== view;
      el(`view-${key}`).setAttribute("aria-pressed", String(key === view));
    }
    for (const id of ["dna-reset", "dna-zoom-in", "dna-zoom-out"])
      el(id).hidden = view !== "3d";
    this.dna?.setActive(view === "3d");
    this.browser?.setActive(view === "browser");
    void this.syncViewers();
  }
  private async syncViewers() {
    if (this.disposed || !this.selected || !this.analysis.supported_modes.includes(this.mode as "up" | "down")) return;
    const d = this.selected,
      cre = designElement(d),
      color = this.color(d);
    if (this.dna)
      this.dna.update(this.record, cre, color, this.showDesign ? d : null);
    else if (this.view === "3d" && !this.dnaLoading) {
      this.dnaLoading = true;
      try {
        const { GenomeViewer } = await import("../lib/genome-viewer");
        if (this.disposed) return;
        this.dna = new GenomeViewer();
        if (this.selected)
          this.dna.update(
            this.record,
            designElement(this.selected),
            this.color(this.selected),
            this.showDesign ? this.selected : null,
          );
        this.dna.setActive(this.view === "3d");
      } catch (error) {
        el("dna-loading").hidden = true;
        el("dna-fallback").hidden = false;
        el("dna-fallback-message").textContent =
          "3D is unavailable. The genomic map and interval browsers remain available.";
        console.error("DNA viewer:", error);
      } finally {
        this.dnaLoading = false;
      }
    }
    if (this.browser)
      this.browser.update(
        this.record,
        cre,
        color,
        this.showDesign ? d : undefined,
      );
    else if (this.view === "browser" && !this.browserLoading) {
      this.browserLoading = true;
      try {
        const { GenomeBrowser } = await import("../lib/genome-browser");
        if (this.disposed) return;
        this.browser = new GenomeBrowser();
        if (this.selected)
          this.browser.update(
            this.record,
            designElement(this.selected),
            this.color(this.selected),
            this.showDesign ? this.selected : undefined,
          );
        this.browser.setActive(this.view === "browser");
      } catch (error) {
        el("igv-error").hidden = false;
        el("igv-error-message").textContent =
          "The genome browser could not be loaded. Retry to reconnect.";
        console.error(error);
      } finally {
        this.browserLoading = false;
      }
    }
  }
  private drawMap(d: PipelineDesign) {
    const { start, end } = this.record.locus.tad;
    const x = (p: number) => 45 + ((p - start) / (end - start)) * 590;
    const tss = this.record.gene.tss,
      enhancer = (d.element.start + d.element.end) / 2;
    const a = (d.anchor_a.start + d.anchor_a.end) / 2,
      b = (d.anchor_b.start + d.anchor_b.end) / 2;
    const svg = el("contact-map");
    svg.setAttribute("viewBox", "0 0 680 285");
    const from = this.showDesign ? a : enhancer,
      to = this.showDesign ? b : tss;
    svg.innerHTML = `<title>${esc(this.name(d))} and ${esc(this.record.gene.symbol)}${this.showDesign ? " with candidate anchor regions" : ""}</title><desc>Genomic coordinates within the 2 Mb analysis window. The dashed arc is schematic, not a measured contact.</desc><text x="35" y="28" class="axis-label">2 Mb analysis window · GRCh38</text><path d="M${x(from)} 175 C${x(from)} 48 ${x(to)} 48 ${x(to)} 175" fill="none" stroke="${this.showDesign ? "#52683f" : this.color(d)}" stroke-width="2" stroke-dasharray="5 5"/><line x1="45" y1="175" x2="635" y2="175" stroke="#a9a899"/>${[
      0, 1, 2, 3, 4,
    ]
      .map((i) => {
        const p = start + ((end - start) * i) / 4;
        return `<line x1="${x(p)}" y1="175" x2="${x(p)}" y2="181" stroke="#858879"/><text x="${x(p)}" y="263" text-anchor="${i === 0 ? "start" : i === 4 ? "end" : "middle"}" class="axis-label">${(p / 1e6).toFixed(2)} Mb</text>`;
      })
      .join(
        "",
      )}<rect x="${x(tss) - 4}" y="168" width="8" height="14" fill="#2a6fa8"/><circle cx="${x(enhancer)}" cy="175" r="5" fill="${this.color(d)}"/><text x="${x(tss)}" y="223" text-anchor="middle" class="axis-label">${esc(this.record.gene.symbol)} TSS</text><text x="${x(enhancer)}" y="238" text-anchor="middle" class="axis-label">${esc(this.name(d))}</text>${
      this.showDesign
        ? [
            [d.anchor_a, "A"],
            [d.anchor_b, "B"],
          ]
            .map(([region, label]) => {
              const r = region as PipelineDesign["anchor_a"];
              return `<g data-anchor="${label}"><rect x="${x(r.start)}" y="163" width="${Math.max(3, x(r.end) - x(r.start))}" height="24" fill="#52683f44" stroke="#52683f"/><text x="${x((r.start + r.end) / 2)}" y="152" text-anchor="middle" class="axis-label">${label}</text></g>`;
            })
            .join("")
        : ""
    }`;
    svg.dataset.designVisible = String(this.showDesign);
    el("map-note").textContent = this.showDesign
      ? "A/B: candidate anchor regions. Dashed arc: proposed tether, not measured contact."
      : "Dashed arc: schematic enhancer–target relationship. No measured contact matrix is displayed.";
    document.querySelector(".map-legend")!.innerHTML =
      `<span><i class="promoter-key"></i>Target TSS</span><span><i class="element-key"></i>Enhancer region</span>${this.showDesign ? "<span>A / B · anchor regions</span>" : ""}`;
  }
  dispose() {
    this.disposed = true;
    this.abort.abort();
    this.dna?.dispose();
    this.browser?.dispose();
  }
}

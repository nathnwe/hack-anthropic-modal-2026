import igv, { type Browser, type CreateOpt, type TrackLoad } from "igv";
import type { PipelineDesign } from "./pipeline-results";
import type { Region } from "./records";
import { type GeneRecord, type CRE, type ReferenceGene } from "./records";
import { genomeReference, jsonRequest, ENSEMBL } from "./live-genes";

import { annotationFeatures } from "./annotations";
import { contextualTracks, signalTrackConfig } from "./genome-tracks";
import { trackNameGutter } from "./igv-labels";
import { escapeHTML as esc } from "./records";
import { referenceGenes, type ReferenceSnapshot } from "./reference-genes";
let referenceSnapshot: Promise<ReferenceSnapshot> | undefined;

const annotationCache = new Map<
  string,
  Promise<ReturnType<typeof annotationFeatures>>
>();
async function annotations(chr: string, start: number, end: number) {
  const chrom = chr.replace(/^chr/, "").replace(/^M$/, "MT");
  if (!/^(?:[1-9]|1\d|2[0-2]|X|Y|MT)$/.test(chrom)) return [];
  const lo = Math.max(1, Math.floor(start / 10000) * 10000 + 1),
    hi = Math.ceil(end / 10000) * 10000;
  if (hi - lo > 5000000)
    throw new Error("Zoom in to load transcript annotations.");
  const key = `${chrom}:${lo}-${hi}`;
  if (!annotationCache.has(key)) {
    const request = jsonRequest(
      `${ENSEMBL}/overlap/region/homo_sapiens/${key}?feature=transcript;feature=exon;feature=cds;content-type=application/json`,
    ).then((data) => annotationFeatures(Array.isArray(data) ? data : []));
    annotationCache.set(key, request);
    // Cache failures too; an explicit Retry clears them. This avoids repeated
    // annotation timeouts blocking every change of interval or anchor overlay.
    void request.catch(() => {});
    if (annotationCache.size > 32)
      annotationCache.delete(annotationCache.keys().next().value!);
  }
  return annotationCache.get(key)!;
}
const byId = (id: string) => document.getElementById(id)!;
const locus = (chrom: string, start: number, end: number, pad: number) =>
  `${chrom}:${Math.max(1, Math.floor(start - pad) + 1)}-${Math.ceil(end + pad)}`;
export class GenomeBrowser {
  private browsers: Partial<Record<"gene" | "enhancer", Browser>> = {};
  private desired: {
    record: GeneRecord;
    cre: CRE | null;
    color: string;
    design?: PipelineDesign;
  } | null = null;
  private active = false;
  private busy = false;
  private disposed = false;
  private key = "";
  private referenceKey = "";
  private gutterDisposers: (() => void)[] = [];
  private contextStarted = new WeakSet<Browser>();
  private pendingTracks = new Set<Promise<unknown>>();
  private reportedErrors = new Set<Browser>();
  private targets: Partial<Record<"gene" | "enhancer", string>> = {};
  private observer: ResizeObserver;
  private width = 0;
  private abort = new AbortController();
  constructor() {
    for (const pane of ["gene", "enhancer"] as const)
      for (const action of ["in", "out", "fit"] as const)
        byId(`igv-${pane}-${action}`).addEventListener(
          "click",
          () => {
            const browser = this.browsers[pane];
            if (!browser) return;
            if (action === "fit" && this.targets[pane])
              void browser.search(this.targets[pane]!);
            else if (action === "in") browser.zoomIn();
            else if (action === "out") browser.zoomOut();
          },
          { signal: this.abort.signal },
        );
    byId("igv-retry").addEventListener(
      "click",
      () => {
        if (this.busy) return;
        if (this.pendingTracks.size) return;
        annotationCache.clear();
        this.reportedErrors.clear();
        this.gutterDisposers.splice(0).forEach((dispose) => dispose());
        Object.values(this.browsers).forEach((b) => igv.removeBrowser(b));
        this.browsers = {};
        this.referenceKey = "";
        this.key = "";
        void this.refresh();
      },
      { signal: this.abort.signal },
    );
    this.observer = new ResizeObserver(([entry]) => {
      const w = Math.round(entry.contentRect.width);
      if (this.active && w > 0 && w !== this.width) {
        this.width = w;
        Object.values(this.browsers).forEach((b) => b.visibilityChange());
      }
    });
    this.observer.observe(byId("browser-panel"));
  }
  update(
    record: GeneRecord,
    cre: CRE | null,
    color: string,
    design?: PipelineDesign,
  ) {
    this.desired = { record, cre, color, design };
    if (this.active) void this.refresh();
  }
  setActive(active: boolean) {
    this.active = active;
    if (active) {
      Object.values(this.browsers).forEach((b) => b.visibilityChange());
      void this.refresh();
    }
  }
  private async create(pane: "gene" | "enhancer", target: string) {
    const config = {
      reference: {
        id: "hg38",
        name: "Human GRCh38",
        format: "chromsizes",
        url: `${import.meta.env.BASE_URL.replace(/\/$/, "")}/data/reference/hg38.chrom.sizes`,
      },
      loadDefaultGenomes: false,
      locus: target,
      showNavigation: false,
      showIdeogram: false,
      showRuler: true,
      showSVGButton: false,
      showSequence: false,
      showTrackLabels: false,
      showTrackDragHandles: false,
      showGearColumn: false,
      tracks: [],
    };
    const browser = await igv.createBrowser(
      byId(`igv-${pane}`),
      config as unknown as CreateOpt,
    );
    if (this.disposed) {
      igv.removeBrowser(browser);
      return null;
    }
    this.browsers[pane] = browser;
    // IGV also raises network errors during pan/zoom, after loadTrack has
    // resolved. Keep those errors visible in our inline status instead of
    // letting the library move keyboard focus to its floating alert dialog.
    const alerts = (browser as unknown as {
      alert: { present: (error: unknown, callback?: () => void) => void };
    }).alert;
    const originalAlert = alerts.present.bind(alerts);
    alerts.present = (error, callback) => {
      const message = error instanceof Error ? error.message : String(error);
      if (!/resource|fetch|network|timed out/i.test(message)) {
        originalAlert(error, callback);
        return;
      }
      if (this.disposed || this.browsers[pane] !== browser) return;
      this.reportedErrors.add(browser);
      const track = this.desired && contextualTracks(this.desired.record)
        .find((item) => message.includes(item.accession));
      byId(`igv-${pane}-status`).textContent =
        `${track?.name || "A remote track"} unavailable.`;
      byId("igv-error").hidden = false;
      byId("igv-error-message").textContent =
        "A remote track could not load. Retry to reconnect.";
      byId("igv-retry").hidden = false;
      callback?.();
    };
    this.gutterDisposers.push(
      trackNameGutter(
        byId(`igv-${pane}`),
        byId(`igv-${pane}-labels`),
        (name) => {
          const assay =
            this.desired &&
            contextualTracks(this.desired.record).find((t) => t.name === name);
          if (assay)
            return {
              color: assay.color,
              detail: "Fold / control",
              url: assay.file_url,
            };
          if (name === "Anchor A")
            return { color: "#52683f", detail: "Near enhancer" };
          if (name === "Anchor B")
            return {
              color: "#52683f",
              detail: this.desired?.design?.mode === "down" ? "Decoy region" : `Near ${this.desired?.record.gene.symbol || "gene"}`,
            };
          if (name === "Enhancer") return { color: this.desired?.color };
          if (name === "Target gene" || name === "TSS")
            return { color: "#2a6fa8" };
          return { color: "#626657" };
        },
      ),
    );
    const label = () => {
      byId(`igv-${pane}-locus`).textContent = [browser.currentLoci()]
        .flat()
        .map((value) =>
          value.replace(
            /:([\d.]+)-([\d.]+)/,
            (_, start, end) =>
              `:${Math.floor(Number(start)).toLocaleString("en-US")}–${Math.ceil(Number(end)).toLocaleString("en-US")}`,
          ),
        )
        .join(" · ");
    };
    browser.on("locuschange", label);
    label();
    return browser;
  }
  private startContextTracks(pane: "gene" | "enhancer", browser: Browser) {
    if (this.contextStarted.has(browser) || !this.desired) return;
    this.contextStarted.add(browser);
    const tracks = contextualTracks(this.desired.record);
    const sources = byId("igv-assay-sources");
    sources.hidden = !tracks.length;
    sources.innerHTML = tracks.length
      ? `<p>K562 · GRCh38. <a href="${import.meta.env.BASE_URL.replace(/\/$/, "")}/data/reference/k562-refseq.json" target="_blank">RefSeq gene models</a> are cached for the six demo windows (one representative transcript per gene; not necessarily canonical). Outside those windows, annotations use Ensembl. Experimental context; these tracks do not change the model score. Each track auto-scales independently in each pane.</p><p>Signal: fold enrichment over the processing pipeline’s control/background; pooled biological replicates. Track names link to the original files.</p><ul>${tracks
          .map(
            (t) =>
              `<li><a href="${esc(t.experiment_url)}" target="_blank" rel="noopener noreferrer">${esc(t.name)} · ${esc(t.experiment)}</a> — ${esc(t.lab)}; replicates ${t.biological_replicates.join(", ")}.<br><a href="${esc(t.file_url)}" target="_blank" rel="noopener noreferrer">${esc(t.accession)}</a>. ENCODE experiment flags: ${esc(
                Object.entries(t.experiment_audits)
                  .map(
                    ([level, categories]) =>
                      `${level.toLowerCase()}: ${categories.join(", ")}`,
                  )
                  .join("; "),
              )}.</li>`,
          )
          .join("")}</ul>`
      : "";
    const failed = new Set<string>();
    let remaining = tracks.length + 1;
    const status = byId(`igv-${pane}-status`);
    status.textContent = "Loading tracks…";
    const current = () => !this.disposed && this.browsers[pane] === browser;
    const finish = () => {
      if (!current()) return;
      status.textContent = failed.size
        ? `${[...failed].join(", ")} unavailable.`
        : remaining
          ? "Loading tracks…"
          : "";
      if (failed.size) {
        this.reportedErrors.add(browser);
        byId("igv-error").hidden = false;
        byId("igv-error-message").textContent = "Some tracks are unavailable.";
        byId("igv-retry").hidden = false;
      }
    };
    const genes = {
      id: "gene-annotations",
      name: this.desired.record.analysis ? "RefSeq genes" : "Genes",
      type: "annotation",
      color: "#626657",
      order: 3,
      displayMode: "EXPANDED",
      height: 54,
      maxHeight: 90,
      visibilityWindow: 1000000,
      reader: {
        readHeader: async () => ({}),
        readFeatures: async (chr: string, start: number, end: number) => {
          try {
            let features;
            if (tracks.length && this.desired?.record.analysis) {
              referenceSnapshot ??= jsonRequest(
                `${import.meta.env.BASE_URL.replace(/\/$/, "")}/data/reference/k562-refseq.json`,
              ) as Promise<ReferenceSnapshot>;
              try {
                features = referenceGenes(
                  await referenceSnapshot,
                  chr,
                  start,
                  end,
                );
              } catch {
                referenceSnapshot = undefined;
              }
            }
            // Outside the supplied windows, retain the live Ensembl lookup.
            features ??= await annotations(chr, start, end);
            failed.delete("Genes");
            finish();
            return features;
          } catch {
            failed.add("Genes");
            finish();
            return [];
          }
        },
      },
    };
    for (const config of [genes, ...tracks.map(signalTrackConfig)]) {
      const load = async () => {
        const attempts = config.id.startsWith("encode-") ? 2 : 1;
        for (let attempt = 0; attempt < attempts; attempt++) {
          if (!current()) return;
          try {
            const track = await browser.loadTrack(
              config as unknown as TrackLoad<"annotation">,
            );
            if (!track) throw new Error("Track was not loaded");
            return;
          } catch (error) {
            if (attempt + 1 === attempts) throw error;
            browser
              .findTracks("id", config.id)
              .forEach((track) => browser.removeTrack(track));
          }
        }
      };
      const task = load()
        .catch((error) => {
          failed.add(config.name);
          console.warn(`Genome track ${config.name}:`, error);
        })
        .finally(() => {
          remaining--;
          this.pendingTracks.delete(task);
          finish();
          if (!remaining && current()) browser.visibilityChange();
        });
      this.pendingTracks.add(task);
    }
  }
  private async show(
    pane: "gene" | "enhancer",
    ref: ReferenceGene,
    start: number,
    end: number,
    color: string,
    anchor?: Region,
    pipelineTss?: number,
  ) {
    const target = locus(
      ref.chrom,
      Math.min(start, anchor?.start ?? start),
      Math.max(end, anchor?.end ?? end),
      Math.max(1500, (end - start) * 0.15),
    );
    this.targets[pane] = target;
    let browser = this.browsers[pane];
    if (!browser) browser = (await this.create(pane, target)) || undefined;
    if (!browser || this.disposed) return;
    const existing = browser.findTracks("id", "selected-interval");
    existing.forEach((t) => browser!.removeTrack(t));
    const track: TrackLoad<"annotation"> & {
      id: string;
      margin: number;
      expandedRowHeight: number;
    } = {
      id: "selected-interval",
      name: pane === "gene" ? "Target gene" : "Enhancer",
      type: "annotation",
      features: [
        {
          chr: ref.chrom,
          start,
          end,
          name: "",
          ...(pane === "gene" ? { strand: ref.strand } : {}),
        },
      ],
      color,
      displayMode: "COLLAPSED",
      height: 26,
      minHeight: 26,
      maxHeight: 26,
      margin: 4,
      expandedRowHeight: 22,
      order: 0,
    };
    await browser.loadTrack(track);
    for (const id of ["pipeline-anchor", "pipeline-tss"])
      browser.findTracks("id", id).forEach((t) => browser!.removeTrack(t));
    if (anchor)
      await browser.loadTrack({
        id: "pipeline-anchor",
        name: `Anchor ${pane === "gene" ? "B" : "A"}`,
        type: "annotation",
        features: [
          {
            chr: anchor.chrom,
            start: anchor.start,
            end: anchor.end,
            name: "",
          },
        ],
        color: "#52683f",
        height: 32,
        minHeight: 32,
        maxHeight: 32,
        margin: 8,
        expandedRowHeight: 24,
        displayMode: "COLLAPSED",
        order: 1,
      } as TrackLoad<"annotation">);
    if (pipelineTss !== undefined)
      await browser.loadTrack({
        id: "pipeline-tss",
        name: "TSS",
        type: "annotation",
        features: [
          {
            chr: ref.chrom,
            start: pipelineTss,
            end: pipelineTss + 1,
            name: "",
          },
        ],
        color: "#2a6fa8",
        height: 24,
        minHeight: 24,
        maxHeight: 24,
        margin: 4,
        expandedRowHeight: 20,
        displayMode: "COLLAPSED",
        order: 2,
      } as TrackLoad<"annotation">);
    await Promise.resolve(browser.search(target));
    byId(`igv-${pane}`).dataset.interval = `${ref.chrom}:${start}-${end}`;
    this.startContextTracks(pane, browser);
  }
  private async refresh() {
    if (this.busy || !this.active || !this.desired || this.disposed) return;
    this.busy = true;
    try {
      while (this.active && this.desired && !this.disposed) {
        const { record, cre, color, design } = this.desired;
        const key = JSON.stringify([
          record.gene.symbol,
          cre?.id,
          color,
          design?.id,
        ]);
        if (key === this.key) break;
        const ref = genomeReference(record);
        byId("igv-error").hidden = this.reportedErrors.size === 0;
        byId("igv-content").hidden = !ref;
        if (!ref) {
          byId("igv-error").hidden = false;
          byId("igv-error-message").textContent =
            "A verified genome assembly is not supplied for this example.";
          byId("igv-retry").hidden = true;
          this.key = key;
          break;
        }
        byId("igv-loading").hidden = false;
        byId("igv-content").setAttribute("aria-busy", "true");
        byId("igv-enhancer-pane").hidden = !cre;
        byId("igv-enhancer-empty").hidden = !!cre;
        byId("igv-gene-title").textContent =
          `Target gene · ${record.gene.symbol}`;
        const referenceKey = JSON.stringify([ref, design?.anchor_b]);
        if (referenceKey !== this.referenceKey) {
          await this.show(
            "gene",
            ref,
            ref.start,
            ref.end,
            "#2a6fa8",
            design?.anchor_b,
            record.analysis ? record.gene.tss : undefined,
          );
          this.referenceKey = referenceKey;
        }
        if (cre) {
          byId("igv-enhancer-title").textContent = "Regulatory element";
          await this.show(
            "enhancer",
            ref,
            cre.start,
            cre.end,
            color,
            design?.anchor_a,
          );
        }
        this.key = key;
        byId("browser-panel").dataset.element = cre?.id || "";
        byId("browser-panel").dataset.design = design?.id || "";
        byId("igv-loading").hidden = true;
      }
    } catch (error) {
      byId("igv-content").hidden = true;
      byId("igv-error").hidden = false;
      byId("igv-error-message").textContent =
        error instanceof Error
          ? error.message
          : "The genome browser could not be loaded.";
      byId("igv-retry").hidden = false;
    } finally {
      this.busy = false;
      byId("igv-loading").hidden = true;
      byId("igv-content").removeAttribute("aria-busy");
    }
  }
  dispose() {
    this.disposed = true;
    this.abort.abort();
    this.observer.disconnect();
    this.gutterDisposers.splice(0).forEach((dispose) => dispose());
    Object.values(this.browsers).forEach((b) => igv.removeBrowser(b));
  }
}

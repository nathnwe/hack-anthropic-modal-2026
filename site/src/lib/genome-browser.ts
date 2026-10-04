import igv, { type Browser, type CreateOpt, type TrackLoad } from "igv";
import {
  elementLabel,
  type GeneRecord,
  type CRE,
  type ReferenceGene,
} from "./records";
import { genomeReference, jsonRequest, ENSEMBL } from "./live-genes";

import { annotationFeatures } from "./annotations";

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
    request.catch(() => annotationCache.delete(key));
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
  } | null = null;
  private active = false;
  private busy = false;
  private disposed = false;
  private key = "";
  private referenceKey = "";
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
  update(record: GeneRecord, cre: CRE | null, color: string) {
    this.desired = { record, cre, color };
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
    const reader = {
      readHeader: async () => ({}),
      readFeatures: async (chr: string, start: number, end: number) => {
        try {
          const data = await annotations(chr, start, end);
          byId(`igv-${pane}-status`).textContent = "";
          return data;
        } catch (error) {
          byId(`igv-${pane}-status`).textContent =
            error instanceof Error
              ? error.message
              : "Annotations are unavailable. Try again.";
          byId("igv-error").hidden = false;
          byId("igv-error-message").textContent =
            "Some annotations could not be loaded. Retry to reconnect.";
          byId("igv-retry").hidden = false;
          return [];
        }
      },
    };
    // IGV 3.8 typings predate the documented twoBit reference and display options.
    const config = {
      reference: {
        id: "hg38",
        name: "Human GRCh38",
        twoBitURL:
          "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit",
      },
      loadDefaultGenomes: false,
      locus: target,
      showNavigation: false,
      showIdeogram: false,
      showRuler: true,
      showSVGButton: false,
      tracks: [
        {
          name: "Ensembl canonical",
          type: "annotation",
          reader,
          displayMode: "EXPANDED",
          height: 90,
          color: "#65765a",
          visibilityWindow: 1000000,
        },
      ],
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
  private async show(
    pane: "gene" | "enhancer",
    ref: ReferenceGene,
    start: number,
    end: number,
    name: string,
    color: string,
  ) {
    const target = locus(
      ref.chrom,
      start,
      end,
      Math.max(1500, (end - start) * 0.15),
    );
    this.targets[pane] = target;
    let browser = this.browsers[pane];
    if (!browser) browser = (await this.create(pane, target)) || undefined;
    if (!browser || this.disposed) return;
    const existing = browser.findTracks("id", "selected-interval");
    existing.forEach((t) => browser!.removeTrack(t));
    const track: TrackLoad<"annotation"> & { id: string } = {
      id: "selected-interval",
      name: "Selected interval",
      type: "annotation",
      features: [
        {
          chr: ref.chrom,
          start,
          end,
          name,
          ...(pane === "gene" ? { strand: ref.strand } : {}),
        },
      ],
      color,
      displayMode: "EXPANDED",
      height: 38,
      order: 0,
    };
    await browser.loadTrack(track);
    await Promise.resolve(browser.search(target));
    byId(`igv-${pane}`).dataset.interval = `${ref.chrom}:${start}-${end}`;
  }
  private async refresh() {
    if (this.busy || !this.active || !this.desired || this.disposed) return;
    this.busy = true;
    try {
      while (this.active && this.desired && !this.disposed) {
        const { record, cre, color } = this.desired;
        const key = JSON.stringify([record.gene.symbol, cre?.id, color]);
        if (key === this.key) break;
        const ref = genomeReference(record);
        byId("igv-error").hidden = true;
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
        const referenceKey = JSON.stringify(ref);
        if (referenceKey !== this.referenceKey) {
          await this.show(
            "gene",
            ref,
            ref.start,
            ref.end,
            record.gene.symbol,
            "#2a6fa8",
          );
          this.referenceKey = referenceKey;
        }
        if (cre) {
          byId("igv-enhancer-title").textContent =
            `Regulatory element · ${elementLabel(record, cre)}`;
          await this.show(
            "enhancer",
            ref,
            cre.start,
            cre.end,
            elementLabel(record, cre),
            color,
          );
        }
        this.key = key;
        byId("browser-panel").dataset.element = cre?.id || "";
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
    Object.values(this.browsers).forEach((b) => igv.removeBrowser(b));
  }
}

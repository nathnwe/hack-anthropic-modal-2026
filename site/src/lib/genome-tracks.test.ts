import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import {
  contextualTracks,
  signalTrackConfig,
  k562TrackCatalogue,
} from "./genome-tracks.ts";
import type { GeneRecord } from "./records.ts";
const load = (): GeneRecord =>
  JSON.parse(
    readFileSync(
      new URL("../../public/data/ldlr.json", import.meta.url),
      "utf8",
    ),
  );

test("K562 experimental signals require a verified matching assembly, cell type and non-placeholder record", () => {
  const record = load();
  assert.equal(contextualTracks(record).length, 4);
  record.locus.cell_type = "HepG2";
  assert.deepEqual(contextualTracks(record), []);
  record.locus.cell_type = "K562";
  record.illustrative = true;
  assert.deepEqual(contextualTracks(record), []);
  record.illustrative = false;
  delete record.analysis;
  assert.deepEqual(contextualTracks(record), []);
  record.caseStudy = { assembly: "hg19" } as GeneRecord["caseStudy"];
  assert.deepEqual(contextualTracks(record), []);
});

test("signal configuration preserves the four pinned ENCODE signals, units and quality provenance", () => {
  const tracks = k562TrackCatalogue.tracks;
  assert.deepEqual(
    tracks.map((t) => t.name),
    ["H3K4me3", "H3K27ac", "H3K27me3", "ATAC-seq"],
  );
  const configs = tracks.map(signalTrackConfig);
  for (const [i, t] of tracks.entries()) {
    assert.equal(t.biosample, "K562");
    assert.equal(t.assembly, "GRCh38");
    assert.equal(t.status, "released");
    assert.equal(t.output_type, "fold change over control");
    assert.ok(t.biological_replicates.length > 1);
    assert.equal(
      t.assay,
      t.name === "ATAC-seq" ? "ATAC-seq" : "Histone ChIP-seq",
    );
    assert.equal(t.target, t.name === "ATAC-seq" ? null : t.name);
    assert.ok(t.url.endsWith(`/${t.accession}.bigWig`));
    assert.equal(configs[i].url, t.url);
    assert.equal(configs[i].min, 0);
    assert.equal(configs[i].windowFunction, "mean");
    assert.ok(configs[i].description.includes(t.output_type));
    assert.match(t.md5sum, /^[0-9a-f]{32}$/);
    assert.ok(Object.keys(t.experiment_audits).length);
  }
  assert.ok(
    tracks
      .find((t) => t.name === "H3K27ac")!
      .experiment_audits.ERROR?.includes("extremely low read depth"),
  );
  assert.ok(
    tracks
      .find((t) => t.name === "H3K27me3")!
      .experiment_audits.NOT_COMPLIANT?.includes("insufficient read depth"),
  );
});

test("cached RefSeq models preserve exon bounds and distinguish uncovered from empty windows", async () => {
  const { referenceGenes } = await import("./reference-genes.ts");
  const snapshot = JSON.parse(
    readFileSync(
      new URL("../../public/data/reference/k562-refseq.json", import.meta.url),
      "utf8",
    ),
  );
  assert.equal(snapshot.assembly, "GRCh38");
  assert.equal(snapshot.windows.length, 6);
  for (const w of snapshot.windows) {
    assert.match(w.source_sha256, /^[0-9a-f]{64}$/);
    assert.ok(w.source_url.includes("ncbiRefSeqCurated"));
    assert.ok(w.features.some((f: { name: string }) => f.name === w.target));
    const identities = new Set();
    for (const f of w.features) {
      assert.ok(!identities.has(`${f.name}:${f.strand}`));
      identities.add(`${f.name}:${f.strand}`);
      assert.equal(f.chr, w.chrom);
      assert.ok(f.start < f.end);
      assert.ok(f.cdStart >= f.start && f.cdEnd <= f.end);
      assert.match(f.id, /^[NX][MR]_/);
      for (const e of f.exons)
        assert.ok(e.start >= f.start && e.end <= f.end && e.end > e.start);
      assert.ok(
        referenceGenes(
          snapshot,
          w.chrom,
          Math.max(f.start, w.start),
          Math.min(f.end, w.end),
        )?.some((r) => r.id === f.id),
      );
    }
    assert.equal(referenceGenes(snapshot, w.chrom, w.start - 1, w.end), null);
  }
  assert.deepEqual(
    referenceGenes(
      {
        assembly: "GRCh38",
        windows: [{ chrom: "chr1", start: 0, end: 10, features: [] }],
      },
      "chr1",
      1,
      2,
    ),
    [],
  );
  assert.equal(
    referenceGenes(
      { ...snapshot, assembly: "hg19" },
      "chr19",
      10500000,
      10510000,
    ),
    null,
  );
});

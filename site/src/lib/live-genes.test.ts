import test from "node:test";
import assert from "node:assert/strict";
import {
  referenceFromEnsembl,
  liveGeneRecord,
  suggestGenes,
  type EnsemblGene,
} from "./live-genes.ts";
import { annotationFeatures, type Feature } from "./annotations.ts";

// Synthetic coordinates used only to exercise coordinate conversion.
const gene: EnsemblGene = {
  assembly_name: "GRCh38",
  species: "homo_sapiens",
  seq_region_name: "8",
  start: 101,
  end: 300,
  strand: 1,
  id: "TEST_GENE",
  display_name: "TEST",
};
test("Ensembl inclusive coordinates become BED without losing either end", () => {
  const plus = referenceFromEnsembl(gene);
  assert.deepEqual(
    [plus.start, plus.end, plus.tss, plus.chrom],
    [100, 300, 100, "chr8"],
  );
  const minus = referenceFromEnsembl({ ...gene, strand: -1 });
  assert.deepEqual(
    [minus.start, minus.end, minus.tss, minus.strand],
    [100, 300, 299, "-"],
  );
  assert.equal(
    referenceFromEnsembl({ ...gene, seq_region_name: "MT" }).chrom,
    "chrM",
  );
});
test("unsupported assemblies or invalid coordinates are rejected", () => {
  for (const change of [
    { assembly_name: "GRCh37" },
    { species: "mus_musculus" },
    { start: 0 },
    { end: 50 },
    { strand: 0 },
  ])
    assert.throws(() => referenceFromEnsembl({ ...gene, ...change }));
});
test("canonical exons and CDS map to the right transcript in IGV coordinates", () => {
  const base: Feature = {
    feature_type: "transcript",
    id: "tx",
    seq_region_name: "8",
    assembly_name: "GRCh38",
    start: 101,
    end: 300,
    strand: -1,
    is_canonical: 1,
  };
  const rows: Feature[] = [
    base,
    { ...base, id: "other", is_canonical: 0 },
    { ...base, feature_type: "exon", Parent: "tx", start: 251, end: 300 },
    { ...base, feature_type: "exon", Parent: "tx", start: 101, end: 150 },
    { ...base, feature_type: "cds", Parent: "tx", start: 121, end: 140 },
    { ...base, feature_type: "exon", Parent: "other", start: 160, end: 190 },
  ];
  const features = annotationFeatures(rows);
  assert.equal(features.length, 1);
  assert.deepEqual(features[0].exons, [
    { start: 100, end: 150 },
    { start: 250, end: 300 },
  ]);
  assert.deepEqual(
    [
      features[0].start,
      features[0].end,
      features[0].cdStart,
      features[0].cdEnd,
      features[0].strand,
    ],
    [100, 300, 120, 140, "-"],
  );
});
test("live lookup uses canonical TSS and never fabricates regulatory predictions", async (t) => {
  const calls: string[] = [];
  t.mock.method(globalThis, "fetch", async (url: string) => {
    calls.push(url);
    return new Response(
      JSON.stringify(
        calls.length === 1
          ? { ...gene, strand: -1, canonical_transcript: "TEST_TX.3" }
          : { ...gene, start: 120, end: 280, strand: -1 },
      ),
      { status: 200 },
    );
  });
  const record = await liveGeneRecord("TEST");
  assert.equal(record.gene.tss, 279);
  assert.equal(record.reference?.end, 300);
  assert.deepEqual(record.locus.cres, []);
  assert.deepEqual(record.staples, []);
  assert.ok(calls[1].includes("/lookup/id/TEST_TX?"));
});
test("suggestions filter species, deduplicate symbols and prefer an exact match", async (t) => {
  t.mock.method(
    globalThis,
    "fetch",
    async () =>
      new Response(
        JSON.stringify({
          hits: [
            { symbol: "TEST_OTHER", name: "other", taxid: 9606 },
            { symbol: "TEST_QUERY", name: "exact", taxid: 9606 },
            { symbol: "TEST_QUERY", name: "duplicate", taxid: 9606 },
            { symbol: "Mouse", taxid: 10090 },
          ],
        }),
        { status: 200 },
      ),
  );
  assert.deepEqual(
    (await suggestGenes("TEST_QUERY")).map((x) => x.symbol),
    ["TEST_QUERY", "TEST_OTHER"],
  );
});

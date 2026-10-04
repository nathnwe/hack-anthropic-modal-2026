import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { Vector3 } from "three";
import {
  dnaCenter,
  dnaFrame,
  detailIntervals,
  locusWindow,
  RISE_NM,
} from "./dna-geometry.ts";
import { k562Record } from "./k562.ts";
import { sortElements } from "./records.ts";
const data = JSON.parse(
  readFileSync(
    new URL("../../public/data/case-studies/myc-k562.json", import.meta.url),
    "utf8",
  ),
);
const record = () => k562Record(data);

test("all 56 K562 source predictions and conditions survive the adapter", () => {
  const r = record();
  assert.equal(r.locus.cres.length, 56);
  assert.equal(new Set(r.locus.cres.map((c) => c.id)).size, 56);
  assert.equal(
    new Set(r.locus.cres.map((c) => `${c.start}:${c.end}`)).size,
    53,
  );
  assert.equal(new Set(r.locus.cres.map((c) => c.source!.accession)).size, 8);
  assert.equal(r.gene.tss, 127736230);
  assert.equal(r.caseStudy!.reference.transcript, "ENST00000621592.8");
  assert.equal(r.locus.contacts.length, 0);
  assert.equal(r.staples.length, 0);
  assert.equal(r.gene.direction, null);
  assert.ok(
    r.locus.cres.some((c) => c.source!.annotation.includes("vorinostat")),
  );
  for (const [i, c] of r.locus.cres.entries()) {
    assert.equal(c.score_rE2G, data.links[i].score);
    assert.equal(c.source!.csv_row, data.links[i].csv_row);
  }
});

test("score and distance sorts select genuinely different first candidates", () => {
  const r = record();
  assert.notEqual(
    sortElements(r, "distance")[0].id,
    sortElements(r, "re2g")[0].id,
  );
  for (const criterion of ["re2g", "distance"]) {
    const sorted = sortElements(r, criterion);
    const metric = (c: (typeof sorted)[number]) =>
      criterion === "re2g"
        ? -(c.score_rE2G ?? -Infinity)
        : Math.abs((c.start + c.end) / 2 - r.gene.tss);
    assert.ok(sorted.every((c, i) => !i || metric(c) >= metric(sorted[i - 1])));
  }
});

test("every K562 candidate has a complete finite continuous view and length-preserving motion", () => {
  const r = record();
  for (const cre of r.locus.cres) {
    const w = locusWindow(r, cre)!;
    assert.ok(w, cre.id);
    assert.ok(
      w.start <= Math.min(cre.start, r.gene.tss) &&
        w.end >= Math.max(cre.end, r.gene.tss),
    );
    let previous = Infinity;
    for (const p of [0, 0.25, 0.5, 0.75, 1]) {
      const a = dnaCenter(w.tss, p, w),
        b = dnaCenter(w.midpoint, p, w);
      assert.ok(a.distanceTo(b) <= previous + 1e-8);
      previous = a.distanceTo(b);
      // Numerically integrate the centreline, independent of its analytic formula.
      let length = 0,
        prev = dnaCenter(w.start, p, w);
      const n = 10000;
      for (let i = 1; i <= n; i++) {
        const point = dnaCenter(w.start + ((w.end - w.start) * i) / n, p, w);
        length += point.distanceTo(prev);
        prev = point;
      }
      const expected = (w.end - w.start) * RISE_NM;
      assert.ok(
        Math.abs(length - expected) / expected < 1e-6,
        `${cre.id}: ${length} vs ${expected}`,
      );
      for (const position of [w.start, w.tss, w.midpoint, w.end]) {
        const f = dnaFrame(position, p, w);
        assert.ok(f.center.toArray().every(Number.isFinite));
        assert.ok(Math.abs(f.normal.dot(f.binormal)) < 1e-10);
        assert.ok(Math.abs(f.binormal.length() - 1) < 1e-10);
      }
    }
    if (!w.overlap) assert.ok(Math.abs(previous - 3) < 1e-6);
  }
});

test("one bp measures 0.34 nm at both short and megabase scales and either direction", () => {
  for (const distance of [1000, 50000, 2000000])
    for (const direction of [-1, 1]) {
      const r = record();
      r.gene.tss = 3000000;
      r.locus.tad = { chrom: "chr8", start: 0, end: 6000000 };
      const c = {
        ...r.locus.cres[0],
        start: r.gene.tss + direction * distance - 20,
        end: r.gene.tss + direction * distance + 21,
      };
      const w = locusWindow(r, c)!;
      assert.equal(
        Math.abs(w.midpoint - w.tss),
        Math.abs(direction * distance + 0.5),
      );
      for (const p of [0, 0.5, 1]) {
        const spacing = dnaCenter(w.midpoint, p, w).distanceTo(
          dnaCenter(w.midpoint + 1, p, w),
        );
        assert.ok(Math.abs(spacing - RISE_NM) < 1e-5);
      }
    }
});

test("molecular zoom includes both contact ends without overlap double-counting", () => {
  const r = record();
  for (const c of r.locus.cres) {
    const w = locusWindow(r, c)!;
    if (w.overlap) continue;
    const center = dnaCenter(w.tss, 1, w)
      .add(dnaCenter(w.midpoint, 1, w))
      .multiplyScalar(0.5);
    const ranges = detailIntervals(w, 1, center, 40);
    for (const position of [w.tss, w.midpoint])
      assert.ok(
        ranges.some((r) => r.start <= position && r.end >= position),
        c.id,
      );
    for (let i = 1; i < ranges.length; i++)
      assert.ok(ranges[i].start > ranges[i - 1].end);
  }
  assert.deepEqual(
    detailIntervals(
      locusWindow(r, r.locus.cres[0])!,
      0,
      new Vector3(1e9, 1e9, 1e9),
      20,
    ),
    [],
  );
});

test("overlap is viewable without a fabricated loop, and invalid coordinates are rejected", () => {
  const r = record(),
    c = { ...r.locus.cres[0], start: r.gene.tss - 100, end: r.gene.tss + 100 };
  const w = locusWindow(r, c)!;
  assert.ok(w.overlap);
  assert.deepEqual(dnaCenter(w.start, 0, w), dnaCenter(w.start, 1, w));
  r.gene.chrom = "chr9";
  assert.equal(locusWindow(r, c), null);
  r.gene.chrom = "chr8";
  r.gene.tss += 0.5;
  assert.equal(locusWindow(r, c), null);
});

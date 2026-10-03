import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import {
  dnaCenter,
  helixFrames,
  locusWindow,
  helixDisplay,
  helixSamples,
  displayCenter,
  displayFrames,
} from "./dna-geometry.ts";
import type { GeneRecord } from "./records.ts";
const fixture = (): GeneRecord =>
  JSON.parse(
    readFileSync(
      new URL("../../public/data/myc.json", import.meta.url),
      "utf8",
    ),
  );

test("genomic highlights follow source positions in either genomic order", () => {
  const r = fixture(),
    cre = r.locus.cres[0];
  const window = locusWindow(r, cre)!;
  assert.equal(
    window.start + window.promoter * (window.end - window.start),
    r.gene.tss,
  );
  assert.equal(
    window.start + window.element * (window.end - window.start),
    (cre.start + cre.end) / 2,
  );
  r.gene.tss = 1_400_000;
  const upstream = locusWindow(r, cre)!;
  assert.ok(upstream.element < upstream.promoter);
  assert.ok(
    upstream.start >= r.locus.tad.start && upstream.end <= r.locus.tad.end,
  );
});

test("missing, overlapping or out-of-domain anchors cannot fabricate a 3D preview", () => {
  const r = fixture(),
    cre = r.locus.cres[0];
  r.gene.tss = cre.start;
  assert.equal(locusWindow(r, cre), null);
  r.gene.tss = r.locus.tad.end + 1;
  assert.equal(locusWindow(r, cre), null);
  r.gene.tss = NaN;
  assert.equal(locusWindow(r, cre), null);
});

test("contact preview brings the supplied anchors closer without disconnecting the DNA", () => {
  const r = fixture(),
    window = locusWindow(r, r.locus.cres[0])!;
  let previous = Infinity;
  for (const progress of [0, 0.25, 0.5, 0.75, 1]) {
    const distance = dnaCenter(window.promoter, progress, window).distanceTo(
      dnaCenter(window.element, progress, window),
    );
    assert.ok(distance < previous);
    previous = distance;
    for (const t of [window.promoter, window.element]) {
      assert.ok(
        dnaCenter(t - 1e-6, progress, window).distanceTo(
          dnaCenter(t + 1e-6, progress, window),
        ) < 0.001,
      );
    }
  }
  assert.ok(previous > 2, "keep the highlighted helix centerlines distinct");
});

test("helix frames remain finite and orthonormal throughout the deformation", () => {
  const r = fixture(),
    window = locusWindow(r, r.locus.cres[0])!;
  for (const progress of [0, 0.25, 0.5, 0.75, 1]) {
    const frames = helixFrames(window, progress, 300);
    for (let i = 0; i < frames.centers.length; i++) {
      assert.ok(frames.centers[i].toArray().every(Number.isFinite));
      assert.ok(Math.abs(frames.normals[i].length() - 1) < 1e-6);
      assert.ok(Math.abs(frames.binormals[i].length() - 1) < 1e-6);
      assert.ok(Math.abs(frames.normals[i].dot(frames.binormals[i])) < 1e-6);
    }
  }
});

// Synthetic coordinates exercise software scale limits, not biological examples.
test("short and megabase-scale views account for every displayed or omitted base", () => {
  for (const separation of [10, 100, 1_000, 50_000, 500_000, 2_000_000]) {
    for (const direction of [-1, 1]) {
      const r = fixture();
      r.gene.tss = 3_000_000;
      r.locus.tad.start = 0;
      r.locus.tad.end = 6_000_000;
      const midpoint = r.gene.tss + direction * separation;
      const cre = {
        ...r.locus.cres[0],
        start: midpoint - 1,
        end: midpoint + 1,
      };
      const display = helixDisplay(locusWindow(r, cre)!);
      const samples = helixSamples(display);
      assert.equal(display.separationBp, separation);
      assert.equal(
        samples.flat().filter((s) => s.rung).length,
        display.shownBp,
      );
      assert.ok(
        display.shownBp <= 180,
        "bounded geometry even at megabase scale",
      );
      assert.equal(
        display.shownBp + display.omittedBp,
        display.segments.at(-1)!.end - display.segments[0].start,
      );
      for (let j = 0; j < samples.length; j++) {
        assert.ok(
          samples[j].every(
            (s) =>
              s.position >= display.segments[j].start &&
              s.position < display.segments[j].end,
          ),
        );
        const rungs = samples[j].filter((s) => s.rung);
        for (let i = 1; i < rungs.length; i++)
          assert.equal(rungs[i].position - rungs[i - 1].position, 1);
      }
      assert.equal(display.omittedBp === 0, separation <= 100);
      let previousDistance = Infinity;
      for (const p of [0, 0.25, 0.5, 0.75, 1]) {
        const distance = displayCenter(
          display.shape.promoter,
          p,
          display,
        ).distanceTo(displayCenter(display.shape.element, p, display));
        assert.ok(distance < previousDistance);
        previousDistance = distance;
        const frames = displayFrames(display, p, 300);
        assert.ok(
          frames.centers.every((v) => v.toArray().every(Number.isFinite)),
        );
        assert.ok(frames.normals.every((v) => Math.abs(v.length() - 1) < 1e-6));
        if (display.gap)
          assert.ok(distance > 3, "local duplexes remain visually distinct");
      }
    }
  }
});

test("odd-width elements and clipped boundaries preserve exact separation and omission counts", () => {
  const r = fixture();
  r.gene.tss = r.locus.tad.start;
  const cre = {
    ...r.locus.cres[0],
    start: r.gene.tss + 1_000,
    end: r.gene.tss + 1_003,
  };
  const display = helixDisplay(locusWindow(r, cre)!);
  assert.equal(display.separationBp, 1_001.5);
  assert.equal(display.segments[0].start, r.locus.tad.start);
  assert.equal(
    display.omittedBp,
    display.segments[1].start - display.segments[0].end,
  );
  assert.ok(
    display.segments.every(
      (s) => Number.isInteger(s.start) && Number.isInteger(s.end),
    ),
  );
});

test("mismatched chromosomes and fractional genomic coordinates cannot create a cutaway", () => {
  const r = fixture();
  r.gene.chrom = "SYNTHETIC_OTHER_CHROMOSOME";
  assert.equal(locusWindow(r, r.locus.cres[0]), null);
  r.gene.chrom = r.locus.tad.chrom.replace(/^chr/, "");
  assert.ok(locusWindow(r, r.locus.cres[0]));
  r.gene.tss += 0.5;
  assert.equal(locusWindow(r, r.locus.cres[0]), null);
});

import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { dnaCenter, helixFrames, locusWindow } from "./dna-geometry.ts";
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

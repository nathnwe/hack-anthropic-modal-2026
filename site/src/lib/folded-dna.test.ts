import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import {
  foldedDNA,
  foldedFrame,
  foldedAnchor,
  hazeCenter,
} from "./folded-dna.ts";
import { locusWindow } from "./dna-geometry.ts";
import { k562Record, K562_VIEWER_IDS } from "./k562.ts";
const data = JSON.parse(
  readFileSync(
    new URL("../../public/data/case-studies/myc-k562.json", import.meta.url),
    "utf8",
  ),
);
const record = k562Record(data, K562_VIEWER_IDS);

test("three sourced examples have exact, distinct omission counts and bounded local windows", () => {
  assert.equal(record.locus.cres.length, 3);
  const expected: Record<string, number> = {
    "K562-49": 162394,
    "K562-51": 6006,
    "K562-55": 1845947,
  };
  const displays = record.locus.cres
    .map((cre) => {
      const d = foldedDNA(locusWindow(record, cre)!)!;
      assert.equal(d.omittedBp, expected[cre.id]);
      assert.equal(d.shownBp, 256);
      assert.equal(d.shownBp + d.omittedBp, d.end - d.start);
      assert.equal(d.arms[0].end + d.omittedBp, d.arms[1].start);
      assert.ok(
        d.arms.some(
          (a) => a.start <= record.gene.tss && a.end > record.gene.tss,
        ),
      );
      assert.ok(
        d.arms.some(
          (a) =>
            a.start <= (cre.start + cre.end) / 2 &&
            a.end > (cre.start + cre.end) / 2,
        ),
      );
      return d;
    })
    .sort((a, b) => a.omittedBp - b.omittedBp);
  assert.ok(
    displays.every((d, i) => !i || d.hazeWidth > displays[i - 1].hazeWidth),
  );
  assert.ok(
    displays.every((d, i) => !i || d.separation > displays[i - 1].separation),
  );
});

test("fade is progressive at each omitted boundary, while marked sites remain opaque", () => {
  for (const cre of record.locus.cres) {
    const d = foldedDNA(locusWindow(record, cre)!)!;
    for (const a of d.arms) {
      const atDistance = (bp: number) =>
        a.side < 0 ? a.end - bp : a.start + bp;
      const fades = [0, 0.5, 8, 16, 24, 34].map(
        (bp) => foldedFrame(atDistance(bp), 0, d).fade,
      );
      assert.equal(fades[0], 0);
      assert.equal(fades.at(-1), 1);
      assert.ok(fades.every((f, i) => !i || f > fades[i - 1]));
      assert.equal(foldedFrame(a.anchor, 0, d).fade, 1);
    }
    assert.throws(
      () => foldedFrame((d.arms[0].end + d.arms[1].start) / 2, 0, d),
      /omitted interval/,
    );
  }
});

test("fluid motion keeps the promoter fixed and preserves local spacing and exact accounting", () => {
  for (const cre of record.locus.cres) {
    const d = foldedDNA(locusWindow(record, cre)!)!;
    let previous = Infinity;
    for (const p of [0, 0.25, 0.5, 0.75, 1]) {
      const gap = foldedAnchor(d.arms[0].anchor, p, d).distanceTo(
        foldedAnchor(d.arms[1].anchor, p, d),
      );
      assert.ok(gap < previous);
      previous = gap;
      assert.ok(hazeCenter(d, p).toArray().every(Number.isFinite));
      assert.deepEqual(
        foldedAnchor(record.gene.tss, p, d).toArray(),
        [0, -8, 0],
      );
      for (const a of d.arms)
        for (let bp = a.start + 0.5; bp < a.end; bp++) {
          const f = foldedFrame(bp, p, d);
          assert.ok(f.center.toArray().every(Number.isFinite));
          assert.ok(
            Math.abs(f.normal.length() - 1) < 1e-6 &&
              Math.abs(f.binormal.length() - 1) < 1e-6,
          );
          assert.ok(Math.abs(f.normal.dot(f.binormal)) < 1e-6);
          if (bp + 1 < a.end) {
            const spacing = f.center.distanceTo(
              foldedFrame(bp + 1, p, d).center,
            );
            assert.ok(Math.abs(spacing - 0.34) < 0.004);
          }
        }
    }
    assert.ok(previous > 3 && previous < 4);
    assert.equal(d.shownBp + d.omittedBp, d.end - d.start);
  }
});

test("enhancer follows a curved path and the DNA bends instead of moving as a rigid body", () => {
  const d = foldedDNA(locusWindow(record, record.locus.cres[0])!)!;
  const a = d.arms.find((a) => a.role === "element")!;
  const start = foldedAnchor(a.anchor, 0, d),
    end = foldedAnchor(a.anchor, 1, d);
  const midpoint = foldedAnchor(a.anchor, 0.5, d);
  assert.ok(midpoint.distanceTo(start.clone().lerp(end, 0.5)) > 3);
  const chord = (p: number) =>
    foldedAnchor(a.start + 2, p, d).distanceTo(foldedAnchor(a.end - 2, p, d));
  assert.ok(Math.abs(chord(0) - chord(1)) > 2);
  // Smooth start and settling, with no time-dependent jitter while paused.
  assert.ok(foldedAnchor(a.anchor, 0.001, d).distanceTo(start) < 1e-5);
  assert.ok(foldedAnchor(a.anchor, 0.999, d).distanceTo(end) < 1e-5);
  const paused = foldedFrame(a.anchor + 30, 0.4, d).center.toArray();
  foldedFrame(a.anchor + 30, 0.9, d);
  assert.deepEqual(foldedFrame(a.anchor + 30, 0.4, d).center.toArray(), paused);
});

test("upstream elements, half-base midpoints and clipped flanks retain accounting", () => {
  const r = k562Record(data);
  r.gene.tss = 127899000;
  const cre = { ...r.locus.cres[0], start: 127737000, end: 127737501 };
  const w = locusWindow(r, cre)!;
  const d = foldedDNA(w)!;
  assert.equal(d.arms[0].anchor, 127737250.5);
  assert.equal(d.arms[1].anchor, r.gene.tss);
  for (const p of [0, 0.5, 1]) {
    assert.deepEqual(foldedAnchor(r.gene.tss, p, d).toArray(), [0, -8, 0]);
    assert.ok(foldedAnchor(cre.start + 250.5, p, d).y > -8);
  }
  const clipped = foldedDNA({
    ...w,
    start: Math.floor(d.arms[0].anchor) - 10,
  })!;
  assert.equal(clipped.arms[0].start, 127737240);
  assert.equal(
    clipped.omittedBp + clipped.shownBp,
    clipped.end - clipped.start,
  );
  assert.equal(foldedDNA({ ...w, overlap: true }), null);
});

test("visible local helices stay separate throughout the three contact previews", () => {
  for (const cre of record.locus.cres) {
    const d = foldedDNA(locusWindow(record, cre)!)!;
    for (let step = 0; step <= 20; step++) {
      const samples = d.arms.map((a) =>
        Array.from({ length: 64 }, (_, i) =>
          foldedFrame(a.start + i * 2 + 0.5, step / 20, d),
        ).filter((f) => f.fade > 0.5),
      );
      for (const a of samples[0])
        for (const b of samples[1])
          assert.ok(
            a.center.distanceTo(b.center) > 2.4,
            "Visible double helices intersect",
          );
    }
  }
});

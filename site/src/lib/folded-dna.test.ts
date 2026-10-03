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

test("folded motion preserves local molecular geometry and exact accounting", () => {
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
            const initial = foldedFrame(bp, 0, d).center.distanceTo(
              foldedFrame(bp + 1, 0, d).center,
            );
            assert.ok(Math.abs(spacing - initial) < 1e-8);
            assert.ok(Math.abs(spacing - 0.34) < 0.004);
          }
        }
    }
    assert.ok(previous > 3 && previous < 4);
    assert.equal(d.shownBp + d.omittedBp, d.end - d.start);
  }
});

test("upstream elements, half-base midpoints and clipped flanks retain accounting", () => {
  const r = k562Record(data);
  r.gene.tss = 127899000;
  const cre = { ...r.locus.cres[0], start: 127737000, end: 127737501 };
  const w = locusWindow(r, cre)!;
  const d = foldedDNA(w)!;
  assert.equal(d.arms[0].anchor, 127737250.5);
  assert.equal(d.arms[1].anchor, r.gene.tss);
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

import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import {
  rankingColors,
  SCORE_HIGH,
  SCORE_LOW,
  SCORE_MISSING,
} from "./element-colors.ts";
import type { GeneRecord } from "./records.ts";
const fixture = (): GeneRecord =>
  JSON.parse(
    readFileSync(
      new URL("../../public/data/myc.json", import.meta.url),
      "utf8",
    ),
  );

test("colours follow numerical score deltas, not equal steps between ranks", () => {
  const record = fixture();
  const original = record.locus.cres[0];
  record.locus.cres = [0.3, 0.36, 0.9].map((score, i) => ({
    ...original,
    id: `test-${i}`,
    score_rE2G: score,
  }));
  const scale = rankingColors(record, "re2g");
  assert.equal(scale.color(record.locus.cres[0]), SCORE_LOW);
  assert.equal(scale.color(record.locus.cres[2]), SCORE_HIGH);
  // 0.36 is 10% of the way across this range, not its midpoint.
  assert.equal(scale.color(record.locus.cres[1]), "#e9d097");
});
test("a distance sort reverses colour direction; equal scores have equal colours", () => {
  const record = fixture();
  record.locus.cres.forEach((c) => (c.score_rE2G = 0.5));
  const same = rankingColors(record, "re2g");
  assert.equal(
    same.color(record.locus.cres[0]),
    same.color(record.locus.cres[1]),
  );
  const scale = rankingColors(record, "distance");
  const order = [...record.locus.cres].sort(
    (a, b) =>
      Math.abs((a.start + a.end) / 2 - record.gene.tss) -
      Math.abs((b.start + b.end) / 2 - record.gene.tss),
  );
  assert.equal(scale.color(order[0]), SCORE_HIGH);
  assert.equal(scale.color(order.at(-1)!), SCORE_LOW);
});
test("missing contacts remain neutral while recorded zero is scored", () => {
  const record = fixture();
  record.locus.contacts = [
    { a: record.gene.tss, b: record.locus.cres[0].start, strength: 0 },
  ];
  const scale = rankingColors(record, "contact");
  assert.equal(scale.min, 0);
  assert.notEqual(scale.color(record.locus.cres[0]), SCORE_MISSING);
  assert.equal(scale.color(record.locus.cres[1]), SCORE_MISSING);
  record.locus.cres = [];
  assert.equal(rankingColors(record, "re2g").min, null);
});

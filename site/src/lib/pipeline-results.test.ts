import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { createHash } from "node:crypto";
import type { GeneRecord } from "./records.ts";
import {
  assertPipelineRecord,
  shortlist,
  rankedEnhancers,
  designElement,
  designColor,
  responseLabel,
  riskLabel,
} from "./pipeline-results.ts";
import { locusWindow } from "./dna-geometry.ts";
import { foldedDNA, foldedAnchor, foldedFrame } from "./folded-dna.ts";
import { selectRecord, type RecordEntry } from "./analysis-routing.ts";

const expected = { ldlr: 1, nsd1: 3, runx1: 1, nf1: 4, mecp2_up: 2, mecp2_down: 1, myc_down: 1 };
const record = (key: string): GeneRecord =>
  JSON.parse(
    readFileSync(
      new URL(`../../public/data/${key}.json`, import.meta.url),
      "utf8",
    ),
  );

test("seven exports preserve source bytes, native ranks, sweeps and risk semantics without fake contract scores", () => {
  for (const [key, count] of Object.entries(expected)) {
    const r = record(key),
      a = r.analysis!;
    assertPipelineRecord(r);
    const bytes = readFileSync(
      new URL(`../../public/data/${a.source_snapshot}`, import.meta.url),
    );
    assert.equal(
      a.source_sha256,
      createHash("sha256").update(bytes).digest("hex"),
    );
    // JSON serialization normalizes the source's rounded -0.0 to 0.
    const source = JSON.parse(JSON.stringify(JSON.parse(bytes.toString())));
    assert.deepEqual(a.source, source);
    assert.equal(shortlist(a).length, count);
    assert.deepEqual(
      shortlist(a).map((d) => d.rank),
      source.ranked.map((d: { rank: number }) => d.rank),
    );
    for (const row of [...source.ranked, ...source.best_below_floor]) {
      const d = a.designs.find(
        (d) =>
          d.anchor_a_kb === row.anchor_a_kb &&
          d.anchor_b_kb === row.anchor_b_kb,
      )!;
      for (const key of Object.keys(row))
        assert.deepEqual(d[key as keyof typeof d], row[key]);
      assert.equal(designElement(d).score_rE2G, null);
      assert.equal(d.element.start, row.enhancer_kb * 1000);
      assert.equal(d.anchor_a.start, row.anchor_a_kb * 1000);
      assert.equal(d.anchor_b.start, row.anchor_b_kb * 1000);
    }
    assert.deepEqual(r.locus.cres, []);
    assert.deepEqual(r.locus.contacts, []);
    assert.deepEqual(r.ranking, []);
    assert.equal(r.staples.length, count);
    assert.ok(
      r.staples.every(
        (d) =>
          d.mode === source.direction &&
          d.predicted_delta === `${source.direction === "down" ? "decrease" : "increase"}_at_reference_strength`,
      ),
    );
    assert.equal(
      a.omitted_nonranked,
      source.n_staples_scored - a.designs.length,
    );
  }
});

test("NSD1 reference response cannot hide its opposite-direction sweep or promote its blocked design", () => {
  const r = record("nsd1"),
    a = r.analysis!,
    d = shortlist(a)[0];
  assert.equal(d.delta_target, 0.4002);
  assert.deepEqual(d.sweep, [0, -0.0003, -0.0034, -0.0082, -0.0052]);
  assert.equal(responseLabel(d), "Increase only at the reference setting");
  const blocked = a.designs.find((d) => d.risk.tier === "BLOCKED")!;
  assert.equal(blocked.risk.worst_gene, "DDX41");
  assert.equal(blocked.status, "blocked");
  assert.ok(!r.staples.some((s) => s.id === blocked.id));
  for (const sort of ["pipeline", "distance", "risk"] as const)
    assert.ok(!shortlist(a, sort).some((d) => d.id === blocked.id));
  blocked.status = "shortlisted";
  blocked.rank = 1;
  blocked.delta_target = 10;
  assert.throws(() => assertPipelineRecord(r), /Rejected design/);
});

test("CLEAR remains an incomplete screen rather than a safety verdict; absent sweeps remain absent", () => {
  const r = record("runx1"),
    a = r.analysis!,
    d = shortlist(a)[0];
  assert.equal(d.risk.tier, "CLEAR");
  assert.equal(riskLabel(d.risk), "No other ClinGen genes");
  assert.deepEqual(r.risks[0].flags, []);
  assert.equal(responseLabel(d), "Increase across the tested strengths");
  assert.ok(
    a.designs
      .filter((d) => d.status !== "shortlisted")
      .every((d) => d.sweep === null),
  );
  assert.equal(a.omitted_nonranked, 4);
});

test("sorts use native metrics, retain ranks, and colour actual relative deltas", () => {
  const a = record("nf1").analysis!,
    original = a.designs.map((d) => d.id);
  for (const sort of ["pipeline", "distance", "risk"] as const) {
    const rows = shortlist(a, sort);
    const values = rows.map((d) =>
      sort === "pipeline"
        ? d.rank!
        : sort === "distance"
          ? d.separation_kb
          : d.risk.score,
    );
    assert.ok(values.every((v, i) => !i || v >= values[i - 1]));
    if (sort === "pipeline") {
      assert.equal(designColor(rows, rows[0], sort), "#ad5641");
      assert.equal(designColor(rows, rows.at(-1)!, sort), "#f0dea0");
    }
  }
  assert.deepEqual(
    a.designs.map((d) => d.id),
    original,
  );
  a.designs = [];
  assert.deepEqual(shortlist(a), []);
});

test("all pipeline anchor pairs have finite four-window geometry with exact bp accounting and actual anchor tethering", () => {
  for (const key of Object.keys(expected)) {
    const r = record(key);
    for (const d of r.analysis!.designs) {
      const cre = designElement(d),
        locus = locusWindow(r, cre)!;
      locus.start = Math.min(locus.start, d.anchor_a.start, d.anchor_b.start);
      locus.end = Math.max(locus.end, d.anchor_a.end, d.anchor_b.end);
      const anchors = {
        gene: (d.anchor_b.start + d.anchor_b.end) / 2,
        element: (d.anchor_a.start + d.anchor_a.end) / 2,
      };
      const folded = foldedDNA(locus, anchors);
      if (locus.overlap) {
        assert.equal(folded, null, "Overlapping target/element keeps the static full-span fallback");
        continue;
      }
      assert.ok(folded, `${d.id} lacks an anchor-aware molecular view`);
      assert.equal(
        folded.shownBp + folded.omittedBp,
        folded.end - folded.start,
      );
      // Check the union independently: no source base may be duplicated when
      // a decoy sits beyond the enhancer and the windows change arms.
      const windows = [[r.gene.tss - 36, r.gene.tss + 36], [locus.midpoint - 36, locus.midpoint + 36],
        [anchors.gene - 64, anchors.gene + 64], [anchors.element - 64, anchors.element + 64]].sort((a,b) => a[0]-b[0]);
      let covered = 0, end = -Infinity;
      for (const [lo, hi] of windows) { covered += Math.max(0, hi - Math.max(end, lo)); end = Math.max(end, hi); }
      assert.equal(folded.shownBp, covered);
      for (const arm of folded.arms) {
        const [a, b] = arm.display!.windows;
        if (b)
          assert.throws(
            () => foldedFrame((a.end + b.start) / 2, 0, folded),
            /compressed local gap/,
          );
      }
      for (const p of [0, 0.5, 1]) {
        for (const point of [
          r.gene.tss,
          locus.midpoint,
          anchors.gene,
          anchors.element,
        ])
          assert.ok(
            foldedAnchor(point, p, folded).toArray().every(Number.isFinite),
          );
        // The lower anchor remains fixed, while enhancer-side DNA moves towards it.
        assert.deepEqual(
          foldedAnchor(anchors.gene, p, folded).toArray(),
          [0, -8, 0],
        );
        for (const arm of folded.arms)
          for (const w of arm.display!.windows) {
            const mid = (w.start + w.end) / 2;
            const spacing = foldedFrame(mid, p, folded).center.distanceTo(
              foldedFrame(mid + 1, p, folded).center,
            );
            assert.ok(Math.abs(spacing - 0.34) < 0.004);
          }
      }
      assert.ok(
        Math.abs(
          foldedAnchor(anchors.element, 1, folded).distanceTo(
            foldedAnchor(anchors.gene, 1, folded),
          ) - 3.5,
        ) < 1e-8,
      );
      assert.ok(
        foldedAnchor(r.gene.tss, 1, folded).distanceTo(
          foldedAnchor(anchors.gene, 1, folded),
        ) > 1,
      );
    }
  }
});

test("new downregulation records preserve decoys, negative effects, p99 and collateral without reversing the goal", () => {
  for (const key of ["mecp2_down", "myc_down"]) {
    const r = record(key), a = r.analysis!, d = rankedEnhancers(a)[0];
    assertPipelineRecord(r);
    assert.equal(r.gene.direction, "too_much");
    assert.deepEqual(a.supported_modes, ["down"]);
    assert.equal(d.mode, "down");
    assert.ok(d.delta_target <= -0.10);
    assert.equal(d.delta_at_p99, d.sweep!.at(-1));
    assert.equal(responseLabel(d), "Decrease across the tested strengths");
    assert.equal(d.anchor_b.start, d.decoy_kb! * 1000);
    // Distinct decoys for the same enhancer must order by reduction, not raw
    // signed number (which would favour the smallest reduction).
    a.designs = [structuredClone(d), {...structuredClone(d), id: "stronger-decoy", rank: 2, delta_target: -0.2}];
    const rows = rankedEnhancers(a);
    assert.equal(rows[0].id, "stronger-decoy");
    assert.equal(designColor(rows, rows[0], "pipeline"), "#ad5641");
    a.designs[0].delta_target = 0.2;
    assert.throws(() => assertPipelineRecord(r), /Rejected design/);
  }
  const down = record("mecp2_down").analysis!.designs[0];
  assert.equal(down.risk.tier, "ELEVATED");
  assert.equal(down.risk.worst_gene, "FLNA");
  assert.equal(down.unflagged_collateral_over_tau.find((c) => c.gene === "TEX28")!.delta, 1.0002);
  const up = record("mecp2_up");
  const rejected = up.analysis!.designs.find((d) => d.risk.tier === "BLOCKED")!;
  assert.equal(rejected.risk.worst_gene, "AVPR2");
  assert.ok(!up.staples.some((s) => s.id === rejected.id));
  assert.equal(up.analysis!.source.dropped.blocked_by_risk, 1);
  assert.equal(record("nsd1").analysis!.source.dropped.blocked_by_risk, 0);
});

test("routing distinguishes both MECP2 directions and preserves explicit legacy MYC links", () => {
  const entries: RecordEntry[] = Object.keys(expected).map((key) => {
    const r = record(key);
    return {symbol:r.gene.symbol, file:`${key}.json`, direction:r.gene.direction, analysis:true, modes:r.analysis!.supported_modes};
  });
  entries.push({symbol:"MYC",file:"myc.json",direction:"too_much",analysis:false});
  assert.equal(selectRecord(entries,"MECP2","up","analysis")!.file,"mecp2_up.json");
  assert.equal(selectRecord(entries,"MECP2","down","analysis")!.file,"mecp2_down.json");
  assert.equal(selectRecord(entries,"MYC","down","analysis")!.file,"myc_down.json");
  assert.equal(selectRecord(entries,"MYC","up",null)!.file,"myc.json");
  assert.equal(selectRecord(entries,"MYC","up","sample")!.file,"myc.json");
  assert.equal(selectRecord(entries,"MYC","up","analysis")!.file,"myc_down.json");
  assert.equal(selectRecord(entries,"BRCA1","up","reference"),undefined);
});

test("enhancer ordering follows reference delta, never native rank or rejected high scores", () => {
  const a = record("nf1").analysis!;
  const eligible = shortlist(a);
  const originalIds = a.designs.map((d) => d.id);
  // Deliberately disagree with the saved rank to detect a rank-based UI sort.
  eligible.forEach((d, i) => {
    d.delta_target = (i + 1) / 10;
  });
  const rejected = a.designs.find((d) => d.status !== "shortlisted")!;
  rejected.delta_target = 100;
  const rows = rankedEnhancers(a);
  assert.deepEqual(
    rows.map((d) => d.id),
    eligible.toReversed().map((d) => d.id),
  );
  assert.ok(!rows.some((d) => d.id === rejected.id));
  assert.deepEqual(
    a.designs.map((d) => d.id),
    originalIds,
  );
  assert.deepEqual(
    rows.map((d) => d.rank),
    eligible.toReversed().map((d) => d.rank),
  );
});

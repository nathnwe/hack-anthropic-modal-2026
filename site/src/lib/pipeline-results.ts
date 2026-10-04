import type { CRE, GeneRecord, ReferenceGene, Region } from "./records.ts";

export type PipelineRisk = {
  tier: "CLEAR" | "LOW" | "ELEVATED" | "BLOCKED";
  score: number;
  worst_gene: string | null;
  n_flagged_in_window: number;
  reason: string;
};
export type PipelineDesign = {
  mode: "up" | "down";
  delta_at_p99: number | null;
  decoy_kb?: number;
  id: string;
  name: string;
  target: string;
  rank: number | null;
  status: "shortlisted" | "blocked" | "below_floor";
  enhancer_kb: number;
  anchor_a_kb: number;
  anchor_b_kb: number;
  anchor_a: Region;
  anchor_b: Region;
  element: Region & { id: string; type: "enhancer" };
  separation_kb: number;
  delta_target: number;
  q: number;
  baseline_T: number;
  target_tpm: number;
  target_flag: string | null;
  sweep: number[] | null;
  sign_stable: boolean;
  confidence: "CONFIDENT" | "TENTATIVE" | "ABSTAIN";
  risk: PipelineRisk;
  unflagged_collateral_over_tau: { gene: string; delta: number }[];
};
export type PipelineAnalysis = {
  version: 1;
  run_id: string;
  source_path: string;
  source_sha256: string;
  source_snapshot: string;
  source_revision: string;
  coordinate_system: string;
  window_kind: "analysis_window";
  reference: ReferenceGene;
  annotation_retrieved: string;
  metric: "relative_regulatory_input_change";
  supported_modes: ("up" | "down")[];
  sweep_percentiles: number[];
  designs: PipelineDesign[];
  omitted_nonranked: number;
  limitations: string[];
  source: {
    window: string;
    n_staples_scored: number;
    n_meeting_dose_goal: number;
    dropped: Record<string, number>;
    tests_evaluated: Record<string, string>;
    config: { MIN_EFFECT: number; MAX_CONFIDENT_DISTANCE_KB: number; [key: string]: unknown };
    flagged_genes_in_window: Record<string, string>;
    anchor_space: {
      excluded_genes: string[];
      excluded_flagged: string[];
      permitted_anchors: number;
    };
    provenance: {
      retrieved: string;
      hic: {
        file: string;
        experiment: string;
        cell: string;
        assembly: string;
        balancing: string;
      };
      dnase: { experiment: string; peaks: string; signal: string };
      genes: string;
      rna: string;
      ctcf: string;
      dosage: string;
      not_covered: string[];
    };
  };
};
export type DesignSort = "pipeline" | "distance" | "risk";
export function shortlist(
  analysis: PipelineAnalysis,
  sort: DesignSort = "pipeline",
) {
  const rows = analysis.designs.filter(
    (d) =>
      d.status === "shortlisted" &&
      d.risk.tier !== "BLOCKED" &&
      d.confidence !== "ABSTAIN",
  );
  return rows.sort(
    (a, b) =>
      (sort === "distance"
        ? a.separation_kb - b.separation_kb
        : sort === "risk"
          ? a.risk.score - b.risk.score
          : 0) ||
      (a.rank ?? Infinity) - (b.rank ?? Infinity) ||
      a.id.localeCompare(b.id),
  );
}
export function designElement(design: PipelineDesign): CRE {
  return {
    id: design.element.id,
    start: design.element.start,
    end: design.element.end,
    type: "enhancer",
    score_rE2G: null,
  };
}
// The hit is presented by its enhancer; its single associated design supplies
// the reference-setting response. Never substitute an intrinsic enhancer score.
export function rankedEnhancers(analysis: PipelineAnalysis) {
  return shortlist(analysis).sort(
    (a, b) =>
      directedEffect(b) - directedEffect(a) ||
      a.rank! - b.rank! ||
      a.id.localeCompare(b.id),
  );
}
export const directedEffect = (d: PipelineDesign) =>
  (d.mode === "down" ? -1 : 1) * d.delta_target;
export function riskLabel(risk: PipelineRisk) {
  return {
    CLEAR: "No other ClinGen genes",
    LOW: "Low screened impact",
    ELEVATED: "Elevated screened impact",
    BLOCKED: "Rejected by risk gate",
  }[risk.tier];
}
export function responseLabel(d: PipelineDesign) {
  const change = d.mode === "down" ? "Decrease" : "Increase";
  if (d.status === "blocked") return "Rejected by the risk gate";
  if (d.status === "below_floor")
    return directedEffect(d) < 0
      ? "Opposite to the requested direction"
      : "Below the model’s 10% threshold";
  if (d.sweep && d.sweep.every((v) => (d.mode === "down" ? -v : v) <= 0))
    return `${change} only at the reference setting`;
  if (!d.sign_stable) return `${change} depends on contact strength`;
  return `${change} across the tested strengths`;
}
export function confidenceReasons(a: PipelineAnalysis, d: PipelineDesign) {
  const reasons: string[] = [];
  if (d.separation_kb > a.source.config.MAX_CONFIDENT_DISTANCE_KB)
    reasons.push(
      `Anchor separation exceeds the ${a.source.config.MAX_CONFIDENT_DISTANCE_KB} kb confidence limit.`,
    );
  if (!d.sign_stable)
    reasons.push(
      d.sweep
        ? `A ${d.mode === "down" ? "decrease" : "increase"} is not sustained at every tested contact strength.`
        : "Contact-strength sweep was not run for this design.",
    );
  if (!reasons.length)
    reasons.push(
      "Distance and direction-stability checks pass. This is model confidence, not experimental validation.",
    );
  return reasons;
}
export function percent(value: number, digits = 1) {
  const rounded = Number((value * 100).toFixed(digits));
  return `${rounded > 0 ? "+" : ""}${Object.is(rounded, -0) ? 0 : rounded}%`;
}
export function designColor(
  rows: PipelineDesign[],
  d: PipelineDesign,
  sort: DesignSort,
) {
  if (d.status !== "shortlisted") return "#81796e";
  const value = (v: PipelineDesign) =>
    sort === "distance"
      ? v.separation_kb
      : sort === "risk"
        ? v.risk.score
        : directedEffect(v);
  const values = rows.map(value),
    min = Math.min(...values),
    max = Math.max(...values);
  let t = min === max ? 1 : (value(d) - min) / (max - min);
  if (min !== max && sort !== "pipeline") t = 1 - t;
  return (
    "#" +
    [0, 1, 2]
      .map((i) => {
        const lo = parseInt("#f0dea0".slice(1 + i * 2, 3 + i * 2), 16),
          hi = parseInt("#ad5641".slice(1 + i * 2, 3 + i * 2), 16);
        return Math.round(lo + (hi - lo) * t)
          .toString(16)
          .padStart(2, "0");
      })
      .join("")
  );
}
export function assertPipelineRecord(record: GeneRecord) {
  const a = record.analysis;
  if (
    !a ||
    a.version !== 1 ||
    a.metric !== "relative_regulatory_input_change" ||
    a.reference.assembly !== "GRCh38"
  )
    throw new Error("Unsupported pipeline analysis.");
  if (record.gene.tss !== a.reference.tss)
    throw new Error("Pipeline reference TSS does not match the record.");
  if (new Set(a.designs.map((d) => d.id)).size !== a.designs.length)
    throw new Error("Duplicate design identifiers.");
  for (const d of a.designs) {
    if (
      d.target !== record.gene.symbol ||
      !a.supported_modes.includes(d.mode) ||
      !Number.isFinite(d.delta_target) ||
      !Number.isFinite(d.risk.score)
    )
      throw new Error("Invalid design metrics.");
    for (const r of [d.anchor_a, d.anchor_b, d.element]) {
      if (
        r.chrom !== record.gene.chrom ||
        !Number.isSafeInteger(r.start) ||
        !Number.isSafeInteger(r.end) ||
        r.end - r.start !== 5000 ||
        r.start < record.locus.tad.start ||
        r.end > record.locus.tad.end
      )
        throw new Error("Invalid 5 kb design interval.");
    }
    if (
      Math.abs(d.anchor_a.start - d.anchor_b.start) !==
      d.separation_kb * 1000
    )
      throw new Error("Anchor separation does not match coordinates.");
    if (
      d.sweep &&
      (d.sweep.length !== a.sweep_percentiles.length ||
        !d.sweep.every(Number.isFinite))
    )
      throw new Error("Invalid contact-strength sweep.");
    if (d.delta_at_p99 !== (d.sweep?.at(-1) ?? null))
      throw new Error("p99 does not match the supplied sweep.");
    if (d.mode === "down" && d.decoy_kb! * 1000 !== d.anchor_b.start)
      throw new Error("Decoy does not match anchor B.");
    if (
      d.status === "shortlisted" &&
      (d.rank === null ||
        d.risk.tier === "BLOCKED" ||
        d.confidence === "ABSTAIN" ||
        directedEffect(d) < a.source.config.MIN_EFFECT)
    )
      throw new Error("Rejected design in shortlist.");
  }
  if (shortlist(a).length !== a.source.n_meeting_dose_goal)
    throw new Error("Shortlist count mismatch.");
}

import { Vector3 } from "three";
import { RISE_NM, type LocusWindow } from "./dna-geometry.ts";

export type FoldedArm = {
  start: number;
  end: number;
  anchor: number;
  side: -1 | 1;
  role: "gene" | "element";
  display?: {
    start: number;
    end: number;
    anchor: number;
    windows: { start: number; end: number; offset: number }[];
  };
};
type ArmPose = { centers: Vector3[]; tangents: Vector3[] };
type FoldVariation = {
  phase: number;
  curl: number;
  frequency: number;
  lean: number;
  offsetX: number;
  depth: number;
};
export type FoldedDNA = {
  arms: FoldedArm[];
  omittedBp: number;
  shownBp: number;
  start: number;
  end: number;
  separation: number;
  hazeWidth: number;
  variation: FoldVariation;
  pose?: { progress: number; arms: ArmPose[] };
};
const FRAMES = 512;
const clamp = (x: number) => Math.max(0, Math.min(1, x));
const smooth = (x: number) => {
  const t = clamp(x);
  return t * t * (3 - 2 * t);
};
const ease = (t: number) => t * t * t * (t * (t * 6 - 15) + 10);

// Stable visual variation, seeded by source coordinates, not structural evidence.
// Returning to a hit restores its fold instead of generating a new random pose.
function foldVariation(locus: LocusWindow): FoldVariation {
  let seed = 2166136261;
  for (const char of `${locus.tss}:${locus.elementStart}:${locus.elementEnd}`)
    seed = Math.imul(seed ^ char.charCodeAt(0), 16777619) >>> 0;
  const random = () => {
    seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0;
    return seed / 4294967296;
  };
  const values = Array.from({ length: 6 }, random);
  return {
    phase: (values[5] - 0.5) * 3.6,
    curl: 0.75 + values[1] * 0.4,
    frequency: 0.2 + values[2] * 0.12,
    lean: (values[3] - 0.5) * 0.5,
    offsetX: -5 + values[4] * 5,
    depth: (values[0] - 0.5) * 0.2,
  };
}

// Explanatory composition, not an inferred TAD or molecular-dynamics trajectory.
// Genomic order and counts remain exact even when the enhancer is upstream.
export function foldedDNA(
  locus: LocusWindow,
  anchors?: { gene: number; element: number },
): FoldedDNA | null {
  if (locus.overlap) return null;
  if (anchors) {
    // Four actual local DNA windows: promoter, enhancer, and the two anchor-bin
    // centres. The gaps between windows have no molecular geometry. Display
    // offsets compress gaps only; one rung still represents one local bp.
    const interleaved = Math.max(Math.min(locus.tss, anchors.gene), Math.min(locus.midpoint, anchors.element))
      <= Math.min(Math.max(locus.tss, anchors.gene), Math.max(locus.midpoint, anchors.element));
    // A decoy can lie beyond the enhancer, so the TSS and enhancer belong to
    // the same genomic arm. Assign local windows by genomic proximity instead
    // of duplicating overlapping sequence on two arms. Existing promoter
    // staples retain their original two-arm composition.
    const sites = [locus.tss, locus.midpoint];
    const nearB = (p: number) => Math.abs(p - anchors.gene) <= Math.abs(p - anchors.element);
    const arms = (
      [
        ["gene", interleaved ? sites.filter(nearB) : [locus.tss], anchors.gene],
        ["element", interleaved ? sites.filter((p) => !nearB(p)) : [locus.midpoint], anchors.element],
      ] as const
    )
      .map(([role, features, anchor]) => {
        const spans = [
          ...features.map((feature) => ({
            start: Math.max(locus.start, Math.floor(feature) - 36),
            end: Math.min(locus.end, Math.floor(feature) + 36),
          })),
          {
            start: Math.max(locus.start, Math.floor(anchor) - 64),
            end: Math.min(locus.end, Math.floor(anchor) + 64),
          },
        ].sort((a, b) => a.start - b.start);
        const ranges: {start: number; end: number}[] = [];
        for (const span of spans) {
          const last = ranges.at(-1);
          if (last && span.start <= last.end) last.end = Math.max(last.end, span.end);
          else ranges.push({...span});
        }
        let offset = 0;
        const windows = ranges.map((r) => {
          const out = { ...r, offset };
          offset += r.end - r.start + 24;
          return out;
        });
        const w = windows.find((w) => anchor >= w.start && anchor <= w.end)!;
        return {
          start: ranges[0].start,
          end: ranges.at(-1)!.end,
          anchor,
          side: 1,
          role,
          display: {
            start: 0,
            end: offset - 24,
            anchor: w.offset + anchor - w.start,
            windows,
          },
        } as FoldedArm;
      })
      .sort((a, b) => a.start - b.start);
    if (arms[0].end >= arms[1].start) return null;
    arms[0].side = -1;
    arms[1].side = 1;
    const shownBp = arms.reduce(
      (sum, a) =>
        sum + a.display!.windows.reduce((n, w) => n + w.end - w.start, 0),
      0,
    );
    const omittedBp = arms[1].end - arms[0].start - shownBp;
    const emphasis = clamp((Math.log10(omittedBp) - 3) / 3.5);
    return {
      arms,
      shownBp,
      omittedBp,
      start: arms[0].start,
      end: arms[1].end,
      separation: 25 + emphasis * 5,
      hazeWidth: 8 + emphasis * 4,
      variation: foldVariation(locus),
    };
  }
  const low = Math.min(locus.tss, locus.midpoint),
    high = Math.max(locus.tss, locus.midpoint);
  const outer = Math.min(64, Math.floor((high - low) / 4));
  const inner = Math.min(96, Math.floor((high - low) / 4));
  if (outer < 8) return null;
  const arms = [low, high].map((anchor, index): FoldedArm => ({
    start: Math.max(
      locus.start,
      Math.floor(anchor) - (index === 0 ? outer : inner),
    ),
    end: Math.min(
      locus.end,
      Math.floor(anchor) + (index === 0 ? inner : outer),
    ),
    anchor,
    side: index === 0 ? -1 : 1,
    role: anchor === locus.tss ? "gene" : "element",
  }));
  const omittedBp = arms[1].start - arms[0].end;
  if (omittedBp <= 0) return null;
  // Logarithmic display emphasis, never a physical separation estimate.
  const emphasis = clamp((Math.log10(omittedBp) - 3) / 3.5);
  return {
    arms,
    omittedBp,
    shownBp: arms.reduce((n, a) => n + a.end - a.start, 0),
    start: arms[0].start,
    end: arms[1].end,
    separation: 17 + emphasis * 5,
    hazeWidth: 8 + emphasis * 4,
    variation: foldVariation(locus),
  };
}

function anchorPoint(arm: FoldedArm, p: number, display: FoldedDNA) {
  if (arm.role === "gene") return new Vector3(0, -8, 0);
  const e = ease(p);
  // A curved approach, easing to rest above the stationary promoter.
  return new Vector3(
    display.variation.offsetX * (1 - e) - 4.5 * Math.sin(Math.PI * e),
    -8 + display.separation * (1 - e) + 3.5 * e,
    0.7 * (1 - e) + 1.3 * Math.sin(Math.PI * e),
  );
}

// q is signed contour distance from the highlighted site, positive toward the
// omitted interval. Unit tangent integration bends the strand without stretching.
function tangent(
  q: number,
  arm: FoldedArm,
  p: number,
  variation: FoldVariation,
) {
  const e = ease(p),
    a = Math.abs(q),
    bend = smooth(a / 7);
  const wave = Math.sin(Math.PI * p) ** 2 * (1 - 0.45 * p);
  // More contour on the right makes a wider, gentler loop at the same bp scale.
  const rightQ = q / 1.5;
  let angle: number;
  if (arm.role === "gene") {
    angle = q >= 0 ? 1.25 * smooth(rightQ / 22) : 0.28 * Math.sin(q * 0.2);
  } else if (q >= 0) {
    const open =
      -0.02 -
      0.35 * smooth(rightQ / 22) +
      0.12 * Math.sin(q * 0.18 + variation.phase);
    const loop =
      1.3 * Math.sin((Math.PI * rightQ) / 21) - 1.0 * smooth(rightQ / 21);
    angle = bend * (open * (1 - e) + loop * e);
  } else {
    const folded =
      variation.curl *
        (1.15 * Math.sin(a * variation.frequency + variation.phase) +
          0.9 * Math.sin(a * 0.13 - variation.phase * 0.4)) +
      variation.lean;
    const settled = -0.22 - 0.22 * Math.sin(a * 0.19);
    angle = bend * (folded * (1 - e) + settled * e);
  }
  // Small travelling bends are deterministic, scrub smoothly and decay to rest.
  angle +=
    (arm.role === "gene" ? 0.045 : 0.19) *
    wave *
    bend *
    Math.sin(q * 0.43 - p * Math.PI * 5);
  const tilt =
    bend *
    (0.07 * Math.sin(q * 0.22) +
      (arm.role === "element"
        ? variation.depth * (1 - e) * Math.sin(q * 0.19 + variation.phase)
        : 0) +
      0.1 * wave * Math.sin(q * 0.31 - p * Math.PI * 4));
  return new Vector3(
    Math.cos(angle) * Math.cos(tilt),
    Math.sin(angle) * Math.cos(tilt),
    Math.sin(tilt),
  );
}

function pose(display: FoldedDNA, p: number) {
  if (display.pose?.progress === p) return display.pose.arms;
  const arms = display.arms.map((arm): ArmPose => {
    const orientation = -arm.side;
    const start = arm.display?.start ?? arm.start,
      end = arm.display?.end ?? arm.end,
      anchor = arm.display?.anchor ?? arm.anchor;
    const step = (end - start) / FRAMES;
    const qAt = (i: number) =>
      (start + step * i - anchor) * RISE_NM * orientation;
    const centers = new Array<Vector3>(FRAMES + 1);
    const tangents = Array.from({ length: FRAMES + 1 }, (_, i) =>
      tangent(qAt(i), arm, p, display.variation).multiplyScalar(orientation),
    );
    const anchorIndex = Math.floor((anchor - start) / step);
    // Integrate outward on both sides of the exact anchor, including half-base sites.
    for (const direction of [-1, 1]) {
      let previousQ = 0;
      let previous = anchorPoint(arm, p, display);
      for (
        let i = direction < 0 ? anchorIndex : anchorIndex + 1;
        i >= 0 && i <= FRAMES;
        i += direction
      ) {
        const q = qAt(i);
        previous = previous
          .clone()
          .addScaledVector(
            tangent((q + previousQ) / 2, arm, p, display.variation),
            q - previousQ,
          );
        centers[i] = previous;
        previousQ = q;
      }
    }
    return { centers, tangents };
  });
  display.pose = { progress: p, arms };
  return arms;
}

export function foldedFrame(
  position: number,
  progress: number,
  display: FoldedDNA,
) {
  const index = display.arms.findIndex(
    (a) => position >= a.start && position <= a.end,
  );
  if (index < 0)
    throw new Error(
      "No molecular geometry exists inside the omitted interval.",
    );
  const arm = display.arms[index],
    p = clamp(progress);
  const local = arm.display?.windows.find(
    (w) => position >= w.start && position <= w.end,
  );
  if (arm.display && !local)
    throw new Error(
      "No molecular geometry exists inside a compressed local gap.",
    );
  const mapped = local ? local.offset + position - local.start : position;
  const start = arm.display?.start ?? arm.start,
    end = arm.display?.end ?? arm.end;
  const frame = clamp((mapped - start) / (end - start)) * FRAMES;
  const lo = Math.min(FRAMES - 1, Math.floor(frame)),
    f = frame - lo;
  const current = pose(display, p)[index];
  const center =
    position === arm.anchor
      ? anchorPoint(arm, p, display)
      : current.centers[lo].clone().lerp(current.centers[lo + 1], f);
  const direction = current.tangents[lo]
    .clone()
    .lerp(current.tangents[lo + 1], f)
    .normalize();
  // Project a consistent reference normal to avoid helix flips as the curve bends.
  const normal = new Vector3(0, 0, 1)
    .addScaledVector(direction, -direction.z)
    .normalize();
  const binormal = new Vector3().crossVectors(direction, normal).normalize();
  const inner = arm.side < 0 ? arm.end - position : position - arm.start;
  const outer = arm.side < 0 ? position - arm.start : arm.end - position;
  const fade = local
    ? Math.min(
        smooth((position - local.start) / 14),
        smooth((local.end - position) / 14),
      )
    : Math.min(smooth(inner / 34), smooth(outer / 14));
  return { center, normal, binormal, fade };
}
export function foldedAnchor(
  position: number,
  progress: number,
  display: FoldedDNA,
) {
  return foldedFrame(position, progress, display).center;
}
export function hazeCenter(display: FoldedDNA, progress: number) {
  const a = foldedAnchor(display.arms[0].end, progress, display);
  const b = foldedAnchor(display.arms[1].start, progress, display);
  return a
    .add(b)
    .multiplyScalar(0.5)
    .add(new Vector3(3, 0, 0));
}

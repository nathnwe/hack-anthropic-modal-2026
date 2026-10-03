import { Vector3 } from "three";
import { RISE_NM, type LocusWindow } from "./dna-geometry.ts";

export type FoldedArm = {
  start: number;
  end: number;
  anchor: number;
  side: -1 | 1;
  role: "gene" | "element";
};
type ArmPose = { centers: Vector3[]; tangents: Vector3[] };
export type FoldedDNA = {
  arms: FoldedArm[];
  omittedBp: number;
  shownBp: number;
  start: number;
  end: number;
  separation: number;
  hazeWidth: number;
  pose?: { progress: number; arms: ArmPose[] };
};
const FRAMES = 512;
const clamp = (x: number) => Math.max(0, Math.min(1, x));
const smooth = (x: number) => {
  const t = clamp(x);
  return t * t * (3 - 2 * t);
};
const ease = (t: number) => t * t * t * (t * (t * 6 - 15) + 10);

// Explanatory composition, not an inferred TAD or molecular-dynamics trajectory.
// Genomic order and counts remain exact even when the enhancer is upstream.
export function foldedDNA(locus: LocusWindow): FoldedDNA | null {
  if (locus.overlap) return null;
  const low = Math.min(locus.tss, locus.midpoint),
    high = Math.max(locus.tss, locus.midpoint);
  const flank = Math.min(64, Math.floor((high - low) / 4));
  if (flank < 8) return null;
  const arms = [low, high].map((anchor, index): FoldedArm => ({
    start: Math.max(locus.start, Math.floor(anchor) - flank),
    end: Math.min(locus.end, Math.floor(anchor) + flank),
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
  };
}

function anchorPoint(arm: FoldedArm, p: number, display: FoldedDNA) {
  if (arm.role === "gene") return new Vector3(0, -8, 0);
  const e = ease(p);
  // A curved approach, easing to rest above the stationary promoter.
  return new Vector3(
    -3 * (1 - e) - 4.5 * Math.sin(Math.PI * e),
    -8 + display.separation * (1 - e) + 3.5 * e,
    0.7 * (1 - e) + 1.3 * Math.sin(Math.PI * e),
  );
}

// q is signed contour distance from the highlighted site, positive toward the
// omitted interval. Unit tangent integration bends the strand without stretching.
function tangent(q: number, arm: FoldedArm, p: number) {
  const e = ease(p),
    a = Math.abs(q),
    bend = smooth(a / 7);
  const wave = Math.sin(Math.PI * p) ** 2 * (1 - 0.45 * p);
  let angle: number;
  if (arm.role === "gene") {
    angle = q >= 0 ? 1.25 * smooth(q / 22) : 0.28 * Math.sin(q * 0.2);
  } else if (q >= 0) {
    const open = -0.3 - 0.9 * smooth(q / 22);
    const loop = 1.8 * Math.sin((Math.PI * q) / 21) - 1.2 * smooth(q / 21);
    angle = bend * (open * (1 - e) + loop * e);
  } else {
    const folded = 1.15 * Math.sin(a * 0.29) + 0.9 * Math.sin(a * 0.13);
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
    const step = (arm.end - arm.start) / FRAMES;
    const qAt = (i: number) =>
      (arm.start + step * i - arm.anchor) * RISE_NM * orientation;
    const centers = new Array<Vector3>(FRAMES + 1);
    const tangents = Array.from({ length: FRAMES + 1 }, (_, i) =>
      tangent(qAt(i), arm, p).multiplyScalar(orientation),
    );
    const anchorIndex = Math.floor((arm.anchor - arm.start) / step);
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
          .addScaledVector(tangent((q + previousQ) / 2, arm, p), q - previousQ);
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
  const frame = clamp((position - arm.start) / (arm.end - arm.start)) * FRAMES;
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
  const fade = Math.min(smooth(inner / 34), smooth(outer / 14));
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

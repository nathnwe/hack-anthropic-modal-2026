import { CatmullRomCurve3, Vector3 } from "three";
import { RISE_NM, type LocusWindow } from "./dna-geometry.ts";

export type FoldedArm = {
  start: number;
  end: number;
  anchor: number;
  side: -1 | 1;
  curve: CatmullRomCurve3;
  origin: Vector3;
  normals: Vector3[];
  binormals: Vector3[];
};
export type FoldedDNA = {
  arms: FoldedArm[];
  omittedBp: number;
  shownBp: number;
  start: number;
  end: number;
  separation: number;
  hazeWidth: number;
};
const FRAMES = 512;
const POINTS = [
  [-10, 7, -2],
  [-16, 4, 0],
  [-16, -3, 1],
  [-11, -6, 3],
  [-7, -3, 1],
  [-10, 2, 0],
  [-7, 7, -2],
  [-2, 9, -1],
  [1, 6, 1],
  [0, 3, 0],
];

// Local, idealised folds for explanatory viewing; no TAD is inferred. Source
// intervals stay exact. Only the unshown genomic interval is spatially compressed.
export function foldedDNA(locus: LocusWindow): FoldedDNA | null {
  if (locus.overlap) return null;
  const low = Math.min(locus.tss, locus.midpoint),
    high = Math.max(locus.tss, locus.midpoint);
  const flank = Math.min(64, Math.floor((high - low) / 4));
  if (flank < 8) return null;
  const windows = [low, high].map((anchor) => ({
    start: Math.max(locus.start, Math.floor(anchor) - flank),
    end: Math.min(locus.end, Math.floor(anchor) + flank),
    anchor,
  }));
  const omittedBp = windows[1].start - windows[0].end;
  if (omittedBp <= 0) return null;
  const arms = windows.map((w, index): FoldedArm => {
    const side = index === 0 ? -1 : 1;
    let points = POINTS.map(([x, y, z]) => new Vector3(x, y, z));
    if (side === 1)
      points = points
        .reverse()
        .map((p) => new Vector3(-p.x, p.y * 0.92, -p.z + Math.sin(p.y * 0.3)));
    const curve = new CatmullRomCurve3(points, false, "centripetal");
    curve.arcLengthDivisions = 2048;
    const scale = ((w.end - w.start) * RISE_NM) / curve.getLength();
    curve.points.forEach((p) => p.multiplyScalar(scale));
    curve.updateArcLengths();
    const frames = curve.computeFrenetFrames(FRAMES, false);
    return {
      ...w,
      side,
      curve,
      origin: curve.getPointAt((w.anchor - w.start) / (w.end - w.start)),
      normals: frames.normals,
      binormals: frames.binormals,
    };
  });
  // Logarithmic visual emphasis, NOT a physical scale or evidence of compaction.
  const emphasis = Math.max(0, Math.min(1, (Math.log10(omittedBp) - 3) / 3.5));
  return {
    arms,
    omittedBp,
    shownBp: arms.reduce((n, a) => n + a.end - a.start, 0),
    start: arms[0].start,
    end: arms[1].end,
    separation: 16 + emphasis * 8,
    hazeWidth: 8 + emphasis * 5,
  };
}

export function foldedFrame(
  position: number,
  progress: number,
  display: FoldedDNA,
) {
  const arm = display.arms.find(
    (a) => position >= a.start && position <= a.end,
  );
  if (!arm)
    throw new Error(
      "No molecular geometry exists inside the omitted interval.",
    );
  const t = Math.max(
    0,
    Math.min(1, (position - arm.start) / (arm.end - arm.start)),
  );
  const p = Math.max(0, Math.min(1, progress));
  const center = arm.curve.getPointAt(t).sub(arm.origin);
  const frame = t * FRAMES,
    lo = Math.min(FRAMES - 1, Math.floor(frame)),
    f = frame - lo;
  const tangent = arm.curve.getTangentAt(t);
  const normal = arm.normals[lo].clone().lerp(arm.normals[lo + 1], f);
  normal.addScaledVector(tangent, -normal.dot(tangent)).normalize();
  const binormal = new Vector3().crossVectors(tangent, normal).normalize();
  const rotation = -arm.side * 0.16 * p;
  const axis = new Vector3(0, 0, 1);
  for (const v of [center, normal, binormal]) v.applyAxisAngle(axis, rotation);
  // Rigid local motion preserves the helix; the haze carries the omitted span.
  center.add(
    new Vector3(
      (arm.side * (display.separation * (1 - p) + 3.4 * p)) / 2,
      -4,
      arm.side * 0.35,
    ),
  );
  const intoGap = arm.side < 0 ? arm.end - position : position - arm.start;
  const fadeT = Math.max(0, Math.min(1, intoGap / 34));
  const fade = fadeT * fadeT * (3 - 2 * fadeT);
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
    .add(new Vector3(0, 1, 0));
}

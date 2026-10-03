import { Vector3 } from "three";
import type { CRE, GeneRecord } from "./records.ts";

export const DNA_LENGTH = 50;
export type LocusWindow = {
  start: number;
  end: number;
  promoter: number;
  element: number;
  elementStart: number;
  elementEnd: number;
};

// Genomic coordinates choose the highlighted locations, never the 3D fold.
export function locusWindow(record: GeneRecord, cre: CRE): LocusWindow | null {
  const { tad } = record.locus;
  const tss = record.gene.tss;
  const midpoint = (cre.start + cre.end) / 2;
  if (
    ![tad.start, tad.end, tss, cre.start, cre.end].every(Number.isFinite) ||
    tad.end <= tad.start ||
    cre.end <= cre.start ||
    tss < tad.start ||
    tss > tad.end ||
    cre.start < tad.start ||
    cre.end > tad.end
  )
    return null;
  const separation = Math.abs(midpoint - tss);
  // Coincident anchors cannot demonstrate bringing two distinct regions together.
  if (separation < 1 || (tss >= cre.start && tss <= cre.end)) return null;
  const padding = Math.max(separation * 0.25, cre.end - cre.start);
  const start = Math.max(tad.start, Math.min(tss, cre.start) - padding);
  const end = Math.min(tad.end, Math.max(tss, cre.end) + padding);
  const at = (n: number) => (n - start) / (end - start);
  return {
    start,
    end,
    promoter: at(tss),
    element: at(midpoint),
    elementStart: at(cre.start),
    elementEnd: at(cre.end),
  };
}

// A designed continuous curve, in arbitrary display units. Both endpoints of
// the highlighted interval approach each other; out-of-plane tails avoid a
// flat crossing. This is explanatory geometry, not molecular dynamics.
export function dnaCenter(
  t: number,
  progress: number,
  locus: LocusWindow,
): Vector3 {
  const a = Math.min(locus.promoter, locus.element);
  const b = Math.max(locus.promoter, locus.element);
  const length = (b - a) * DNA_LENGTH;
  const p = Math.max(0, Math.min(1, progress));
  const angle = 1.1 + 4.52 * p;
  const radius = length / angle;
  const s = Math.max(0, Math.min(1, (t - a) / (b - a)));
  const theta = (s - 0.5) * angle;
  const point = new Vector3(
    radius * Math.sin(theta),
    radius * (Math.cos(theta) - Math.cos(angle / 2)) - 3,
    4 * p * Math.sin(s * Math.PI * 2),
  );
  if (t < a || t > b) {
    const tangent = new Vector3(
      Math.cos(theta),
      -Math.sin(theta),
      ((8 * Math.PI * p) / length) * Math.cos(s * Math.PI * 2),
    ).normalize();
    point.addScaledVector(tangent, (t < a ? t - a : t - b) * DNA_LENGTH);
  }
  return point;
}

export function helixFrames(
  locus: LocusWindow,
  progress: number,
  segments: number,
) {
  const centers: Vector3[] = [],
    normals: Vector3[] = [],
    binormals: Vector3[] = [];
  let previousNormal = new Vector3(0, 0, 1);
  for (let i = 0; i <= segments; i++) {
    const t = i / segments;
    const center = dnaCenter(t, progress, locus);
    const tangent = dnaCenter(Math.min(1, t + 0.0001), progress, locus)
      .sub(dnaCenter(Math.max(0, t - 0.0001), progress, locus))
      .normalize();
    // Parallel transport keeps helix phase continuous as the loop closes.
    const normal = previousNormal
      .clone()
      .addScaledVector(tangent, -previousNormal.dot(tangent))
      .normalize();
    const binormal = new Vector3().crossVectors(tangent, normal).normalize();
    centers.push(center);
    normals.push(normal);
    binormals.push(binormal);
    previousNormal = normal;
  }
  return { centers, normals, binormals };
}

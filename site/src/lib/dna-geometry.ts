import { Vector3 } from "three";
import type { CRE, GeneRecord } from "./records.ts";

// Idealised B-DNA dimensions. World units are nanometres; contour distance,
// diameter and twist use the same scale at every zoom level. No sequence skipped.
export const RISE_NM = 0.34;
export const HELIX_RADIUS_NM = 1;
export const BP_PER_TURN = 10.5;
export type LocusWindow = {
  start: number;
  end: number;
  tss: number;
  midpoint: number;
  elementStart: number;
  elementEnd: number;
  overlap: boolean;
};
export function locusWindow(record: GeneRecord, cre: CRE): LocusWindow | null {
  const { tad } = record.locus;
  const tss = record.gene.tss,
    midpoint = (cre.start + cre.end) / 2;
  if (
    ![tad.start, tad.end, tss, cre.start, cre.end].every(
      Number.isSafeInteger,
    ) ||
    record.gene.chrom.replace(/^chr/i, "") !== tad.chrom.replace(/^chr/i, "") ||
    tad.start < 0 ||
    tad.end <= tad.start ||
    cre.end <= cre.start ||
    tss < tad.start ||
    tss >= tad.end ||
    cre.start < tad.start ||
    cre.end > tad.end
  )
    return null;
  const separation = Math.abs(midpoint - tss);
  const padding = Math.max(80, Math.ceil(separation * 0.08));
  return {
    start: Math.max(tad.start, Math.min(tss, cre.start) - padding),
    end: Math.min(tad.end, Math.max(tss, cre.end) + padding),
    tss,
    midpoint,
    elementStart: cre.start,
    elementEnd: cre.end,
    overlap: tss >= cre.start && tss < cre.end,
  };
}
type Curve = {
  low: number;
  high: number;
  length: number;
  scale: number;
  total: number;
  gap: number;
  bulge: number;
  height: number;
  arc: Float64Array;
  xs: Float64Array;
  ys: Float64Array;
};
const SEGMENTS = 8192;
const parameters = new WeakMap<
  LocusWindow,
  { progress: number; value: Curve }
>();
const closedGaps = new WeakMap<LocusWindow, number>();
function point(u: number, gap: number, bulge: number, height: number) {
  const sine = Math.sin(Math.PI * u);
  return [(2 * u - 1) * (gap / 2 + bulge * sine * sine), height * sine];
}
function table(gap: number, bulge: number, height: number) {
  const arc = new Float64Array(SEGMENTS + 1),
    xs = new Float64Array(SEGMENTS + 1),
    ys = new Float64Array(SEGMENTS + 1);
  for (let i = 0; i <= SEGMENTS; i++) {
    [xs[i], ys[i]] = point(i / SEGMENTS, gap, bulge, height);
    if (i)
      arc[i] = arc[i - 1] + Math.hypot(xs[i] - xs[i - 1], ys[i] - ys[i - 1]);
  }
  return { arc, xs, ys, total: arc[SEGMENTS] };
}
export function loopParameters(locus: LocusWindow, progress: number): Curve {
  const cached = parameters.get(locus);
  if (cached?.progress === progress) return cached.value;
  const low = Math.min(locus.tss, locus.midpoint),
    high = Math.max(locus.tss, locus.midpoint),
    length = Math.max(RISE_NM, (high - low) * RISE_NM);
  let closed = closedGaps.get(locus);
  if (closed === undefined) {
    closed = 3 / length;
    for (let i = 0; i < 6; i++)
      closed = (3 * table(closed, 0.65, 0.7).total) / length;
    closedGaps.set(locus, closed);
  }
  const p = locus.overlap ? 0 : Math.max(0, Math.min(1, progress));
  const gap = 1 + (Math.min(1, closed) - 1) * p,
    bulge = 0.65 * p,
    height = 0.12 + 0.58 * p;
  const data = table(gap, bulge, height);
  const value = {
    low,
    high,
    length,
    scale: length / data.total,
    gap,
    bulge,
    height,
    ...data,
  };
  parameters.set(locus, { progress, value });
  return value;
}
export function dnaFrame(
  position: number,
  progress: number,
  locus: LocusWindow,
) {
  const c = loopParameters(locus, progress);
  if (c.high - c.low < 1)
    return {
      center: new Vector3((position - locus.tss) * RISE_NM, 0, 0),
      normal: new Vector3(0, 0, 1),
      binormal: new Vector3(0, -1, 0),
    };
  const distance = Math.max(
    0,
    Math.min(c.total, ((position - c.low) / (c.high - c.low)) * c.total),
  );
  let lo = 0,
    hi = SEGMENTS;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (c.arc[mid] < distance) lo = mid;
    else hi = mid;
  }
  const fraction = (distance - c.arc[lo]) / (c.arc[hi] - c.arc[lo]);
  const u = (lo + fraction) / SEGMENTS;
  // Interpolate on the same length table used for bp positions. This avoids
  // sub-bp stretching from interpolating the nonlinear parameter instead.
  const x = c.xs[lo] + fraction * (c.xs[hi] - c.xs[lo]);
  const y = c.ys[lo] + fraction * (c.ys[hi] - c.ys[lo]);
  const sine = Math.sin(Math.PI * u),
    cosine = Math.cos(Math.PI * u);
  const dx =
    c.gap +
    2 * c.bulge * sine * sine +
    (2 * u - 1) * 2 * c.bulge * sine * cosine * Math.PI;
  const dy = c.height * Math.PI * cosine;
  const tangent = new Vector3(dx, dy, 0).normalize();
  const center = new Vector3(x * c.scale, y * c.scale, 0);
  if (position < c.low)
    center.addScaledVector(tangent, (position - c.low) * RISE_NM);
  if (position > c.high)
    center.addScaledVector(tangent, (position - c.high) * RISE_NM);
  const normal = new Vector3(0, 0, 1);
  return {
    center,
    normal,
    binormal: new Vector3().crossVectors(tangent, normal),
  };
}
export const dnaCenter = (
  position: number,
  progress: number,
  locus: LocusWindow,
) => dnaFrame(position, progress, locus).center;

// Search the whole continuous path for portions inside the molecular viewport.
// Both arms of a closed loop remain present. Coarse geometry is retained outside it.
export function detailIntervals(
  locus: LocusWindow,
  progress: number,
  target: Vector3,
  radiusNm: number,
) {
  const c = loopParameters(locus, progress);
  const positions = [c.low, c.high];
  for (let i = 0; i < SEGMENTS; i++) {
    const ax = c.xs[i] * c.scale,
      ay = c.ys[i] * c.scale,
      dx = (c.xs[i + 1] - c.xs[i]) * c.scale,
      dy = (c.ys[i + 1] - c.ys[i]) * c.scale;
    const f = Math.max(
      0,
      Math.min(
        1,
        ((target.x - ax) * dx + (target.y - ay) * dy) / (dx * dx + dy * dy),
      ),
    );
    if (
      Math.hypot(ax + f * dx - target.x, ay + f * dy - target.y, target.z) <=
      radiusNm
    )
      positions.push(
        c.low +
          ((c.arc[i] + f * (c.arc[i + 1] - c.arc[i])) / c.total) *
            (c.high - c.low),
      );
  }
  for (const [anchor, side] of [
    [c.low, -1],
    [c.high, 1],
  ]) {
    const frame = dnaFrame(anchor, progress, locus);
    const tangent = new Vector3().crossVectors(frame.normal, frame.binormal);
    const projection = target.clone().sub(frame.center).dot(tangent) / RISE_NM;
    positions.push(
      anchor + (side < 0 ? Math.min(0, projection) : Math.max(0, projection)),
    );
  }
  const ranges = positions
    .filter(
      (p) =>
        p >= locus.start &&
        p <= locus.end &&
        dnaCenter(p, progress, locus).distanceTo(target) <= radiusNm * 1.2,
    )
    .map((p) => ({
      start: Math.max(locus.start, Math.floor(p - radiusNm / RISE_NM)),
      end: Math.min(locus.end, Math.ceil(p + radiusNm / RISE_NM)),
    }))
    .sort((a, b) => a.start - b.start);
  const merged: { start: number; end: number }[] = [];
  for (const r of ranges) {
    const last = merged.at(-1);
    if (last && r.start <= last.end) last.end = Math.max(last.end, r.end);
    else merged.push({ ...r });
  }
  return merged;
}

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
  tss: number;
  midpoint: number;
};

// Genomic coordinates choose the highlighted locations, never the 3D fold.
export function locusWindow(record: GeneRecord, cre: CRE): LocusWindow | null {
  const { tad } = record.locus;
  const tss = record.gene.tss;
  const midpoint = cre.start + (cre.end - cre.start) / 2;
  if (
    ![tad.start, tad.end, tss, cre.start, cre.end].every(
      Number.isSafeInteger,
    ) ||
    record.gene.chrom.replace(/^chr/i, "").toUpperCase() !==
      tad.chrom.replace(/^chr/i, "").toUpperCase() ||
    tad.start < 0 ||
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
    tss,
    midpoint,
  };
}

export type HelixSegment = {
  start: number;
  end: number;
  displayStart: number;
};
export type HelixDisplay = {
  segments: HelixSegment[];
  span: number;
  omittedBp: number;
  separationBp: number;
  shownBp: number;
  shape: LocusWindow;
  gap: { start: number; end: number } | null;
};

// A viewing choice, not a biological cutoff. Keep one rung per displayed bp;
// distant stretches are explicitly omitted instead of rescaling their bases.
const FLANK_BP = 36;
const MAX_CONTINUOUS_BP = 180;
const DISPLAY_GAP = 48;

export function helixDisplay(locus: LocusWindow): HelixDisplay {
  const low = Math.min(locus.tss, locus.midpoint);
  const high = Math.max(locus.tss, locus.midpoint);
  const start = Math.max(Math.ceil(locus.start), Math.floor(low) - FLANK_BP);
  const end = Math.min(Math.floor(locus.end), Math.ceil(high) + FLANK_BP);
  const segments: HelixSegment[] = [{ start, end, displayStart: 0 }];
  let omittedBp = 0;
  let gap: HelixDisplay["gap"] = null;
  if (end - start > MAX_CONTINUOUS_BP) {
    const leftEnd = Math.floor(low) + FLANK_BP;
    const rightStart = Math.ceil(high) - FLANK_BP;
    omittedBp = rightStart - leftEnd;
    segments[0].end = leftEnd;
    const gapStart = leftEnd - start;
    segments.push({
      start: rightStart,
      end,
      displayStart: gapStart + DISPLAY_GAP,
    });
    gap = { start: gapStart, end: gapStart + DISPLAY_GAP };
  }
  const last = segments.at(-1)!;
  const span = last.displayStart + last.end - last.start;
  const at = (position: number) => {
    const segment = segments.find(
      (s) => position >= s.start && position <= s.end,
    )!;
    return (segment.displayStart + position - segment.start) / span;
  };
  const shape = {
    ...locus,
    promoter: at(locus.tss),
    element: at(locus.midpoint),
  };
  return {
    segments,
    span,
    omittedBp,
    separationBp: high - low,
    shownBp: segments.reduce((sum, s) => sum + s.end - s.start, 0),
    shape,
    gap: gap ? { start: gap.start / span, end: gap.end / span } : null,
  };
}

export function displayCenter(
  t: number,
  progress: number,
  display: HelixDisplay,
): Vector3 {
  const p = Math.max(0, Math.min(1, progress));
  if (!display.gap) return dnaCenter(t, p, display.shape);
  const a = Math.min(display.shape.promoter, display.shape.element);
  const b = Math.max(display.shape.promoter, display.shape.element);
  const separation = 14 - 10.6 * p;
  const local = (position: number, side: number) => {
    const offset = (position - (side < 0 ? a : b)) * display.span;
    return new Vector3(
      side * (separation / 2 + 0.0018 * offset * offset),
      -side * offset * 0.3,
      side * 0.6,
    );
  };
  if (t <= display.gap.start) return local(t, -1);
  if (t >= display.gap.end) return local(t, 1);
  // Only the dotted omission guide follows this arch; no molecular bonds cross it.
  const left = local(display.gap.start, -1);
  const right = local(display.gap.end, 1);
  const u = (t - display.gap.start) / (display.gap.end - display.gap.start);
  return left.lerp(right, u).add(new Vector3(0, 5 * Math.sin(Math.PI * u), 0));
}

export function displayFrames(
  display: HelixDisplay,
  progress: number,
  segments: number,
) {
  return transportFrames((t) => displayCenter(t, progress, display), segments);
}

export function helixSamples(display: HelixDisplay) {
  return display.segments.map((segment, index) => {
    const pairs = segment.end - segment.start;
    return Array.from({ length: pairs * 2 - 1 }, (_, i) => ({
      position: segment.start + 0.5 + i / 2,
      t: (segment.displayStart + 0.5 + i / 2) / display.span,
      rung: i % 2 === 0,
      // Fade only the ends bordering omitted sequence. Never join them by a bond.
      fade: display.gap
        ? Math.min(1, (index === 0 ? pairs - 0.5 - i / 2 : 0.5 + i / 2) / 7)
        : 1,
    }));
  });
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
  return transportFrames((t) => dnaCenter(t, progress, locus), segments);
}

function transportFrames(centerAt: (t: number) => Vector3, segments: number) {
  const centers: Vector3[] = [],
    normals: Vector3[] = [],
    binormals: Vector3[] = [];
  let previousNormal = new Vector3(0, 0, 1);
  for (let i = 0; i <= segments; i++) {
    const t = i / segments;
    const center = centerAt(t);
    const tangent = centerAt(Math.min(1, t + 0.0001))
      .sub(centerAt(Math.max(0, t - 0.0001)))
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

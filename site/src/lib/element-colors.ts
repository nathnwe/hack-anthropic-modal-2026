import { contactStrength, type CRE, type GeneRecord } from "./records.ts";
export const SCORE_LOW = "#f0dea0";
export const SCORE_HIGH = "#ad5641";
export const SCORE_MISSING = "#aaa99b";
export function rankingValue(
  record: GeneRecord,
  cre: CRE,
  criterion: string,
): number | null {
  if (criterion === "re2g")
    return Number.isFinite(cre.score_rE2G) ? cre.score_rE2G : null;
  if (criterion === "contact") return contactStrength(record, cre);
  if (criterion === "distance")
    return Math.abs((cre.start + cre.end) / 2 - record.gene.tss);
  return null;
}
export function rankingColors(record: GeneRecord, criterion: string) {
  const values = record.locus.cres
    .map((c) => rankingValue(record, c, criterion))
    .filter((v): v is number => v !== null && Number.isFinite(v));
  const min = values.length ? Math.min(...values) : null,
    max = values.length ? Math.max(...values) : null;
  const color = (cre: CRE) => {
    const value = rankingValue(record, cre, criterion);
    if (value === null || min === null || max === null) return SCORE_MISSING;
    let t = max === min ? 0.5 : (value - min) / (max - min);
    if (criterion === "distance") t = 1 - t;
    const channels = [0, 1, 2].map((i) => {
      const low = parseInt(SCORE_LOW.slice(1 + 2 * i, 3 + 2 * i), 16),
        high = parseInt(SCORE_HIGH.slice(1 + 2 * i, 3 + 2 * i), 16);
      return Math.round(low + (high - low) * t)
        .toString(16)
        .padStart(2, "0");
    });
    return "#" + channels.join("");
  };
  return { min, max, color };
}

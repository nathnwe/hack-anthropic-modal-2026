import type { Mode } from "./records.ts";

export type RecordEntry = {
  symbol: string;
  file: string;
  direction: string | null;
  analysis?: boolean;
  modes?: Mode[];
};

// Keep explicit legacy samples separate from the new pipeline analyses.
export function selectRecord(
  entries: RecordEntry[], symbol: string, mode: Mode | null, source: string | null,
) {
  const matches = entries.filter((entry) => entry.symbol === symbol);
  const legacy = matches.find((entry) => !entry.analysis);
  if (source === "sample" && legacy) return legacy;
  const analysis = matches.filter((entry) => entry.analysis);
  return analysis.find((entry) => !mode || entry.modes?.includes(mode)) ||
    (source === "analysis" ? analysis[0] : legacy || analysis[0]);
}

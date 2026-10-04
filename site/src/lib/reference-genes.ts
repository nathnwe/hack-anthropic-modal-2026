export type GeneModel = {
  chr: string;
  start: number;
  end: number;
  name: string;
  id: string;
  strand: string;
  cdStart: number;
  cdEnd: number;
  exons: { start: number; end: number }[];
};
export type ReferenceSnapshot = {
  assembly: string;
  windows: {
    chrom: string;
    start: number;
    end: number;
    features: GeneModel[];
  }[];
};
/** null means the requested interval is outside cached coverage; [] is a covered empty interval. */
export function referenceGenes(
  snapshot: ReferenceSnapshot,
  chr: string,
  start: number,
  end: number,
): GeneModel[] | null {
  if (snapshot.assembly !== "GRCh38") return null;
  const window = snapshot.windows.find(
    (w) => w.chrom === chr && w.start <= start && w.end >= end,
  );
  return window
    ? window.features.filter((f) => f.end > start && f.start < end)
    : null;
}

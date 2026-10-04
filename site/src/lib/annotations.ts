export type Feature = {
  feature_type: string;
  id: string;
  external_name?: string;
  Parent?: string;
  seq_region_name: string;
  assembly_name: string;
  start: number;
  end: number;
  strand: number;
  is_canonical?: number;
};
export function annotationFeatures(rows: Feature[]) {
  const valid = rows.filter(
    (r) =>
      r.assembly_name === "GRCh38" &&
      Number.isSafeInteger(r.start) &&
      Number.isSafeInteger(r.end) &&
      r.start > 0 &&
      r.end >= r.start,
  );
  return valid
    .filter((r) => r.feature_type === "transcript" && r.is_canonical === 1)
    .map((r) => {
      const exons = valid
        .filter((e) => e.feature_type === "exon" && e.Parent === r.id)
        .map((e) => ({ start: e.start - 1, end: e.end }))
        .sort((a, b) => a.start - b.start);
      const cds = valid.filter(
        (e) => e.feature_type === "cds" && e.Parent === r.id,
      );
      return {
        chr: "chr" + (r.seq_region_name === "MT" ? "M" : r.seq_region_name),
        start: r.start - 1,
        end: r.end,
        name: r.external_name || r.id,
        strand: r.strand === 1 ? "+" : "-",
        exons,
        cdStart: cds.length ? Math.min(...cds.map((e) => e.start - 1)) : r.end,
        cdEnd: cds.length ? Math.max(...cds.map((e) => e.end)) : r.end,
      };
    });
}

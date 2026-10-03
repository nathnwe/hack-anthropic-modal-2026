import type { GeneRecord } from "./records.ts";

export type K562Dataset = {
  kind: string;
  symbol: string;
  cell_type: string;
  assembly: string;
  coordinate_system: string;
  source_file: string;
  source_sha256: string;
  reference: NonNullable<GeneRecord["caseStudy"]>["reference"];
  links: {
    id: string;
    csv_row: number;
    chrom: string;
    start: number;
    end: number;
    score: number;
    model: string;
    evidence: string;
    accession: string;
    url: string;
    annotation: string;
  }[];
};

// A frontend-only adapter. The source file is a prediction dataset, not a
// pipeline contract record. In particular it has no TAD, contacts or staples.
export function k562Record(data: K562Dataset): GeneRecord {
  if (
    data.kind !== "re2g-case-study" ||
    data.assembly !== "GRCh38" ||
    data.symbol !== "MYC" ||
    !data.links.length ||
    !Number.isSafeInteger(data.reference.tss)
  )
    throw new Error("Invalid K562 case study.");
  for (const r of data.links)
    if (
      r.chrom !== "chr8" ||
      !Number.isSafeInteger(r.start) ||
      !Number.isSafeInteger(r.end) ||
      r.end <= r.start ||
      !Number.isFinite(r.score)
    )
      throw new Error("Invalid K562 interval.");
  const start = Math.min(
    data.reference.start,
    ...data.links.map((r) => r.start),
  );
  const end = Math.max(data.reference.end, ...data.links.map((r) => r.end));
  return {
    illustrative: false,
    caseStudy: {
      assembly: data.assembly,
      coordinate_system: data.coordinate_system,
      reference: data.reference,
      source_file: data.source_file,
    },
    gene: {
      symbol: "MYC",
      chrom: "chr8",
      tss: data.reference.tss,
      disease: "K562 · GRCh38 · ENCODE-rE2G predictions",
      direction: null,
      evidence: [
        "ENCODE-rE2G: https://www.nature.com/articles/s41586-026-10781-4",
        `Reference transcript ${data.reference.transcript}: ${data.reference.url}`,
      ],
    },
    locus: {
      cell_type: "K562",
      tad: { chrom: "chr8", start: Math.max(0, start - 1000), end: end + 1000 },
      cres: data.links.map((r) => ({
        id: r.id,
        start: r.start,
        end: r.end,
        type: "enhancer",
        score_rE2G: r.score,
        source: {
          csv_row: r.csv_row,
          accession: r.accession,
          url: r.url,
          model: r.model,
          annotation: r.annotation,
          evidence: r.evidence,
        },
      })),
      contacts: [],
    },
    staples: [],
    risks: [],
    ranking: [],
  };
}

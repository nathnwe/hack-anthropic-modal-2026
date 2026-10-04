import type { GeneRecord, ReferenceGene } from "./records.ts";
export const ENSEMBL = "https://rest.ensembl.org";
export type GeneSuggestion = { symbol: string; name: string };
export async function jsonRequest(url: string, signal?: AbortSignal) {
  const response = await fetch(url, {
    signal: signal
      ? AbortSignal.any([signal, AbortSignal.timeout(15000)])
      : AbortSignal.timeout(15000),
  });
  if (!response.ok)
    throw new Error(
      response.status === 429
        ? "The annotation service is busy. Please try again shortly."
        : `The annotation service returned HTTP ${response.status}.`,
    );
  return response.json();
}
const suggestions = new Map<string, GeneSuggestion[]>();
export async function suggestGenes(
  query: string,
  signal?: AbortSignal,
): Promise<GeneSuggestion[]> {
  const term = query.trim();
  if (term.length < 2) return [];
  if (suggestions.has(term.toUpperCase()))
    return suggestions.get(term.toUpperCase())!;
  const tokens = term
    .split(/\s+/)
    .map((t) => t.replace(/[^a-zA-Z0-9_.-]/g, "").replace(/-/g, "\\-"))
    .filter(Boolean);
  if (!tokens.length) return [];
  const q =
    tokens.length === 1
      ? `symbol:${tokens[0]}* OR alias:${tokens[0]}*`
      : `name:(${tokens.map((t) => t + "*").join(" AND ")})`;
  const params = new URLSearchParams({
    q,
    species: "human",
    fields: "symbol,name,taxid",
    size: "12",
  });
  const data = await jsonRequest(
    `https://mygene.info/v3/query?${params}`,
    signal,
  );
  const seen = new Set<string>();
  const hits: GeneSuggestion[] = (Array.isArray(data.hits) ? data.hits : [])
    .filter((hit: { symbol?: string; taxid?: number }) => {
      if (hit.taxid !== 9606 || !hit.symbol || seen.has(hit.symbol))
        return false;
      seen.add(hit.symbol);
      return true;
    })
    .map((hit: GeneSuggestion) => ({
      symbol: hit.symbol,
      name: hit.name || hit.symbol,
    }));
  hits.sort(
    (a, b) =>
      Number(b.symbol.toUpperCase() === term.toUpperCase()) -
      Number(a.symbol.toUpperCase() === term.toUpperCase()),
  );
  suggestions.set(term.toUpperCase(), hits.slice(0, 8));
  return suggestions.get(term.toUpperCase())!;
}
export type EnsemblGene = {
  assembly_name: string;
  species: string;
  seq_region_name: string;
  start: number;
  end: number;
  strand: number;
  id: string;
  display_name: string;
  description?: string;
  canonical_transcript?: string;
};
export function referenceFromEnsembl(data: EnsemblGene): ReferenceGene {
  if (
    data.assembly_name !== "GRCh38" ||
    data.species !== "homo_sapiens" ||
    !/^(?:[1-9]|1\d|2[0-2]|X|Y|MT)$/.test(data.seq_region_name) ||
    !Number.isSafeInteger(data.start) ||
    !Number.isSafeInteger(data.end) ||
    data.start < 1 ||
    data.end < data.start ||
    ![1, -1].includes(data.strand)
  )
    throw new Error(
      "This result does not have supported human GRCh38 coordinates.",
    );
  return {
    assembly: "GRCh38",
    chrom:
      data.seq_region_name === "MT" ? "chrM" : `chr${data.seq_region_name}`,
    start: data.start - 1,
    end: data.end,
    tss: data.strand === 1 ? data.start - 1 : data.end - 1,
    strand: data.strand === 1 ? "+" : "-",
    id: data.id,
    name: (data.description || data.display_name).replace(
      /\s*\[Source:.*$/,
      "",
    ),
    url: `https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=${encodeURIComponent(data.id)}`,
  };
}
export async function liveGeneRecord(symbol: string): Promise<GeneRecord> {
  if (!/^[a-zA-Z0-9_.-]{1,60}$/.test(symbol))
    throw new Error("Enter a gene symbol, or choose a search suggestion.");
  let data: EnsemblGene;
  try {
    data = await jsonRequest(
      `${ENSEMBL}/lookup/symbol/homo_sapiens/${encodeURIComponent(symbol)}?content-type=application/json`,
    );
  } catch (error) {
    throw new Error(
      `Could not resolve “${symbol}” in the live gene reference. Choose a suggestion or try again. ${error instanceof Error ? error.message : ""}`,
    );
  }
  const reference = referenceFromEnsembl(data);
  if (data.canonical_transcript) {
    const transcript = await jsonRequest(
      `${ENSEMBL}/lookup/id/${encodeURIComponent(data.canonical_transcript.split(".")[0])}?content-type=application/json`,
    );
    if (
      transcript.assembly_name === "GRCh38" &&
      transcript.seq_region_name === data.seq_region_name &&
      Number.isSafeInteger(transcript.start) &&
      Number.isSafeInteger(transcript.end) &&
      transcript.start >= data.start &&
      transcript.end <= data.end &&
      transcript.end >= transcript.start &&
      transcript.strand === data.strand
    )
      reference.tss =
        transcript.strand === 1 ? transcript.start - 1 : transcript.end - 1;
    else
      throw new Error(
        "The canonical transcript coordinates do not match this reference gene.",
      );
  }
  return {
    reference,
    illustrative: false,
    gene: {
      symbol: data.display_name || symbol.toUpperCase(),
      chrom: reference.chrom,
      tss: reference.tss,
      disease: reference.name,
      direction: null,
      evidence: [reference.url],
    },
    locus: {
      cell_type: "",
      tad: {
        chrom: reference.chrom,
        start: Math.max(0, reference.start - 10000),
        end: reference.end + 10000,
      },
      cres: [],
      contacts: [],
    },
    staples: [],
    risks: [],
    ranking: [],
  };
}
export function genomeReference(record: GeneRecord): ReferenceGene | null {
  if (record.illustrative) return null;
  if (record.reference) return record.reference;
  const c = record.caseStudy;
  if (!c || c.assembly !== "GRCh38") return null;
  return {
    assembly: "GRCh38",
    chrom: record.gene.chrom,
    start: c.reference.start,
    end: c.reference.end,
    tss: c.reference.tss,
    strand: c.reference.strand === "-" ? "-" : "+",
    id: c.reference.transcript,
    name: record.gene.symbol,
    url: c.reference.url,
  };
}

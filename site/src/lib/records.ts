export type Mode = "up" | "down" | "off";
export type Region = { chrom: string; start: number; end: number };
export type CRE = {
  id: string;
  start: number;
  end: number;
  type: "enhancer" | "silencer" | "insulator" | "promoter";
  score_rE2G: number;
};
export type Contact = { a: number; b: number; strength: number };
export type GeneRecord = {
  illustrative?: boolean;
  gene: {
    symbol: string;
    chrom: string;
    tss: number;
    disease: string;
    direction: "too_much" | "too_little";
    evidence: string[];
  };
  locus: { cell_type: string; tad: Region; cres: CRE[]; contacts: Contact[] };
  staples: {
    id: string;
    anchor_a: Region;
    anchor_b: Region;
    mode: Mode;
    predicted_delta: string;
  }[];
  risks: {
    staple_id: string;
    flags: {
      type: string;
      gene?: string;
      severity: "low" | "medium" | "high";
      detail: string;
    }[];
  }[];
  ranking: {
    staple_id: string;
    score: number;
    tier: "A" | "B" | "C";
    rationale: string[];
    method_recommendation: string;
  }[];
};

export const modeLabels: Record<Mode, string> = {
  up: "Upregulate",
  down: "Downregulate",
  off: "Silence",
};
export const isMode = (value: string | null): value is Mode =>
  value === "up" || value === "down" || value === "off";
export const formatNumber = (n: number) => n.toLocaleString("en-GB");
export const regionLabel = (r: Region) =>
  `${r.chrom.replace(/^chr/i, "chr").startsWith("chr") ? r.chrom : "chr" + r.chrom}:${formatNumber(r.start)}–${formatNumber(r.end)}`;
export const isPlaceholder = (text: string) =>
  /placeholder|illustrative/i.test(text);

// The contract supplies point contacts, without bin resolution. Only an exact
// TSS endpoint and an endpoint inside a CRE count; do not invent a tolerance.
// Endpoints are treated as inclusive for display; assembly/convention is unknown.
export function contactsFor(record: GeneRecord, cre: CRE): Contact[] {
  const inCRE = (position: number) =>
    position >= cre.start && position <= cre.end;
  return record.locus.contacts.filter(
    (c) =>
      (c.a === record.gene.tss && inCRE(c.b)) ||
      (c.b === record.gene.tss && inCRE(c.a)),
  );
}
export function contactStrength(record: GeneRecord, cre: CRE): number | null {
  const values = contactsFor(record, cre).map((c) => c.strength);
  return values.length ? Math.max(...values) : null;
}
export function sortElements(record: GeneRecord, criterion: string): CRE[] {
  const result = [...record.locus.cres];
  if (criterion === "contact")
    return result.sort(
      (a, b) =>
        (contactStrength(record, b) ?? -Infinity) -
          (contactStrength(record, a) ?? -Infinity) || a.id.localeCompare(b.id),
    );
  if (criterion === "re2g")
    return result.sort(
      (a, b) => b.score_rE2G - a.score_rE2G || a.id.localeCompare(b.id),
    );
  return result;
}
// A new ranking always starts the viewer on its first available element.
export function selectedElementFor(
  record: GeneRecord,
  criterion: string,
): CRE | null {
  return sortElements(record, criterion)[0] ?? null;
}
export function candidatesFor(record: GeneRecord, mode: Mode) {
  const score = (id: string) =>
    record.ranking.find((r) => r.staple_id === id)?.score ?? -Infinity;
  return record.staples
    .filter((s) => s.mode === mode)
    .sort((a, b) => score(b.id) - score(a.id) || a.id.localeCompare(b.id));
}
export function evidenceLink(text: string): string | null {
  if (isPlaceholder(text)) return null;
  const url = text
    .match(/https?:\/\/[^\s<>"']+/i)?.[0]
    ?.replace(/[.,;)]+$/, "");
  if (url) return url;
  const pmid = text.match(/\bPMID\s*:?\s*(\d+)\b/i)?.[1];
  return pmid ? `https://pubmed.ncbi.nlm.nih.gov/${pmid}/` : null;
}
export function escapeHTML(value: unknown): string {
  return String(value).replace(
    /[&<>"']/g,
    (c) =>
      ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[
        c
      ]!,
  );
}

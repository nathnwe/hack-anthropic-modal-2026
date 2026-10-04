import catalogue from "../data/k562-tracks.json" with { type: "json" };
import type { GeneRecord } from "./records.ts";

/** Assay data are specific to both the biosample and the coordinate assembly. */
export function contextualTracks(record: GeneRecord) {
  const assembly =
    record.analysis?.reference.assembly ??
    record.caseStudy?.assembly ??
    record.reference?.assembly;
  return !record.illustrative &&
    record.locus.cell_type === "K562" &&
    assembly === "GRCh38"
    ? catalogue.tracks
    : [];
}
export function signalTrackConfig(
  track: (typeof catalogue.tracks)[number],
  index: number,
) {
  return {
    id: track.id,
    name: track.name,
    type: "wig",
    format: "bigwig",
    url: track.url,
    color: track.color,
    height: 48,
    min: 0,
    autoscale: true,
    windowFunction: "mean",
    graphType: "bar",
    order: 10 + index,
    timeout: 30000,
    description: `K562 · GRCh38 · ${track.output_type} · pooled replicates ${track.biological_replicates.join(", ")} · ${track.accession}`,
  };
}
export { catalogue as k562TrackCatalogue };

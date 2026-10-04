// Export the pipeline's precomputed K562 runs without changing the shared schema.
// Refresh reference annotations explicitly; ordinary exports/builds are offline.
import { readFile, writeFile, mkdir } from "node:fs/promises";
import { createHash } from "node:crypto";
import { execFileSync } from "node:child_process";
import { fileURLToPath } from "node:url";

const root = new URL("../../", import.meta.url);
const snapshotURL = new URL("./reference-annotations.json", import.meta.url);
const sourceRoot = new URL("site/public/data/pipeline-source/", root);
const refIndex = process.argv.indexOf("--source-ref");
if (refIndex !== -1) {
  const ref = process.argv[refIndex + 1];
  if (!ref || ref.startsWith("-")) throw new Error("--source-ref requires a fetched Git ref");
  const git = (...args) => execFileSync("git", args, {cwd: fileURLToPath(root)});
  const revision = git("rev-parse", "--verify", `${ref}^{commit}`).toString().trim();
  const paths = git("ls-tree", "-r", "--name-only", revision, "data/derived").toString().trim().split("\n")
    .filter((p) => /^data\/derived\/(demo_(ldlr|nsd1|runx1|nf1|mecp2_up|mecp2_down|myc_down)|website_results|screen_results|screen\/(up|down)_[A-Z0-9]+)\.json$/.test(p));
  const manifest = {revision, repository: "https://github.com/nathnwe/hack-anthropic-modal-2026", files: {}};
  for (const path of paths) {
    const bytes = git("show", `${revision}:${path}`);
    const file = path.replace("data/derived/", "");
    const destination = new URL(file, sourceRoot);
    await mkdir(new URL(".", destination), {recursive:true});
    await writeFile(destination, bytes);
    manifest.files[path] = {file, sha256:createHash("sha256").update(bytes).digest("hex")};
  }
  await writeFile(new URL("manifest.json", sourceRoot), JSON.stringify(manifest,null,2)+"\n");
}
const sourceManifest = JSON.parse(await readFile(new URL("manifest.json", sourceRoot), "utf8"));
const targets = {
  ldlr: {
    disease: "Familial hypercholesterolaemia",
    pmids: ["33740630", "15321837"],
  },
  nsd1: {
    disease: "Sotos syndrome",
    pmids: ["11896389", "17565729", "23190751"],
  },
  runx1: {
    disease: "Familial platelet disorder with AML predisposition",
    pmids: ["31723833", "18478040", "19357396", "33616470"],
  },
  nf1: {
    disease: "Neurofibromatosis type 1",
    pmids: ["1757093", "1302608", "34427956", "34750850", "28230061"],
  },
  mecp2_up: { disease: "Rett syndrome (loss of function)", pmids: ["31206249", "10508514", "16647997"] },
  mecp2_down: { disease: "MECP2 duplication syndrome (gain)", pmids: ["16080119", "32043567", "29618507"] },
  myc_down: { disease: "Oncogene, overexpressed in cancer", pmids: ["27708057"] },
};
const runs = await Promise.all(
  Object.keys(targets).map(async (key) => {
    const text = await readFile(
      new URL(`demo_${key}.json`, sourceRoot),
      "utf8",
    );
    const hash = createHash("sha256").update(text).digest("hex");
    if (hash !== sourceManifest.files[`data/derived/demo_${key}.json`]?.sha256)
      throw new Error(`Source snapshot checksum mismatch: ${key}`);
    return { key, text, run: JSON.parse(text) };
  }),
);
let snapshots;
if (process.argv.includes("--refresh-annotations")) {
  snapshots = {};
  for (const { key, run } of runs) {
    const [chrom, range] = run.window.split(":");
    const [start, end] = range.split("-").map(Number);
    const url = `https://api.genome.ucsc.edu/getData/track?genome=hg38;track=ncbiRefSeqCurated;chrom=${chrom};start=${start};end=${end}`;
    const response = await fetch(url, { signal: AbortSignal.timeout(60000) });
    if (!response.ok) throw new Error(`UCSC ${key}: HTTP ${response.status}`);
    const json = await response.json();
    const transcripts = json.ncbiRefSeqCurated.filter(
      (t) => t.name2 === run.targets[0] && /^(NM_|NR_)/.test(t.name),
    );
    if (!transcripts.length)
      throw new Error(`No sourced transcripts for ${key}`);
    snapshots[run.targets[0].toLowerCase()] = { url, retrieved: new Date().toISOString(), transcripts };
    console.log(`Retrieved ${key}: ${transcripts.length} transcripts`);
  }
  await writeFile(snapshotURL, JSON.stringify(snapshots, null, 2) + "\n");
} else snapshots = JSON.parse(await readFile(snapshotURL, "utf8"));

for (const { key, text, run } of runs) {
  const symbol = run.targets[0], mode = run.direction,
    annotation = snapshots[symbol.toLowerCase()];
  if (!["up", "down"].includes(mode)) throw new Error("Unsupported source direction");
  const sign = mode === "down" ? -1 : 1;
  const ts = annotation.transcripts;
  // Match pipeline.real_data.fetch_genes: the TSS shared by the most isoforms,
  // first encountered on a tie. Preserve its minus-strand txEnd convention.
  const counts = new Map();
  for (const t of ts) {
    const p = t.strand === "+" ? t.txStart : t.txEnd;
    counts.set(p, (counts.get(p) || 0) + 1);
  }
  const tss = [...counts].sort((a, b) => b[1] - a[1])[0][0];
  const transcript = ts.find(
    (t) => (t.strand === "+" ? t.txStart : t.txEnd) === tss,
  );
  const [chrom, range] = run.window.split(":");
  const [start, end] = range.split("-").map(Number);
  const bin = run.provenance.bin_bp;
  const reference = {
    assembly: "GRCh38",
    chrom,
    start: Math.min(...ts.map((t) => t.txStart)),
    end: Math.max(...ts.map((t) => t.txEnd)),
    tss,
    strand: transcript.strand,
    id: transcript.name,
    name: symbol,
    url: annotation.url,
  };
  if (tss < start || tss >= end)
    throw new Error(`${symbol} TSS outside analysis window`);
  const region = (kb) => ({ chrom, start: kb * 1000, end: kb * 1000 + bin });
  const designs = [...run.ranked, ...run.best_below_floor].map((row) => {
    const id = `${key}-a${row.anchor_a_kb * 1000}-b${row.anchor_b_kb * 1000}`;
    const status =
      row.risk.tier === "BLOCKED"
        ? "blocked"
        : row.rank !== null
          ? "shortlisted"
          : "below_floor";
    return {
      ...row,
      mode,
      id,
      status,
      anchor_a: region(row.anchor_a_kb),
      anchor_b: region(row.anchor_b_kb),
      element: {
        id: `${key}-e${row.enhancer_kb * 1000}`,
        ...region(row.enhancer_kb),
        type: "enhancer",
      },
    };
  });
  if (new Set(designs.map((d) => d.id)).size !== designs.length)
    throw new Error("Duplicate design IDs");
  for (const d of designs) {
    for (const r of [d.element, d.anchor_a, d.anchor_b])
      if (r.start < start || r.end > end)
        throw new Error(`${d.id}: region outside window`);
    if (
      Math.abs(d.anchor_a.start - d.anchor_b.start) !==
      d.separation_kb * 1000
    )
      throw new Error("Separation mismatch");
    if (
      d.status === "shortlisted" &&
      (sign * d.delta_target < run.config.MIN_EFFECT || d.risk.tier === "BLOCKED")
    )
      throw new Error("Invalid shortlist");
  }
  const accepted = designs.filter((d) => d.status === "shortlisted");
  const record = {
    illustrative: false,
    gene: {
      symbol,
      chrom,
      tss,
      disease: targets[key].disease,
      direction: mode === "down" ? "too_much" : "too_little",
      evidence: [
        `Pipeline results: https://github.com/nathnwe/hack-anthropic-modal-2026/blob/${sourceManifest.revision}/docs/results-report.md`,
        "ClinGen dosage sensitivity: https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv",
        ...targets[key].pmids.map((p) => `PMID: ${p}`),
      ],
    },
    // Compatibility interval only: this run did not call TAD boundaries.
    locus: {
      cell_type: "K562",
      tad: { chrom, start, end },
      scope_kind: "analysis_window",
      cres: [],
      contacts: [],
    },
    staples: accepted.map((d) => ({
      id: d.id,
      anchor_a: d.anchor_a,
      anchor_b: d.anchor_b,
      mode,
      predicted_delta: `${mode === "down" ? "decrease" : "increase"}_at_reference_strength`,
    })),
    risks: accepted.map((d) => ({
      staple_id: d.id,
      flags:
        d.risk.tier === "CLEAR"
          ? []
          : [
              {
                type: "clingen_regulatory_input",
                severity: d.risk.tier === "ELEVATED" ? "medium" : "low",
                ...(d.risk.worst_gene ? { gene: d.risk.worst_gene } : {}),
                detail: d.risk.reason,
              },
            ],
    })),
    // The run supplies neither an arbitrary 0–100 composite nor a method recommendation.
    ranking: [],
    analysis: {
      version: 1,
      run_id: `${key}-k562-${run.provenance.retrieved}`,
      source_path: `data/derived/demo_${key}.json`,
      source_snapshot: `pipeline-source/demo_${key}.json`,
      source_revision: sourceManifest.revision,
      source_sha256: createHash("sha256").update(text).digest("hex"),
      coordinate_system:
        "0-based, half-open 5 kb bins; reference TSS follows RefSeq txStart/txEnd boundary convention",
      window_kind: "analysis_window",
      reference,
      annotation_retrieved: annotation.retrieved,
      metric: "relative_regulatory_input_change",
      supported_modes: [mode],
      sweep_percentiles: [50, 75, 90, 95, 99],
      designs,
      source: run,
      omitted_nonranked: Math.max(
        0,
        run.n_staples_scored - run.ranked.length - run.best_below_floor.length,
      ),
      limitations: [
        "Simulated regulatory input, not expression change. No engineered contact has been experimentally validated.",
        "K562 demonstration; disease-relevant tissue transfer has not been tested.",
        "Risk screening covers ClinGen dosage sensitivity only. COSMIC and DepMap were not assessed.",
        "The strongest natural contact at this separation is the reference setting, above the p99 sweep setting.",
        "5 kb candidate anchor regions; guide sequences, PAM availability and binding off-targets have not been assessed.",
        ...(mode === "down" ? ["Decoy-mediated lowering relies on the model's contact-conservation assumption; this mechanism is unvalidated. The enhancer is tethered to a decoy, not the target promoter."] : []),
      ],
    },
  };
  await writeFile(
    new URL(`site/public/data/${key}.json`, root),
    JSON.stringify(record, null, 2) + "\n",
  );
  console.log(
    `Exported ${symbol}: ${accepted.length} shortlisted, ${designs.length - accepted.length} retained non-ranked designs`,
  );
}

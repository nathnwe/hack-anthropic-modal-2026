# Rewire Bio — independent frontend

An Astro static site developed on `feat/site-alt`, independently of `feat/site`,
and integrated into `main`. The landing page, workbench, results and API pages share
a parchment/olive identity, original interwoven DNA SVG mark and responsive CSS.
The frontend serves static regulatory records plus live human reference annotations; it does not run backend inference.

## Run locally

From the repository root, using the existing conda environment:

```sh
conda activate genome3d
npm ci --prefix site
npm run dev --prefix site -- --port 4322 --host 127.0.0.1
```

Open `http://127.0.0.1:4322/hack-anthropic-modal-2026/`.
The GitHub Pages base path is already configured. This worktree uses port 4322;
the other frontend can keep using 4321. Astro 7 may run the dev server in the
background; use its `astro dev status` command from `site/` to inspect it.

## Public deployment

Public website: [Rewire Bio](https://leandre-tappenden.github.io/rewire-bio-site/).

The [publishing repository](https://github.com/Leandre-Tappenden/rewire-bio-site)
builds this repository's `main` branch and deploys it under Leandre's account,
where Pages settings are controlled independently of the shared project.
Application development stays in this repository's `site/` directory.

After pushing website changes to `main`, open
[Publish Rewire Bio](https://github.com/Leandre-Tappenden/rewire-bio-site/actions/workflows/pages.yml)
and choose **Run workflow**, or use the signed-in GitHub CLI:

```sh
gh workflow run pages.yml --repo Leandre-Tappenden/rewire-bio-site --ref main
```

Source pushes do not automatically trigger the separate publishing repository.
The publishing workflow checks, tests and builds the latest `main` with Node 24
before deploying, and records the source revision in its build summary. Failed
builds leave the previous published version intact. No personal access token is
stored in either workflow.

The source repository still checks site changes automatically. Its original
Pages destination is now a manual-only option, because that repository's Pages
setup has not been completed. Use the publishing workflow linked above for the
active public website.

The alternate build uses `REWIRE_SITE_URL=https://leandre-tappenden.github.io`
and `REWIRE_BASE_PATH=/rewire-bio-site`. Without those variables, local development
keeps the original `/hack-anthropic-modal-2026` path on port 4322.

LDLR, NSD1, RUNX1, NF1 and MECP2 have precomputed K562 upregulation analyses; MECP2 and MYC also have decoy-based downregulation analyses. MYC retains its separate earlier rE2G case study; SHANK3 and PMP22 remain labelled fixtures. All 3D folding is illustrative.

## Pages

- `/`: full-viewport wordmark, About, Team and Contact anchors, scroll-linked DNA
  thread. Reduced-motion preference renders the thread without animation.
- `/workbench/`: live human gene autocomplete (MyGene.info) and an up/down query.
  Selecting a suggestion resolves GRCh38 coordinates through Ensembl. Example
  buttons explicitly choose the supplied samples and their corresponding direction.
- `/results/?gene=SHANK3&intent=up`: ranked elements beside a linked contact
  viewer. The default viewer is an interactive molecular-style 3D schematic with
  rotation, zoom, expansion and a controllable contact preview. The supplied
  contact map and stacked IGV genome browsers are separate tabs. The viewer starts on the top-ranked element, resets when the sort changes,
  and follows element selection. Only the separate + button opens details.
  Legacy candidate staples expand on demand. Pipeline analyses instead show the
  selected enhancer’s **Stapling strategy** beneath the list: green A/B target
  regions, compact suggested approaches and an optional anchor overlay.
- `/api/`: static JSON endpoints and contract documentation.

All paths above are under `/hack-anthropic-modal-2026` locally and
`/rewire-bio-site` on the public deployment linked above.

## Data and scientific boundaries

The three legacy files (`myc.json`, `shank3.json`, `pmp22.json`) in `public/data/`
are exact copies of the contract fixtures and carry `illustrative: true`. Placeholder warnings cover disease
annotations, cell type, coordinates, scores, tiers, risk flags and recommendations.
Regulatory predictions remain precomputed. Live reference lookups provide gene
annotations, never inferred enhancers, contacts or staples.

- Composite rankings and A/B/C tiers belong to **staples**, not CREs. Staples are
  sorted by supplied composite score, but shown with ordinal tiers.
- CREs may be sorted by supplied contact strength, raw rE2G score or TSS distance.
  Colours interpolate from pale yellow to rust using the actual min/max values
  among the displayed hits (nearest is darkest for distance). Ties have identical
  colours; missing values are neutral. The same colour follows each hit through
  the list, 3D viewer, contact map and IGV interval. rE2G is not labelled potency
  or interpreted as a probability of engineering success.
- Composite, potency and achievability sorts are disabled for CREs because the
  contract does not provide those element-level metrics.
- A contact is associated with a CRE only when one endpoint equals the gene TSS
  and the other is inside the CRE interval, using inclusive endpoints. No bin
  width or proximity tolerance is fabricated. Missing evidence is distinct from
  zero and sorts last. Maximum supplied strength is used for multiple matches.
- The contact plot uses genomic position on the x axis. Arc height is schematic;
  it is not a probability or expression-effect measurement. Out-of-domain
  contacts are omitted with a visible count. Mark types and labels supplement colour.
- Directions filter staples; they never relabel an incompatible design. An empty
  result is valid. Unresolved genes, invalid modes and failed loads have error states.
- `gene.evidence` is labelled gene-level context, not proof about a specific CRE.
  Existing URLs and PMIDs are linked; placeholder strings are visibly flagged.
- The schema does not define genome assembly, coordinate convention, contact
  resolution or tier semantics. The frontend does not supply these missing facts.

To integrate new records, put schema-conforming JSON in `public/data/` and rebuild.
The manifest is generated from those files. Do not remove placeholder status
without verified replacement data. The method library in `src/data/strategies.ts`
has primary-paper links and is separate from locus-level predictions. Its
attractions/issues are qualitative design considerations.

## K562 pipeline integration

The default workbench examples are **LDLR (1 design), NSD1 (3), RUNX1 (1),
NF1 (4), MECP2-up (2), MECP2-down (1), and MYC-down (1)**: seven runs, six genes, 13 shortlisted designs. Search and the existing direction selector resolve the matching run. Unsupported directions retain the unavailable state; explicit earlier MYC samples still load the rE2G dataset. These are precomputed results, not live inference.

Refresh the pinned source snapshot after fetching the shared repository:

```sh
git fetch origin
npm run export:pipeline --prefix site -- --source-ref origin/main
```

`site/scripts/export-pipeline.mjs` imports `data/derived/demo_*.json` from the supplied Git revision into `public/data/pipeline-source/`, records that revision and source SHA-256 hashes, then converts the seven completed demos. Running without `--source-ref` regenerates from the pinned local snapshot offline. The partial ten-record screen and combined website export are also preserved there as source JSON; incomplete screening summaries are not promoted to full demo analyses. It uses the checked-in UCSC
RefSeq snapshot `scripts/reference-annotations.json` to reproduce the pipeline's
most-frequent-isoform TSS convention (including the minus-strand txEnd boundary).
Only an intentional `--refresh-annotations` fetches new UCSC annotations; ordinary
exports and builds need no annotation network access. Review annotation changes
against the pipeline before accepting a refresh.

The shared schema is unchanged. Gene, locus, shortlisted staples and compatible
risk flags use its existing fields. An explicit `analysis.version: 1` extension
preserves the native designs, source ranks, sweeps, numeric model output and risk
labels. `locus.tad` is a compatibility interval marked `scope_kind: analysis_window`,
never presented as a called TAD. `cres`, `contacts` and `ranking` remain empty when
the export lacks rE2G scores, contact matrices, an arbitrary composite or a method
recommendation. Nothing fills those fields with zero or invented values.

- **Enhancer regions** are the unit of selection, ordered only by their associated
  design's reference-setting ΔT/T. Clicking selects; only + expands the score.
  ΔT/T is simulated **regulatory input**, not expression or an intrinsic enhancer
  activity score. Colour encodes the relative ΔT/T across the displayed shortlist.
- The selected enhancer's single exported anchor pair appears directly below the
  list. Anchor coordinates are 5 kb candidate regions; no guide sequences, PAMs
  or exact binding sites were designed. Native source ranks remain in JSON.
- Model confidence, neighbouring-gene screening and contact-strength sweeps expand
  within that design panel. The stronger reference setting is identified separately.
  Settings are neither probabilities nor experimental doses. In particular, NSD1
  enhancer 1's design is positive only at the reference setting.
- CLEAR is displayed as **No other ClinGen genes**, never low risk or safe. LOW,
  ELEVATED and BLOCKED remain distinct. Compatibility flags map LOW to low and
  ELEVATED to medium; BLOCKED never enters the shortlist/staples. Retained rejected
  and below-floor designs, and counts of omitted designs, remain in View JSON.
- **Visualise staple design** defaults off. The 3D view then retains the earlier two-window
  enhancer–TSS close-up. Turning it on adds local windows around anchor A/B centres
  and joins those centres with the proposed tether. These can differ from the
  enhancer and TSS. Intervening sequence is compressed with exact excluded bp
  counts. A/B denote 5 kb regions, not guide sites or atomic binding positions.
  The animation is illustrative, not a measured structure or molecular simulation.
- Downregulation retains negative ΔT/T and sorts by reduction magnitude. Anchor B is labelled **Decoy region**, never near the target promoter. In 3D, the design overlay previews the actual enhancer–decoy pair; without it the contact animation is disabled for these runs. Interleaved local windows are partitioned by anchor proximity to avoid duplicating source sequence. The supplied p99 value appears in the existing enhancer disclosure and the full sweep remains in Model checks.
- The same toggle adds/removes the actual A/B intervals in the genomic map and both
  IGV panes. The map uses true genomic coordinates and a schematic dashed arc,
  without fabricating a Hi-C heatmap. IGV panes use a checked-in GRCh38
  chromosome-size reference and do not show nucleotide sequence. Cached RefSeq
  gene models cover the six pipeline analysis windows, with live Ensembl
  annotations outside those windows (see below). The sourced pipeline RefSeq
  TSS is always shown, independently of anchor visibility.
- CLOuD9, LADL and BPCL appear as compact expandable literature options below the
  staple target regions. They are not per-locus recommendations: compatibility and
  binding off-targets have not been assessed. There is no bottom strategies table.

Source-hash, shortlist, risk-gate, confidence, coordinate and geometry checks run
with `npm test`. Contract validation also covers the generated records. The
frontend changes stay inside `site/`; pipeline/scoring code and shared contracts
are unchanged.

## MYC / K562 case study

The workbench currently displays **three** predictions from the same Extended
export, `ENCFF269DKY`, to keep viewer refinement focused: `K562-51` (6,134.5 bp),
`K562-49` (162,522.5 bp), and `K562-55` (1,846,075 bp) from the reference TSS.
This is a display selection, not a biological shortlist. The full source data
remains available through View JSON. `K562_VIEWER_IDS` controls this selection.

`public/data/case-studies/myc-k562.json` retains all **56** K562 rows from
`MYC_rE2G_all_cell_types_GRCh38.csv`: **53** distinct intervals and **8** source
exports. Rows are not deduplicated, merged or treated as independent validation.
Each has its original CSV row number, accession, URL, condition annotation, model
variant and raw score. The input CSV SHA-256 is included. Regenerate with:

```sh
python site/scripts/import-k562.py input.csv ensembl-lookup.json
```

The second input is the saved public Ensembl response from
`https://rest.ensembl.org/lookup/id/ENSG00000136997?expand=1;content-type=application/json`.
The reference is canonical transcript **ENST00000621592.8**, plus strand,
GRCh38 chr8:127736231–127742951 in Ensembl's 1-based coordinates, converted to
BED `[127736230, 127742951)`. Distances use **127736230** to each element midpoint;
they are reference-annotation distances, not a claim about the rE2G model's
undisclosed target TSS or a measured K562-specific TSS.

This is a separate source dataset, not a replacement shared-contract record.
`k562.ts` adapts it only for the frontend: empty contacts, staples, risks and rankings;
no invented disease or dosage claim. The internal plotting bounds cover the source
intervals and are labelled a **display span, not a TAD**. The CSV supplies no TAD,
contact matrix or measured 3D conformation. Original contract fixtures stay untouched.
The API links both this case study and the original fixtures.

All K562 rows are potential regulatory elements, including sodium butyrate,
vorinostat and DMSO conditions. An absent treatment annotation means **not listed**,
not verified untreated. A score sort compares supplied numbers, not calibrated
cross-model evidence. A distance sort is also available; contact-strength sorting
is disabled because measured contacts are absent. Dashed map arcs are explicitly
predicted regulatory links. Source: [ENCODE-rE2G](https://www.nature.com/articles/s41586-026-10781-4)
and the per-row ENCODE accessions.

## 3D viewer boundaries

The default **3D DNA** view opens directly on PDB-style molecular detail.
It shows two 160 bp local windows around the reference TSS and the selected
regulatory-element midpoint: 64 bp on the outer side and 96 bp toward the gap,
giving the right-hand loop more contour and a wider arc. The exact inner gap is
calculated from those window boundaries. The DNA fades progressively over the
inner 34 bp into a soft haze labelled with the **exact omitted count** and no
secondary caption. No molecular bonds bridge that omitted interval.
Shown bp + omitted bp equals the source span covered by this close-up;
source coordinates and counts never change during the contact preview.

The local centrelines are curved, three-dimensional illustrative folds. The two
visible fragments have one rung per local bp, approximately 0.34 nm centreline
rise and 10.5 bp per helix turn. Beads suggest molecular detail; they are not a
sequence-derived atomic model or a PDB structure. The coloured patches locate
20 bp around the TSS and element midpoint; they do not imply that the entire gene
or regulatory element is shown. Larger source gaps have a broader haze on a
logarithmic **display** scale, with stable molecular magnification across the three
examples. Neither the haze size nor the arrangement measures physical compaction.

The user-provided TAD illustration informs the fluid visual composition only.
**No TAD boundary, chromatin conformation or real contact is inferred.** The CSV
supplies regulatory predictions and genomic coordinates, not these structural
measurements. Following the user's loop sketch, the promoter stays anchored on
the lower strand while the enhancer approaches from above along a curved path.
Each hit gets a distinct starting enhancer fold seeded from its source coordinates.
This is repeatable visual variation, not an inferred structural difference. The
variation settles away during the approach, preserving the common final layout.
The omitted interval occupies the right-hand arc, with soft fades at the outer
window boundaries. Unit-tangent integration bends the centrelines while preserving
local contour spacing. Smooth easing and small, damped travelling bends make the
seven-second transition fluid; every frame is deterministic and can be scrubbed
or paused without jitter. The preview is not a molecular-dynamics simulation and
does not model nucleosomes, forces, temperature, an intervention or expression
effects. The final 3.5 nm centreline separation is an illustrative viewing choice.

The earlier focus/scale shortcuts are removed from the interface. The full-span
renderer remains only as the fallback for overlapping local windows.

The viewer supports rotation, zoom/pan, expansion, keyboard controls, a contact
slider and reduced motion. Overlapping TSS/element intervals remain viewable in
the full model with looping disabled. Invalid coordinates retain the 2D map
fallback. Shader alpha hashing makes the molecular fade gradual, with a procedural
smoke texture and no external image asset. The separate Three.js lazy chunk
exceeds the build's 500 kB advisory threshold.

Nominal helix dimensions: [DNA mechanics and topology](https://pmc.ncbi.nlm.nih.gov/articles/PMC7288220/).

## Live search and genome browsers

The optional live-reference URL is `/results/?gene=BRCA1&source=reference&intent=up`.
Typed queries/suggestions use this route; sample buttons retain the existing
precomputed-record flow. An unknown symbol also attempts live resolution.

- [MyGene.info](https://docs.mygene.info/en/latest/doc/query_service.html) supplies
  debounced human symbol/alias/name suggestions, with cancellation and caching.
- [Ensembl REST](https://rest.ensembl.org/) resolves genes and canonical transcript
  TSS positions. Only human GRCh38 primary chromosomes are supported. Ensembl's
  1-based inclusive intervals convert to 0-based, half-open frontend intervals.
- [IGV.js](https://igv.org/doc/igvjs/) renders two independent panes: the selected
  enhancer above, the target gene below. Both pan and zoom; Fit restores the
  selected interval. Changing a hit keeps the lower pane's navigation intact.
  The selected interval uses its ranking colour; the gene stays blue. Unknown
  enhancer strand is left unspecified.
- Track names sit in a fixed left gutter, clear of the plotted intervals. Green
  Anchor A/B rows match the staple panel and molecular labels. Browser coordinates
  are shown in 1-based notation; stored intervals are 0-based, half-open.
- The checked-in `public/data/reference/hg38.chrom.sizes` supplies chromosome
  lengths; nucleotide sequence is not loaded. RefSeq gene models for the six
  pipeline windows are in `public/data/reference/k562-refseq.json`, sourced from
  UCSC's `ncbiRefSeqCurated` track. Each window records its source URL and payload
  SHA-256. One representative per gene/strand is chosen by preferring NM_ coding
  transcripts, then longest spliced length, genomic span and accession. These
  are **not claimed to be canonical** and do not replace the pipeline's TSS
  convention. Outside the cached coverage, Ensembl supplies live canonical models.
- Non-illustrative records with explicit **K562 / GRCh38** context load the four
  ENCODE signal tracks below in both panes. Other cell types, unknown contexts
  and fixtures do not receive these tracks. They provide experimental context;
  they do not alter the shortlist, scores or proposed designs.
- The IGV code loads only when its tab is opened. API failures have visible
  error/retry states. External APIs and the reference host require a network
  connection and CORS support; no credentials or backend service are required.

### K562 experimental tracks

Pinned released GRCh38 bigWig files from each experiment's default analysis:

| Track | Experiment | Signal file |
| --- | --- | --- |
| H3K4me3 ChIP-seq | ENCSR668LDD | [ENCFF253TOF](https://www.encodeproject.org/files/ENCFF253TOF/) |
| H3K27ac ChIP-seq | ENCSR000AKP | [ENCFF381NDD](https://www.encodeproject.org/files/ENCFF381NDD/) |
| H3K27me3 ChIP-seq | ENCSR000AKQ | [ENCFF242ENK](https://www.encodeproject.org/files/ENCFF242ENK/) |
| ATAC-seq | ENCSR868FGK | [ENCFF102ARJ](https://www.encodeproject.org/files/ENCFF102ARJ/) |

These are pooled biological replicates (1/2 for H3K4me3; 1/2/3 for the others),
reported as **fold change over control/background**, not probabilities. Tracks
stream public byte ranges only when the browser opens, with mean aggregation and
independent vertical autoscaling per track and pane. Heights are therefore not a
shared quantitative scale. Links, output units and ENCODE audit flags are available
under **Tracks & sources**. The source metadata, checksums and pinned cloud URLs
are recorded in `src/data/k562-tracks.json`.

The selected experiments have quality flags: H3K27ac includes an extremely-low-read-depth
error; H3K27me3 includes insufficient-depth noncompliance; H3K4me3 has control-related
warnings; ATAC-seq has a library-complexity warning. They remain visible in the
sources disclosure and should be considered when interpreting the signals.
A failed signal load receives one retry, then a visible error rather than a
fabricated zero track. Cached gene models and selected intervals remain usable
independently of the remote signal service.

Live genes without supplied regulatory records open their gene annotation and
show an empty enhancer pane. Reference annotations do not predict regulatory
activity. Illustrative fixtures lack a verified assembly and cannot be overlaid
on real sequence. The sourced MYC sample has verified GRCh38 intervals. Source
records remain in JSON, although CSV/project-specific detail is removed from
the workbench interface. The 3D molecular geometry remains illustrative.

## Validation

```sh
npm test --prefix site
npm run check --prefix site
npm run build --prefix site
uv run python contracts/validate.py
```

Tests cover exact source counts, provenance, source score preservation, distinct
ranking orders, all 56 locus windows, independent numerical contour integration,
base-pair spacing at kilobase and megabase scales in both genomic directions,
finite frames, both arms of the junction, overlap handling and invalid coordinates.
Existing tests retain fixture contact matching, empty results, direction filtering,
escaping and selection behavior. Folded-view tests cover exact omission accounting,
progressive fading, fixed promoter placement, non-rigid motion, local contour
spacing, smooth endpoints, source order and clipped windows. Browser
QA covers all three selected examples, independent selection/details, score and
distance ranking, live search, IGV interval switching, zoom and mobile layout.

`npm audit --omit=dev` currently reports the upstream `http-cache-semantics`
max-stale advisory, also attributed to Astro. The registry's latest versions
(Astro 7.3.5, http-cache-semantics 4.2.0 at implementation time) do not offer a
compatible patched release. No forced downgrade was applied. This application
builds static files and has no authenticated server or cross-user response cache.
Recheck the dependency before adding server rendering or deploying a server.

## Credits and remaining content

- Astro (MIT) provides the static build; TypeScript and Node's test runner provide
  checks. Styling uses ordinary CSS, with no UI component library.
- Three.js and OrbitControls (MIT) provide the interactive 3D renderer and camera.
- IGV.js (MIT) provides the genome browsers; MyGene.info, Ensembl and UCSC supply
  the public reference data described above.
- DM Sans and IBM Plex Mono are loaded from Google Fonts (SIL Open Font License),
  with system fallbacks when unavailable. Serif emphasis uses system Georgia.
- The logo, conceptual diagrams and DNA thread are original SVG/CSS artwork.
- The landing-page loop panel is adapted from the team's `LoopPanel.astro` on
  `feat/site` (9732ea7), retaining its scroll-driven expansion and blue/copper
  palette with a warm contact glow. A lazy-loaded Three.js scene now renders
  shaded molecular spheres and bonds around a continuously bending double helix.
  Its travelling bends and loop closure are illustrative, not molecular dynamics,
  a measured conformation, an atomic structure or a contour-conserving simulation.
  Rendering stops off-screen, in background tabs, while paused and during video
  playback. Reduced motion shows a still folded model. The locally generated
  `public/images/landing-dna.svg` is the static fallback when WebGL is unavailable.
- “visit the cell” opens the user-supplied film in a full-screen modal with a
  crossfade and subtle zoom. The replacement 40-second, 1080p H.264 MP4 is stored at
  `public/media/rewire-bio-cell.mp4` (no audio track; fast-start metadata retained).
  Source: `WhatsApp Video 2026-10-04 at 12.27.57.mp4`, copied without transcoding.
  The video URL carries a content-version query to refresh earlier cached copies.
  Video loading begins on click. Native playback controls, a Close button and
  Escape support playback and return to the same page position; reopening starts
  the film again. Keyboard focus stays in the modal and returns to the launch
  link on close. Reduced motion removes the transitions. With JavaScript disabled,
  the link opens the video file directly.
- Landing-page scientific context follows the project background and
  [Morgan et al., 2017](https://www.nature.com/articles/ncomms15993).
- Fixture provenance: `contracts/fixtures/`; source of truth: `contracts/schema.json`.
- The method library links the primary studies for CLOuD9, LADL, bivalent dCas,
  BPCL, LDB1, CTCF and ZF tethering.
- Team profiles link to public professional pages. Biographical sources are in
  `src/data/team.ts`; portrait sources and cropping are recorded in
  `public/images/team/CREDITS.md`. Experimental findings and verified gene
  records remain separate from these profiles and must come from the pipeline.

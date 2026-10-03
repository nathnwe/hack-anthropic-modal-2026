# Rewire Bio — independent frontend

An Astro static site developed on `feat/site-alt`, independently of `feat/site`,
and integrated into `main`. The landing page, workbench, results and API pages share
a parchment/olive identity, original interwoven DNA SVG mark and responsive CSS.
The frontend serves static records; it does not run backend inference.

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

Public URL after Pages activation: [Rewire Bio](https://nathnwe.github.io/hack-anthropic-modal-2026/).

One-time setup requires a repository administrator, maintainer, or someone with
permission to manage Pages: open **Settings → Pages → Build and deployment** and
select **GitHub Actions** as the source. Then run **Actions → Deploy site to
GitHub Pages → Run workflow**, selecting `main`.

The workflow in `.github/workflows/pages.yml` automatically checks, builds and
deploys changes to `site/` or that workflow when they reach `main`. Local changes
and feature-branch pushes stay out of the public site until merged. A failed
check or build stops deployment and leaves the previous successful site live.
The workflow uses Node 24 and runs the type checks and frontend tests before
building the static site. The existing contract-validation workflow is separate.

MYC uses sourced K562 rE2G predictions; the other examples remain labelled fixtures. All 3D folding is illustrative.

## Pages

- `/`: full-viewport wordmark, About, Team and Contact anchors, scroll-linked DNA
  thread. Reduced-motion preference renders the thread without animation.
- `/workbench/`: gene lookup and up/down query. Example buttons select a
  fixture and its corresponding direction. Unsupported genes show an error.
- `/results/?gene=SHANK3&intent=up`: ranked elements beside a linked contact
  viewer. The default viewer is an interactive molecular-style 3D schematic with
  rotation, zoom, expansion and a controllable contact preview. The supplied
  contact map remains available as a separate tab. The viewer starts on the top-ranked element, resets when the sort changes,
  and follows element selection independently of whether its details stay open.
  Candidate staples and method tradeoffs expand on demand.
- `/api/`: static JSON endpoints and contract documentation.

All paths above are under `/hack-anthropic-modal-2026` in development and production.

## Data and scientific boundaries

The three files in `public/data/` are exact copies of the contract fixtures.
They all carry `illustrative: true`. Placeholder warnings cover disease
annotations, cell type, coordinates, scores, tiers, risk flags and recommendations.
The workbench currently serves precomputed records only.

- Composite rankings and A/B/C tiers belong to **staples**, not CREs. Staples are
  sorted by supplied composite score, but shown with ordinal tiers.
- CREs may be sorted by supplied contact strength or raw rE2G score. rE2G is not
  labelled potency or interpreted as a probability of engineering success.
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
  result is valid. Unknown genes, invalid modes and failed loads have error states.
- `gene.evidence` is labelled gene-level context, not proof about a specific CRE.
  Existing URLs and PMIDs are linked; placeholder strings are visibly flagged.
- The schema does not define genome assembly, coordinate convention, contact
  resolution or tier semantics. The frontend does not supply these missing facts.

To integrate new records, put schema-conforming JSON in `public/data/` and rebuild.
The manifest is generated from those files. Do not remove placeholder status
without verified replacement data. The method library in `src/data/strategies.ts`
has primary-paper links and is separate from locus-level predictions. Its
attractions/issues are qualitative design considerations.

## MYC / K562 case study

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

The full interval from the reference TSS to the complete selected element, plus
bounded flanks, is represented continuously. **No ellipsis and no omitted bp.**
The distance readout is TSS-to-midpoint; the total model count includes flanks.
World units are nanometres, using idealised B-DNA with 0.34 nm centreline rise per
bp, 10.5 bp per turn and a 1 nm helix radius. These are nominal modelling dimensions,
not sequence-specific atomic coordinates. See [DNA mechanics and topology](https://pmc.ncbi.nlm.nih.gov/articles/PMC7288220/).

A curve indexed by cumulative arc length preserves DNA contour length during the
contact preview. The proposed junction approaches a **chosen illustrative 3 nm
centreline separation**, not an inferred biological distance or staple design.
The loop is explanatory geometry, not a measured conformation or molecular dynamics.
It omits nucleosomes and higher-order chromatin packing. Actual in-cell spatial
separation cannot be recovered from this CSV.

At overview scale a continuous centreline represents the entire DNA; it and the
position markers have a minimum pixel size for visibility. As the camera zooms in,
the visible region resolves into a molecular-style helix with **one rung per bp**.
Only visible molecular detail is allocated; the continuous model remains present
outside the viewport. Full span, gene, element and junction shortcuts connect these
scales. A small full-span inset retains context during close-ups. Blue shows the
reference transcript (when supplied), rust the selected interval; the TSS marker
has priority where they overlap. All counts and genomic coordinates remain fixed
while folding; only the illustrative conformation changes.

The viewer supports drag rotation, zoom/pan, expansion, keyboard controls, a
contact slider and reduced motion. Overlapping TSS/element intervals remain
viewable with looping disabled. Invalid coordinates retain the 2D map fallback.
The separate Three.js lazy chunk exceeds the build's 500 kB advisory threshold.

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
escaping and selection behavior. Browser QA clicks every K562 row and checks
molecular detail, longest-span selection, sorting, map fallback and mobile layout.

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
- DM Sans and IBM Plex Mono are loaded from Google Fonts (SIL Open Font License),
  with system fallbacks when unavailable. Serif emphasis uses system Georgia.
- The logo, conceptual diagrams and DNA thread are original SVG/CSS artwork.
- Fixture provenance: `contracts/fixtures/`; source of truth: `contracts/schema.json`.
- The method library links the primary studies for CLOuD9, LADL, bivalent dCas,
  BPCL, LDB1, CTCF and ZF tethering.
- Real team names, contact details, experimental findings and verified gene
  records remain to be supplied. The GitHub project is the current contact link.

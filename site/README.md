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

The default **Folded DNA** view opens directly on PDB-style molecular detail.
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

**Genomic span** retains the earlier full-length view: a continuous centreline,
with one consistent nanometre scale and no internal omission. DNA resolves into
individual bp as the camera zooms in. Its overview line and position markers have
a minimum pixel size for visibility. A small whole-span inset accompanies its
molecular close-ups. The default remains Folded DNA, including after switching
regulatory elements or resetting the camera. Gene, element and junction shortcuts
return to the folded molecular view.

The viewer supports rotation, zoom/pan, expansion, keyboard controls, a contact
slider and reduced motion. Overlapping TSS/element intervals remain viewable in
the full model with looping disabled. Invalid coordinates retain the 2D map
fallback. Shader alpha hashing makes the molecular fade gradual, with a procedural
smoke texture and no external image asset. The separate Three.js lazy chunk
exceeds the build's 500 kB advisory threshold.

Nominal helix dimensions: [DNA mechanics and topology](https://pmc.ncbi.nlm.nih.gov/articles/PMC7288220/).

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
QA covers all three selected examples, sorting, both viewing scales and mobile layout.

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

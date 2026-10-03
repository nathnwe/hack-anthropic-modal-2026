# Rewire Bio — independent frontend

An Astro static site on `feat/site-alt`, developed independently of `feat/site`.
Only `site/` is changed. The landing page, workbench, results and API pages share
a parchment/olive identity, original interwoven DNA SVG mark and responsive CSS. No deployment
or backend inference is performed by this frontend.

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

## 3D viewer boundaries

The Three.js viewer is loaded on demand on the results page. Genomic positions
identify the two highlighted regions within a local window. The centreline,
helix turns, beads and bonds are generated illustration geometry: no sequence,
atomic coordinates, measured spatial distances or folding prediction is supplied.
The gene marker denotes the supplied transcription start, not a fabricated gene
body. Highlight widths are enlarged for visibility. Bringing the markers together
is a hypothetical contact exploration, independent of expression direction or a
method recommendation; it does not imply the chosen direction has a candidate.

The underlying DNA stays connected through the illustrative deformation. An
invalid, overlapping or out-of-domain pair offers the 2D contact map instead.
The viewer supports drag rotation, zoom/pan, keyboard controls, a contact slider
and reduced-motion preferences. It renders on demand and pauses off screen or in
a hidden tab. A WebGL failure retains access to the contact map. The 3D bundle is
larger than the build tool's 500 kB advisory threshold and is a separate lazy chunk.

## Validation

```sh
npm test --prefix site
npm run check --prefix site
npm run build --prefix site
uv run python contracts/validate.py
```

Unit tests exercise missing-versus-zero contact values, exact endpoint matching,
direction filtering, non-mutating supported sorts, viewer defaults under opposing
rankings, source handling and HTML escaping. DNA tests cover genomic mapping,
invalid anchors, curve continuity, approaching anchors and finite helix frames. Browser QA covers desktop and phone layouts, the full query flow,
unsupported genes, evidence expansion, direction changes, empty results and API
links. Candidate tables have contained horizontal scrolling on phones; the contact map
resizes to keep the selected element visible.

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

No push, pull request or GitHub Pages deployment has been performed.

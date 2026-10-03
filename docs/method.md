# Method

From DNA sequence and ATAC to a ranked shortlist of staple targets.

This is the locked spec. `docs/pipeline-design.md` is the brainstorm that led here and is kept
only as a record. Every equation below is implemented in `pipeline/` and `scoring/`; every number
quoted in a figure or a slide comes from `pipeline/simulate.py`, so the prose and the pictures
cannot drift apart.

Each step has the maths, then a plain-English line, then the data it needs.

**What a staple is.** A pair of dCas proteins joined by a linker, each guided to a chosen spot in
the genome, physically pulling those two spots together. It changes 3D contact. It does not change
the DNA sequence and it does not change chromatin accessibility. That last fact is the whole
reason this method has the shape it does.

---

## Notation

| symbol | meaning | units |
|---|---|---|
| `e` | a candidate regulatory element (enhancer) | bin index |
| `G` | a gene, located at its TSS | bin index |
| `A_e` | activity of element `e` | read counts |
| `C(e, G)` | contact frequency between `e` and `G`'s promoter | normalised |
| `T_G` | total regulatory input to gene `G` | `A · C` |
| `q` | strength of the forced contact | same units as `C` |
| `σ` | width of the staple's influence | bins |
| `τ` | how much a flagged gene may move before we refuse | fraction |

Resolution is 5 kb throughout, matching ABC's published definition and UniversalEPI's output.
The scoring window is 2 Mb centred on the target.

---

## 1. Build one contact map

$$C \;=\; \begin{cases}
\text{KR-normalised Hi-C at 5 kb} & \text{if measured Hi-C exists for this cell type}\\[2pt]
\text{UniversalEPI}(\text{sequence},\ \text{ATAC}) & \text{else if ATAC exists}\\[2pt]
\text{ENCODE Hi-C megamap (ENCSR492BKE)} & \text{else}\\[2pt]
\big(d+1\big)^{-1} & \text{last resort}
\end{cases}$$

**In plain terms.** We need a number for how often any two points in the region touch. Measured
Hi-C is the real thing, so we use it where it exists. Where it does not, UniversalEPI predicts it
from the DNA sequence plus an ATAC-seq track, which is cheap and available for far more cell types.
Failing that we fall back to an average across many cell types, and finally to the fact that nearby
things touch more than distant things.

The rung matters and is recorded with every result:

| contact source | ENCODE-rE2G precision at 70% recall |
|---|---|
| cell-type Hi-C | **53.3%** |
| averaged megamap | 48.0% |
| inverse distance | 44.3% |

**Data.** Hi-C: ENCODE / 4DN, KR-normalised, 5 kb. ATAC: ENCODE DNase or ATAC-seq for the cell
type. Sequence: hg38. Promoters: GENCODE v29, 500 bp centred on the RefSeq TSS with the most
coding isoforms.

> **Figure 1** — `docs/figures/method/fig1_contact_source.png`

---

## 2. Activity-by-Contact, and the two numbers it gives you

Fulco et al. 2019, *Nat Genet* 51:1664–1669.

$$A_e \;=\; \sqrt{\mathrm{DNase}_e \times \mathrm{H3K27ac}_e}
\qquad\text{(we default to }A_e = \mathrm{ATAC}_e\text{; see below)}$$

$$\boxed{\;T_G \;=\; \sum_{e} A_e \cdot C(e, G)\;}
\qquad\qquad
\boxed{\;\mathrm{ABC}_{E,G} \;=\; \frac{A_E \cdot C(E, G)}{\sum_e A_e \cdot C(e, G)}\;}$$

Sum over all elements within the window. Published threshold for calling a link: `ABC ≥ 0.02`.

**In plain terms.** An enhancer's influence on a gene is how loud it is times how often it is in
touch. Loudness is `A`; touching is `C`. Multiply, and you have that one enhancer's contribution.

Then the two numbers, which are easy to confuse and answer different questions:

- **`T_G` is the dose.** Add up every enhancer's contribution and you have the total regulatory
  input the gene is receiving. This is the dial a staple is supposed to turn. We report `ΔT/T`.
- **`ABC_{E,G}` is the attribution.** Divide one enhancer's contribution by the total and you get
  its share of the credit. It is normalised per gene, so it always sums to 100%. Read it after an
  intervention and you learn about redistribution, not increase.

Getting these the wrong way round is the single easiest way to produce a confident wrong answer.

**On H3K27ac.** The published activity term needs a ChIP-seq track that most cell types do not
have. We default to accessibility alone, which is what ENCODE-rE2G's own core model uses, and
accept the precision cost. Pass H3K27ac to recover the published form; the code supports both
(`pipeline/abc.py`).

**What ABC was validated on.** CRISPRi-FlowFISH across more than 3,500 element–gene connections.
ABC was built to predict quantitative effects, which is why it, and not ENCODE-rE2G, carries the
magnitude in this pipeline. rE2G returns a calibrated *probability* that CRISPRi of an element
changes a gene's expression — its score against measured effect size is only `r = −0.45`.

> **Figure 2** — `docs/figures/method/fig2_abc_anatomy.png`

---

## 3. Filter the physical conflicts

$$\mathcal{F} \;=\; \{e\} \,\cup\, \{\mathrm{TSS} \pm 5\,\mathrm{kb}\} \,\cup\, \{\text{gene bodies}\} \,\cup\, \{\text{CTCF motif hits}\}$$

$$\mathcal{A} \;=\; \{\,b \;:\; \mathrm{ATAC}(b) \;\wedge\; b \notin \mathcal{F} \,\}$$

**In plain terms.** The staple has to grab something, and it grabs open chromatin. But parking a
large protein on a working part of the genome breaks that part. So we subtract everything we would
damage — the enhancer itself, the promoter and 5 kb either side, the transcribed body of any gene,
and CTCF boundary sites — and what is left is where a staple may land.

This is a filter, not a score. It removes the impossible rather than ranking the possible, which
is why "tractable" does not need to be an axis in the ranking.

**Data.** ATAC/DNase peaks for the cell type · GENCODE v29 gene models · CTCF motif hits from
JASPAR **MA0139.1**.

> **Figure 3** — `docs/figures/method/fig3_anchor_filter.png`

---

## 4. Simulate the staple

### 4a. The edit is a blob, not a point

$$K(i, j) \;=\; \exp\!\left(-\frac{(i - A)^2 + (j - B)^2}{2\sigma^2}\right), \qquad \sigma \approx 3 \text{ bins } (15\,\text{kb})$$

$$C'[i, j] \;=\; \max\big(\,C[i, j],\; q \cdot K(i, j)\,\big)$$

**In plain terms, and this is the step the whole safety argument rests on.** ABC's denominator is
computed separately for each gene. So if you change a single cell of the contact matrix, only the
gene you aimed at can possibly move — every other gene's sum is untouched. Enhancer hijacking
becomes arithmetically impossible and your off-target analysis returns zeros no matter what. A
model that cannot express the failure mode cannot rule it out.

So the staple deposits a Gaussian blob around the anchor pair instead. Physically this is the
right picture too: tethering two points drags their neighbourhoods along.

**Where σ comes from.** ENCODE-rE2G found statistically significant CRISPRi crosstalk between
elements 1–10 kb apart, `P = 1.79 × 10⁻⁵`. σ ≈ 15 kb is consistent with that and is the most
important unvalidated number in the method.

### 4b. Contact is conserved, not created

$$C'' \;=\; \text{rescale } C' \text{ so that } \sum_j C''[i, j] = \sum_j C[i, j] \ \ \forall i,
\quad\text{symmetrised, iterated}$$

**In plain terms.** A stretch of chromosome has a finite amount of contact to give. Pulling it
towards one place pulls it away from everywhere else. Without this constraint a staple adds
closeness for free, nothing in the window can ever go *down*, and enhancer theft — where your
staple steals an enhancer another gene was using — cannot be represented at all.

### 4c. Sweep the strength rather than guessing it

$$q \in \big\{\,\text{p50},\ \text{p75},\ \text{p90},\ \text{p95},\ \text{p99 of } \{C[i, j] : |i - j| = |A - B|\}\,\big\}$$

**In plain terms.** Nobody has measured how tight a contact a real dCas staple produces. Rather
than invent a number, we run the whole thing across the range of contacts that genuinely occur at
that genomic separation, and report the curve. An unknown parameter becomes an output.

A single headline figure, where one is needed, uses the strongest real contact at that separation
— "as tight as the tightest natural loop this far apart." Describable, checkable, not a free
parameter.

> **Figures 4, 5, 6** — `fig4_point_vs_kernel.png`, `fig5_before_after_diff.png`,
> `fig6_dose_response.png`

---

## 5. Re-score every gene in the window

$$T'_G \;=\; \sum_e A_e \cdot C''(e, G)
\qquad\qquad
\frac{\Delta T}{T}\bigg|_G \;=\; \frac{T'_G - T_G}{T_G}$$

for every gene `G` whose TSS lies within 1 Mb of either anchor.

**In plain terms.** Run step 2 again on the edited map, for every gene in the neighbourhood rather
than only the one we were aiming at. Each gene's regulatory input went up, went down, or did not
move. Collateral gain and enhancer theft are simply other entries in the same table.

**What this is not.** `ΔT/T` is not a percentage change in expression. Neither ABC nor rE2G is
calibrated to expression units. It is a relative change in predicted regulatory input, usable for
ranking candidates against each other and not for promising a number to a biologist.

> **Figure 7** — `docs/figures/method/fig7_window_readout.png`

---

## 6. Safety

Flagged genes come from three fixed lists:

| flag | source |
|---|---|
| oncogene, tumour suppressor | COSMIC Cancer Gene Census |
| triplosensitive, haploinsufficient | ClinGen dosage sensitivity |
| essential | DepMap common-essential |

A candidate fails if **either** test trips, for any flagged gene `F` in the window:

$$\textbf{absolute:}\quad \left|\frac{\Delta T}{T}\right|_F > \tau, \qquad \tau = 0.10$$

$$\textbf{relative:}\quad \left|\frac{\Delta T}{T}\right|_F \;\ge\; \left|\frac{\Delta T}{T}\right|_{\text{target}}$$

**In plain terms.** Some genes are dangerous to disturb in either direction — the ones that cause
cancer when turned up, the ones that cause disease when turned down, the ones cells cannot live
without. Two separate tests, because they fail for different reasons. The first catches a flagged
gene that moved a lot. The second catches a staple that does more to an oncogene than it does to
its own target — which should kill a candidate however small both numbers are.

**Asymmetric trust.** Per `CLAUDE.md`, a language model may *raise* a concern, with a citation. It
may never clear one. Silence from a model is not evidence of safety, and no LLM output enters the
ranking.

> **Figure 8** — `docs/figures/method/fig8_safety_gate.png`

---

## 7. Confidence

Not a judgement call. One question — *does the answer survive the things we do not know?* — asked
four ways, all deterministic:

| test | passes when |
|---|---|
| ensemble abstention | UniversalEPI's K=10 deep ensemble does not return a maximum-confidence fold change of exactly 0 (which it does when its prediction intervals overlap) |
| contact provenance | the map came from measured Hi-C or UniversalEPI, not a power law |
| distance | anchors within 500 kb (every model tested failed three real MYC enhancers at ~2 Mb) |
| sweep stability | the gene ordering does not change as `q` runs from p50 to p99 |

All four → `CONFIDENT`. Abstention → `ABSTAIN`, and the candidate is withheld rather than ranked.
Anything else → `TENTATIVE`.

---

## 8. Ranking

```
1. filter      step 3 — removes the impossible
2. gate        step 6 — SAFE or the candidate is dropped, no exceptions
3. sort        by ΔT/T on the target, descending
4. label       confidence tier shown beside each row, not sorted on
5. report      anchor separation, PAM availability, provenance — raw, for the researcher
```

**Why a gate and not a weighted score.** A weighted sum would let a large predicted effect
outrank a safety concern — that is what weighted sums do. A gate cannot be traded against.
Deliverability stays out of the model entirely and is reported as facts for the researcher to
judge.

The ranking is deterministic and auditable end to end. Given the same inputs it produces the same
order, and every step of that order can be inspected.

---

## Assumptions, stated plainly

These are not measurements. They are choices, and a reviewer should be able to find them in one
place.

1. **σ ≈ 15 kb** — the width of a staple's influence. Supported by CRISPRi crosstalk at 1–10 kb
   (`P = 1.8 × 10⁻⁵`), not measured for a staple.
2. **Row-sum conservation** — that a locus has a fixed contact budget. Physically motivated,
   not verified.
3. **τ = 10%** — the point at which a flagged gene has moved too far. Arbitrary starting value.
4. **Contact substitution** — ENCODE-rE2G takes contact from Hi-C or the megamap; feeding it
   UniversalEPI predictions instead is an unvalidated swap.

---

## Limitations, for the future-work slide

- **No expression calibration.** `ΔT/T` has no units a biologist can act on. Fixing this means
  regressing it against measured CRISPRi effect sizes (Gasperini, Schraivogel).
- **No validation of the counterfactual.** Nobody has applied a staple and measured the result,
  so the central prediction is untested by construction.
- **Activation only.** rE2G is trained on activating elements and silencers are not modelled, so
  direction is baked in. We have no mechanism yet for turning a triplosensitive gene *down*.
- **Genomic off-target is out of scope.** We model the *regulatory* accident: the staple lands
  exactly where we aimed and the wrong gene responds. The staple binding somewhere else entirely
  is guide design, and we do not solve it.
- **Cell-type transfer is weak.** rE2G's precision falls to 53.2% across five held-out cell types.
  Roughly half the links we reason over are wrong outside the training cell type.
- **Enhancers are promiscuous.** 1.57 genes per enhancer on average, up to 34. An intervention is
  rarely as targeted as it looks.
- **Long range fails.** Every model tested missed three real MYC enhancers at ~2 Mb.
- **Sequence models cannot help.** Enformer and relatives largely ignore distal enhancers — about
  90% of prediction variance is driven by the promoter alone (Karollus, Mauermeier & Gagneur 2023,
  *Genome Biology*). A staple is invisible to them, since it changes neither sequence nor
  accessibility.
- **GraphReg is the architecturally correct alternative** and is not yet adopted. It takes the
  contact graph as an explicit input, which is the natural home for an intervention, but it is
  validated for motif ablation rather than edge addition and is trained per cell type.

---

## Credits

ABC — Fulco et al. 2019, *Nat Genet* 51:1664–1669 ·
ENCODE-rE2G — Gschwind et al. 2026, *Nature* 657:179–190 ·
UniversalEPI — Grover et al. 2026, *Nucleic Acids Res* 54(10) gkag485 ·
GraphReg — Karbalayghareh, Sahin & Leslie 2022, *Genome Research* 32:930 ·
COSMIC Cancer Gene Census · ClinGen dosage sensitivity · DepMap · GENCODE · JASPAR · ENCODE.

## Reproducing the numbers

```bash
uv run --with numpy python -m pipeline.simulate
uv run --with matplotlib --with numpy python animation/figures/method_figures.py
```

The demo locus (`pipeline/locus.py`) is **synthetic and marked `illustrative`** — element
positions, activities and gene positions are invented, and the contact matrix is generated from a
distance-decay power law with TAD blocks and CTCF loops. The arithmetic applied to it is the real
pipeline arithmetic. Only the inputs are made up.

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

**Diagonal.** The Hi-C diagonal is self-contact and is imputed from the two neighbouring bins
rather than zeroed. Zeroing silently drops any element sharing a bin with a promoter (5 of 76
genes in the real window; for HBE1 it understated T by 21%).

**Data.** Hi-C: ENCODE / 4DN, KR-normalised, 5 kb (Juicer SCALE balancing where KR is absent, as
in the real K562 file). ATAC: ENCODE DNase or ATAC-seq for the cell
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
\quad\text{alternately scaled and symmetrised, 100 iterations}$$

(Four iterations, as first written, left up to 0.9% row-sum error. A hundred converge to ~1e-16.)

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
anywhere in the window — "as tight as the tightest natural loop this far apart." Describable,
checkable, not a free parameter.

**Known defect, and its handling.** If the pair being stapled *is* that strongest contact, then
`max(C, q·K)` leaves the anchor pixel unchanged and renormalisation slightly *weakens* it
(verified on real data: LCR→HBG2, 49.13 → 47.62). Such a pair is reported as
**already at the yardstick — no-op** and is not scored as though the staple did something. The
yardstick also lies above the sweep's p99, so the headline is the top of the range, not its middle.

> **Figures 4, 5, 6** — `fig4_point_vs_kernel.png`, `fig5_before_after_diff.png`,
> `fig6_dose_response.png`

---

## 5. Re-score every gene in the window

$$T'_G \;=\; \sum_e A_e \cdot C''(e, G)
\qquad\qquad
\frac{\Delta T}{T}\bigg|_G \;=\; \frac{T'_G - T_G}{T_G}$$

for every gene `G` in the 2 Mb window (not only those within 1 Mb of an anchor: renormalisation
moves distant genes by up to ~0.5%, below τ, but it is real model behaviour). A gene with zero
baseline `T` makes `ΔT/T` undefined; the run refuses rather than emitting `inf`.

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

**Fail closed.** A flagged gene whose change is undefined (`NaN`, e.g. from a zero baseline) is
reported as a concern, not waved through, and a target missing from the deltas is an error rather
than a silently disabled relative test. A null staple (target change exactly 0) does not trip the
relative test on untouched flagged genes.

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
| sign stability | the predicted change on the target keeps the same sign at every `q` from p50 to p99 |

(Full gene-ordering stability was the first definition. Over 592 real staples it was satisfied by
none, so it cannot be the test; sign stability is the weaker, implementable one.)

**Which tests are live depends on the run.** In the real-data run UniversalEPI was not used, so
*abstention* is a constant and *provenance* is constant "measured Hi-C": only distance and sign
stability actually vary, and the output says so beside the tier.

All four → `CONFIDENT`. Abstention → `ABSTAIN`, and the candidate is withheld rather than ranked.
Anything else → `TENTATIVE`.

---

## 8. Ranking

```
1. filter      step 3 — removes the impossible
2. gate        step 6 — SAFE or the candidate is dropped, no exceptions
3. floor       ΔT/T on the target must reach MIN_EFFECT (10%, an assumption); a staple that
               lowers the dose, or barely moves it, is not a candidate
4. sort        by ΔT/T on the target, descending (ties on the rounded value break by name)
5. label       confidence tier shown beside each row, not sorted on
6. report      anchor separation, PAM availability, provenance — raw, for the researcher
```

**Why a gate and not a weighted score.** A weighted sum would let a large predicted effect
outrank a safety concern — that is what weighted sums do. A gate cannot be traded against.
Deliverability stays out of the model entirely and is reported as facts for the researcher to
judge.

The ranking is deterministic and auditable end to end. Given the same inputs it produces the same
order, and every step of that order can be inspected.

---

## Real-data run: K562, β-globin locus (chr11:4.7–6.7 Mb)

Run on real inputs, none synthetic: K562 intact Hi-C (ENCODE ENCFF621AIY, range-read from a 33.8 GB
file, 5 kb, SCALE balancing), K562 DNase peaks and signal (ENCSR000EKS), K562 RNA-seq
(ENCSR000AEM), RefSeq genes, JASPAR CTCF motifs, ClinGen dosage genes. Code:
`pipeline/real_data.py`, `pipeline/real_run.py`, `pipeline/validate_expression.py`; outputs in
`data/derived/`. Every step was then re-derived by an independent verifier agent that tried to
break it; the defects they found are fixed or listed below.

**What worked**

| check | result |
|---|---|
| Hi-C is sane | symmetric, finite, distance-decay slope −1.05 (typical −1.0 to −1.3) |
| T tracks measured expression | Spearman **ρ = 0.67** over 76 genes; random-promoter null 95th pct 0.19, beaten in 2,000/2,000 |
| LCR is top ABC element for all five globin genes | yes, but it wins on *activity* (5.09 vs window median 0.72), not contact; HBB margin only 1.3×. A consistency check, not independent confirmation |
| LCR→HBB staple raises HBB | **+43.8%**, from a contact at only the 80th percentile, so it is a real forced contact. Sign agrees with forced-looping experiments (Deng 2012 *Cell*; 2014 *Nature*), which were not done in K562 |

**What did not work, stated plainly**

- **The search found nothing.** Under the strict step-3 rule only 4 of 81 accessible sites are
  permitted anchors (gene density: bodies, TSS padding and enhancers of expressed genes remove 75%
  of the window), so 604 of 608 candidate staples had no valid anchor. None of the 4 reaches the
  10% dose floor. Honest headline: **no staple in this window meets the dose goal under the strict
  rule.** This is a finding about step 3, not about the data.
- **The safety gate was never exercised.** No ClinGen-flagged gene lies in this window (of 76
  coding genes only DCHS1 appears in ClinGen at all, with insufficient evidence), and COSMIC and
  DepMap could not be obtained. Both gate tests are vacuous here, so `safe = True` carries no
  information. The gate's logic is tested on synthetic data and by property tests only.
- **Two of four confidence tests are constants** (abstention: UniversalEPI not run; provenance:
  measured Hi-C).
- **The off-target readout depends on σ; the target effect does not.** LCR→HBB is +39.4%, +43.8%,
  +42.5% at σ = 1.5, 3, 4.5 bins, but the neighbour HBD moves +2.0%, +11.9%, +12.7%.
- **Neighbours 5 kb apart cannot be told apart.** Staple the LCR to HBD and HBB moves +41%, 2.6×
  the intended effect. This is the ≤30 kb design rule seen on real data.
- **The yardstick makes the strongest contacts a no-op.** LCR→HBG2 is already the strongest
  contact at its separation (100th percentile), so the staple does nothing and is reported as such.
- **ABC is unaware of expression.** An unexpressed olfactory receptor (OR51V1, 0 TPM) scores
  T = 746, comparable to HBG2 (758; 8,201 TPM). T measures regional activity, not on/off state.
- **Low-coverage bins** (bin 102 at 5,210 kb, two bins upstream of HBB; bin 168) have row sums
  ~0.5× the median even after balancing, so contacts involving them are biased low. Not masked.

## Assumptions, stated plainly

These are not measurements. They are choices, and a reviewer should be able to find them in one
place.

1. **σ ≈ 15 kb** — the width of a staple's influence. Supported by CRISPRi crosstalk at 1–10 kb
   (`P = 1.8 × 10⁻⁵`), not measured for a staple.
2. **Row-sum conservation** — that a locus has a fixed contact budget. Physically motivated,
   not verified.
3. **τ = 10%** — the point at which a flagged gene has moved too far. Arbitrary starting value.
4. **MIN_EFFECT = 10%** — the dose change below which a staple is not a candidate. Set equal to τ
   for scale, with no independent justification.
5. **Staple strength = strongest real contact at the separation** — a describable yardstick, not a
   measurement, and it makes already-strong pairs a no-op (step 4c).
6. **Diagonal imputed from neighbouring bins** — a choice; zeroing is the alternative and biases T low.
7. **Contact substitution** — ENCODE-rE2G takes contact from Hi-C or the megamap; feeding it
   UniversalEPI predictions instead is an unvalidated swap.

---

## Limitations, for the future-work slide

- **No expression calibration.** `ΔT/T` has no units a biologist can act on. Fixing this means
  regressing it against measured CRISPRi effect sizes (Gasperini, Schraivogel).
- **Step 3 may be too strict at 5 kb.** In a gene-dense locus nearly every accessible site is an
  enhancer of an expressed gene or a promoter, so "accessible but not functional" is almost empty.
  Options, not yet decided: exon-level rather than whole-transcript gene bodies (two genes,
  MMP26 at 287 kb and OR51B5 at 165 kb, forbid ~465 kb alone), a finer resolution, or allowing
  anchors in the flanks of peaks. Dropping the body rule alone only lifts 4 anchors to 8.
- **Published forced looping puts the effector *at* the promoter,** which the TSS exclusion forbids
  by design. Validation therefore has to waive the exclusions; the pipeline's own recommendation
  and the published design are different experiments.
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
ClinGen dosage sensitivity · JASPAR (via UCSC `jaspar2022`) · UCSC Genome Browser API
(`ncbiRefSeqCurated`, `knownGene`) · ENCODE (Hi-C ENCFF621AIY / ENCSR479XDG, DNase ENCSR000EKS,
RNA-seq ENCSR000AEM) · Rao et al. 2014 (in situ Hi-C method) · Deng et al. 2012 *Cell*, 2014
*Nature* (forced looping) · hic-straw and pybigtools (range-reads) · NumPy. Not obtained:
COSMIC Cancer Gene Census (login), DepMap (verification page).

## Reproducing the numbers

```bash
uv run --with numpy python -m pipeline.simulate
uv run --with matplotlib --with numpy python animation/figures/method_figures.py

# real K562 data (hic-straw needs an older macOS SDK on recent Xcode command-line tools)
SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk \
  uv run --with hic-straw --with pybigtools --with numpy python -m pipeline.real_run
```

The demo locus (`pipeline/locus.py`) is **synthetic and marked `illustrative`** — element
positions, activities and gene positions are invented, and the contact matrix is generated from a
distance-decay power law with TAD blocks and CTCF loops. The arithmetic applied to it is the real
pipeline arithmetic. Only the inputs are made up.

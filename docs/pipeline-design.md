> **SUPERSEDED (2026-10-03) by [`docs/method.md`](method.md).** Kept as the record of how the method
> was arrived at, including the composition problem in section 1 and the reasoning that led to
> using ABC for magnitude. Do not implement from this file.

# Pipeline design — technical brainstorm

Written 2026-10-03 from the two paper briefs (`docs/briefs/`). Read `docs/targets.md` for the
target list this consumes. **Status: proposal, not agreed.** Section 1 contains a problem that
changes the architecture; settle it before writing pipeline code.

---

## 1. The problem nobody has stated yet: neither model can score a staple

Both papers were read with the pipeline in mind. The finding that matters:

| | UniversalEPI (Grover 2026) | ENCODE-rE2G (Gschwind 2026) |
|---|---|---|
| Question it answers | How often do these two loci **touch**? | Does this element **regulate** this gene? |
| Output | log Hi-C per ATAC-peak pair, **+ uncertainty** | **Calibrated probability** in [0,1] |
| Input that carries cell type | ATAC-seq | DNase-seq (ATAC weaker) |
| Max distance | 2 Mb (benchmarked to 1 Mb) | trained 1 kb–1 Mb |
| In-silico edits | **Yes** — CTCF inversion and a SNP, both validated | No |
| Effect size | No | **No** — score vs CRISPRi effect r = −0.45 only |

**The gap.** A staple forces two loci together *without changing the DNA sequence and without
changing chromatin accessibility*. UniversalEPI's only inputs are sequence and ATAC. So a staple is
**invisible to UniversalEPI's input layer** — there is no channel to express "I have physically
tethered A to B." You cannot run a staple through it as an edit the way the paper runs a CTCF
inversion.

Worse, UniversalEPI's authors say in terms that predicted strong contacts **are not strictly
functional**: elements can touch before any expression change. So contact alone is not an answer
either.

### The fix: compose the two models through the contact feature

ENCODE-rE2G **takes 3D contact frequency as an input feature** — it is, with distance, the most
informative feature group in the model. That is the seam.

```
                 sequence + ATAC
                        │
                 UniversalEPI            ← contact map, with uncertainty
                        │
          ┌─────────────┴─────────────┐
   reference contact            do(contact(A,B) := stapled)   ← the staple is an
          │                            │                        INTERVENTION ON THE
          └─────────────┬─────────────┘                         FEATURE, not on the input
                        │
                  ENCODE-rE2G           ← P(E regulates G), before and after
                        │
                   Δ probability  →  tier, not a number
```

The staple is modelled as a **do-operation on the contact feature feeding rE2G**, not as a sequence
edit. Each model is used only where it was validated: UniversalEPI for contact (including
counterfactual contact under genuine *sequence* edits, e.g. CTCF changes), rE2G for whether contact
translates into regulation.

### What this buys, and what it costs

Buys:
- A principled, auditable path from "staple here" to "expression probably moves".
- Natural handling of **down**-regulation: lower the contact feature (insulate) and rE2G's score
  falls. We do **not** need silencers — which is essential, because rE2G models activation only.
- UniversalEPI's uncertainty propagates: when the ensemble's intervals overlap, the
  maximum-confidence fold-change is **exactly zero**, which is a free, principled **abstain**.

Costs, all of which go in the limitations panel:
- rE2G normally takes contact from the ENCODE Hi-C Megamap or cell-type Hi-C. **Substituting
  UniversalEPI predictions is an unvalidated swap.** Mitigation: on demo loci, compare UniversalEPI's
  contact against the real Hi-C rE2G would have used, and publish the agreement. If they disagree,
  say so.
- **What value does a staple set the contact to?** We don't know. Nobody has measured the contact
  frequency a dCas staple produces. Proposal: don't pick one. Sweep the contact feature across its
  observed genome-wide range and report the **dose–response curve of rE2G score vs imposed contact**.
  A curve that is flat means the staple cannot help; a curve that is steep means it can. This turns
  an unknown parameter into the actual output, which is both more honest and more useful.
- rE2G gives **probability of an effect, not size of effect**. Δscore is not Δexpression. Tiers only.

**Open question for the team:** is anyone willing to defend a different composition? The alternative
is to treat the whole thing as retrieval — shortlist enhancers by rE2G score and rank by contact —
with no counterfactual at all. That is weaker but much safer. Decide explicitly.

---

## 2. Walking the researcher's workflow, step by step

The user's spec, with what each step actually costs.

### 2.1 "Which chromosome is the gene on?"
Trivial, and already solved: ClinGen's TSV carries GRCh38 coordinates per gene (`docs/targets.md`).
Use GENCODE v29 TSS to match rE2G's gene model — rE2G uses one 500 bp TSS-centred promoter per gene
symbol, from the RefSeq TSS with the most coding isoforms. **Match their definition exactly** or
every downstream join is subtly wrong.

### 2.2 "Which TAD is the gene in?" — harder than it looks

**UniversalEPI has no TAD or boundary output head.** The paper uses TADs only as an analysis
construct. So we must derive them. Three options:

1. **Published TAD calls** for the chosen cell types. Cheapest, but ties us to cell types with calls.
2. **Insulation score computed on the UniversalEPI predicted map.** Keeps everything in one model and
   works for any cell type with ATAC. Not something the paper validates — we would be first.
3. **Skip TADs; use the contact map directly.** Arguably the most honest: a TAD is a summary of the
   contact map, and the staple decision only needs the map. The TAD is for the *human* reading the
   dashboard.

Recommendation: **(3) for the scoring, (1) for the picture.** Show the researcher a TAD because it
is how they think, but never let a TAD call drive the ranking. Note the reason out loud: TAD calls
are method-dependent and a boundary is not a wall — the deck's slide 2 already makes this point.

### 2.3 "What CREs regulate the gene?"
Straight lookup in the encyclopedia: **92,176,227 E–G links, 1,458 biosamples, 369 cell types,
hg38**. Pull, don't re-run.

Priors worth hard-coding as sanity checks — if our locus is wildly outside these, something broke:
- median E–G distance ≈ **37 kb**; **74.6%** of links under 100 kb
- **~3 enhancers per gene (median), 5.91 mean**
- **~1.6 genes per enhancer (median 1, max 34)** ← this is the enhancer-theft risk in 2.5

And the honesty number for the UI: at the authors' 70%-recall threshold, precision is **53–58%**.
**Roughly one in two shortlisted links is a false positive.** That belongs on the dashboard in words,
not buried.

### 2.4 "Which two regions to staple?" — the anchor-placement constraint

The user's constraint — *the staple must not block protein recruitment to the enhancer or promoter* —
is the most concrete and implementable part of the whole spec. It is a geometry problem on
coordinates we already have.

Define an **anchor-exclusion set**, and require both anchors to miss all of it:

| Exclude | Why | Where the coordinates come from |
|---|---|---|
| The 500 bp rE2G element itself | dCas over the enhancer core blocks TF binding — this is literally how CRISPRi works | rE2G element BED |
| Promoter, ±1 kb of TSS | same, for Pol II / GTFs | GENCODE v29 |
| Target gene body | rE2G excludes these too; also risks transcriptional interference | GENCODE v29 |
| **CTCF motif occurrences** | moving or occluding a boundary anchor has effects far beyond our locus | JASPAR **MA0139.1** PSSM scan — the same motif UniversalEPI uses, with orientation |
| Other elements within ~1 kb | CRISPRi at one element measurably reduces H3K27ac at neighbours 1–10 kb away (P = 1.8×10⁻⁵) | rE2G element BED |

Then **prefer** anchors in accessible-but-non-regulatory flanking DNA: present in the ATAC peak set
(so the staple can reach it in that cell type) but carrying no rE2G link and no CTCF motif. There is
a real tension here — accessible regions tend to be functional — so expect this filter to bite, and
report how many candidate anchors survive it. If zero survive, **that is a result**: this locus is
not stapleable, say so rather than relaxing the filter silently.

Sub-question we should not fudge: a dCas9 footprint is ~20 bp of protospacer plus a large protein.
"Does not interfere" is a claim about a few hundred bp, not a few bp. Keep the exclusion zones
generous and state the assumed footprint as a parameter.

### 2.5 "Secondary risks after stapling" — score the window, not the pair

**Revised 2026-10-03 after team input.** The earlier draft listed separate risk flags computed
around the target pair. Better framing: **the unit of analysis is the window.** Compute Δ rE2G score
for *every* (element, gene) pair in the window, before and after the intervention, and read the
whole table. The target gene is just the row you hoped would move.

This is not extra work — it is the natural shape of the data. UniversalEPI already emits all
pairwise contacts among the central 200 peaks, and rE2G is a logistic regression, so scoring every
gene in the window costs almost nothing beyond scoring one.

One sweep subsumes three things that would otherwise be separate checks:

- **enhancer theft** — an enhancer serves 1.57 genes on average, up to 34; its other targets lose out
- **collateral gain** — the new contact reaches a promoter we did not intend
- **boundary collateral** — genes affected by a disturbed CTCF boundary

Two flags survive as genuinely separate, because they are not visible as a row moving:

- **Super-additive overshoot.** Enhancer pairs 3–89 kb apart with high contact combine
  **16–39% above additive**; pairs >100 kb apart only 2–10%. Stapling an enhancer *into* an existing
  cluster can overshoot. For a narrow-window gene like MECP2 that is the failure mode, not a
  footnote. rE2G's additive-ish score will not show this — it has to be flagged from the geometry.
- **Boundary disturbance as out-of-scope.** If an anchor sits near a CTCF boundary, the consequences
  extend outside the modelled window. Flag it as *not covered* rather than scoring it.

**Cell type is a precondition, not a detail.** "The genes nearby" is only defined once a biosample is
fixed — rE2G predictions are per-biosample across 1,458 of them. Pick cell types first.

#### Annotating the rows that moved

A row changing is not a risk until the gene is one you must not perturb. Deterministic sources:
curated oncogene lists; **ClinGen triplosensitivity** (pleasing symmetry — the table defining our
targets also defines our risks: raising a TS gene is dangerous by our own definition); essentiality.

**Coverage limit to state in the UI, not hide:** triplosensitivity is curated for 3 genes and 20
regions. There is no dense resource for "is it dangerous to overexpress this gene." The risk layer
flags known hazards and **cannot certify safety**. An unflagged gene is *unexamined*, not *cleared*.

#### Where the LLM is allowed to act

Literature review over flagged genes is worth doing, but it sits in the causal path of the ranking,
which is where hallucination does real damage — and `CLAUDE.md` already says ranking is
deterministic. One rule resolves it:

> **Asymmetric trust.** The model may *raise* a concern level, with a citation. It may never *lower*
> one. "No concern found" carries zero weight.

A search that finds nothing is evidence of a search, not of safety. Escalation costs little when
wrong and catches what no database holds; de-escalation means a fabricated all-clear can silently
promote a target. Grades move one way, every escalation carries a checkable citation, and the
deterministic floor is never raised by prose.

### 2.8 The down direction needs a mechanism (open gap)

Roughly half the target list is *too much* of a gene, and **a staple only ever creates a contact**.
Neither the original plan nor my first draft said how that becomes less expression. Three options,
not equally plausible:

1. **Competitive sequestration** — staple the enhancer to a decoy so it is occupied elsewhere.
   Expressible as lowering the E–P contact feature, but "an enhancer can be monopolised" is an
   assumption, not a result.
2. **Relocation to a repressive compartment** — CRISPR-GO already moves a locus to a *nuclear body*
   rather than to another locus (deck slide 5). This is the most direct answer and it is already on
   our own method slide; it just was never connected to the too-much targets.
3. **Engineered insulating loop** — force an extrusion barrier between enhancer and promoter. Most
   speculative, furthest-reaching off-target consequences.

Decide early: it determines whether the demo can claim both directions. MECP2 is attractive
*because* it shows both — but only if both have a mechanism behind them.

### 2.6 Aggregation into something a human can read

Per `CLAUDE.md`: deterministic and auditable. Proposal — **no single number**. Four independent axes,
each a letter, shown as a row:

```
  TRACTABLE   CONFIDENT   SAFE   DELIVERABLE
      A           C         B         A
```

- **Tractable** — is there a staple that passes the anchor-exclusion filter and moves the rE2G score?
- **Confident** — UniversalEPI uncertainty; **C if the maximum-confidence FC is zero** (model abstains)
- **Safe** — worst of the four risk flags
- **Deliverable** — cell type reachable, reversibility available (§3)

Why not one score: the four axes fail for *different reasons* and a researcher needs to know which.
A target that is `A/A/A/C` needs a delivery chemist; one that is `A/C/A/A` needs an experiment. A
single 73/100 hides both. This also makes the thing honest by construction — we cannot average away
a C in Safe.

LLM involvement stays where `CLAUDE.md` puts it: writing the rationale sentence from the structured
record. It never touches a letter.

### 2.7 Choosing the stapling method

A lookup table, not a model. The deck (slide 5–6) already enumerates six stapling systems and four
reversibility mechanisms. The decision inputs are: cell type and whether it is dividing, in vivo vs
ex vivo, required speed of reversal, and whether permanence is acceptable. Ex vivo haematopoietic
targets (RUNX1, BCL11A) are the easy case; CNS targets (MECP2, SCN1A) are the hard one and the
honest answer there may be "delivery is the blocker, not the design."

---

## 3. Validation: how we show we know when we're wrong

This is the Track 2 pitch and it should be built first, not last.

1. **Positive control — ZRS → SHH.** A ~800 bp enhancer ~1 Mb from its target, where human genetics
   gives the ground truth, and where the enhancer itself is ClinGen triplosensitive (see
   `docs/targets.md`). **If the pipeline cannot recover ZRS → SHH, it is wrong.** Run it, show the
   result, pass or fail.
2. **Negative control.** Keep one target the pipeline should reject (`TBX5` is a candidate — a
   developmental window likely closed post-natally). A pipeline that scores everything is a pipeline
   that means nothing.
3. **Known failure modes, declared up front.** From the papers, not invented: very long-range links
   are missed (all models failed three real MYC enhancers at ~2 Mb); ubiquitously expressed genes are
   poorly predicted (1.5% vs 7.6% positive rate); small effects (<15%) are systematically
   under-ascertained; everything supervised is K562-trained.
4. **The abstain path.** Maximum-confidence FC = 0 → the pipeline says "I don't know" and shows why.
   Demo this on a real locus. An agent that declines is the whole point of the track.

---

## 4. Compute, and whether Modal is needed

| | |
|---|---|
| UniversalEPI inference | **8.4 min single model / 84.5 min 10-model ensemble** for three whole test chromosomes. One 4 Mb locus is sub-second on GPU (my extrapolation, not a paper number). |
| Model size | 2.5M parameters |
| rE2G | logistic regression — negligible |

**So the honest answer is that we do not need a GPU for a 15-gene demo.** Precompute, cache small
JSON, done.

Where Modal genuinely earns its place is the **contact dose–response sweep** in §1: 10-model
ensemble × N contact values × N candidate anchor pairs × N loci is embarrassingly parallel and is the
one thing that is actually expensive. If we want the Modal prize, that is the defensible use — not
wrapping a sub-second forward pass in a GPU container to look impressive. Say this plainly in the
pitch; judges can tell.

---

## 5. Data to pull first (blocks everything else)

1. rE2G predictions for chosen biosamples — ENCODE portal; **accessions are in Supplementary Table 12,
   which is not in our PDF**. Someone must fetch the Nature supplementary.
2. UniversalEPI precomputed predictions — UCSC track hub for 157 ENCODE ATAC datasets (116 cell lines
   + 41 primary). **The paper gives only the generic hgHubConnect URL**; find the real hub via the
   repo wiki or Zenodo (`10.5281/zenodo.14622040`).
3. GENCODE v29 TSS, matching rE2G's promoter definition.
4. JASPAR MA0139.1 for the CTCF scan.
5. An oncogene list for risk flag 1.

**Licence check before we build:** the UniversalEPI *article* is CC BY-NC 4.0. The **repo's software
licence is not stated in the paper** — check it directly. Non-commercial terms matter even for a
hackathon credit slide.

---

## 6. What I would cut

The plan targets 10–15 genes and 2–3 cell types. Given §1 is unresolved, I would rather ship
**one locus end to end with the abstain path and the ZRS control working** than fifteen loci with a
score we cannot defend. The demo shape in `background.md` already says "one hero gene end to end" —
that instinct was right; the gene list is breadth for the dashboard, not for the science.

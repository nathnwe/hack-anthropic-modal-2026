# Pipeline results for the website

**Status:** 2026-10-04. Real K562 data (Hi-C, DNase, RNA-seq, gene models, ClinGen). Everything here is an
*in silico prediction, unvalidated*. Method: [method.md](method.md). Data: `data/derived/website_results.json`
(also `demo_*.json` per target, `screen/` per screened gene, `screen_tables.md`).

## 1. What to show

For each target gene the pipeline returns a ranked **staple design** (an enhancer-side anchor pulled to a
promoter-side or decoy anchor), the **predicted change in the target's regulatory input (ΔT/T)**, and a
**risk tier** from every other dosage-sensitive gene in the 2 Mb window.

| gene | goal | disease | best design (enhancer → anchor, separation) | ΔT/T at yardstick | ΔT/T at p99 | risk | confidence |
|---|---|---|---|---|---|---|---|
| **NF1** | raise | Neurofibromatosis type 1 | 31,555 → 31,080 kb (470 kb) | +32.2% | **+10.9%** | CLEAR | tentative |
| **NSD1** | raise | Sotos syndrome | 177,365 → 177,125 kb (220 kb) | +40.0% | −0.5% | LOW (1.9%) | tentative |
| **RUNX1** | raise | Familial platelet disorder, AML predisposition | 34,615 → 34,875 kb (255 kb) | +15.2% | +4.2% | CLEAR | **confident** |
| **MECP2** | raise | Rett syndrome (loss) | 154,000 → 154,100 kb (105 kb) | +16.4% | +6.4% | LOW (0.4%) | **confident** |
| **LDLR** | raise | Familial hypercholesterolaemia | 10,460 → 11,080 kb (610 kb) | +10.1% | +6.8% | LOW (0.0%) | tentative |
| **MECP2** | **lower** | MECP2 duplication syndrome (gain) | 154,120 → decoy 154,300 kb (185 kb) | −11.0% | −0.7% | **ELEVATED** (2.2%, FLNA) | confident |
| **MYC** | **lower** | Oncogene, overexpressed in cancer | 127,740 → decoy 127,270 kb (475 kb) | −10.4% | −1.8% | CLEAR | confident |

(MECP2-raise lists the top *confident* design; the top by size is +19.6% at 30 kb, p99 +3.0%, tentative.)

**How to read the two ΔT/T columns (this matters).** The staple's strength is unknown, so we simulate it from
weak to strong. *Yardstick* is the strongest natural contact at that separation anywhere in the window: the
ceiling. *p99* is the 99th percentile of real contacts at that separation: a conservative figure. The headline
sits at the top of the range, not the middle. Only NF1 stays above the 10% goal at p99. **Show both, or show the
conservative one.** NSD1's +40% is almost entirely a yardstick effect (it is ≈0 at every weaker strength).

## 2. The five points the page should make

1. **Raising dose is predicted to work; lowering it is much weaker.** Raising gives +10 to +40% at the
   yardstick, up to +11% at p99. Lowering (pulling an enhancer away to a decoy) gives about −10 to −11% at the
   yardstick and under −2% at p99. The lowering mechanism rests only on the contact-conservation assumption.
2. **The safety gate does real work.** At MECP2 the largest raising design (+43.9%) is **rejected** because it
   would cut the neighbouring haploinsufficient gene **AVPR2 by 20.4%**. A tool that only maximised the target
   would have ranked it first.
3. **A null result is a result.** RUNX1, NF1 and MYC have no other ClinGen dosage-sensitive gene in the
   window, so they score CLEAR: *to our knowledge this region is safe*. That is the limit of one curated list, not proof.
4. **Unflagged collateral is shown, not gated.** The MECP2-lower design also doubles TEX28 and TKTL1 (both below 1 TPM
   in K562, neither on ClinGen), so it is not blocked. The site should show the collateral list.
5. **Confidence is mostly tentative.** Only distance (≤500 kb) and sign-stability across strengths are live tests;
   UniversalEPI abstention is not run and provenance is constant (measured Hi-C).

## 3. Per-gene notes

- **NF1** (chr17:30.2–32.2 Mb, 14 designs scored, 4 reach the goal). Best at p99 of the set. Four designs
  +12 to +32%, all CLEAR, all tentative because they are not sign-stable at the weakest strengths.
- **NSD1** (chr5:176.2–178.2 Mb). Neighbour DDX41 is haploinsufficient; the top designs move it by ≤1.9% (LOW).
  The +40% headline is not supported by weaker strengths; present NSD1 as a strong-yardstick, weak-conservative case.
- **RUNX1** (chr21:33.9–35.9 Mb). One design, sign-stable at every strength tested: the cleanest case.
- **MECP2** (chrX:153.06–155.06 Mb, 9 ClinGen-flagged genes: the densest safety test of the set). Raise: 2 designs
  survive, 1 rejected by the gate (AVPR2). Lower: 130 decoys tried, 1 reaches −10%, with ELEVATED risk (FLNA +2.2%).
  Both Rett (raise) and duplication (lower) are pathogenic directions for this one gene, so a design needs a
  dose *window*, which ΔT/T alone does not provide. K562 expresses MECP2 at 24 TPM; K562 is not neural tissue.
- **LDLR** (chr19:10.1–12.1 Mb). One design, 610 kb (beyond the 500 kb confident limit), barely above the floor. The weak case.
- **MYC** (chr8:127–129 Mb). The window holds only two coding genes (MYC sits in a gene desert), so collateral is
  near-empty and the CLEAR tier is trivially easy to reach. MYC's well-known distal enhancers (K562 enhancers were mapped by
  Fulco et al. 2016, PMID 27708057) can lie beyond a 2 Mb window; we did not compare our design with their coordinates.
  One design of 197 reaches the −10% floor (−10.4%, p99 −1.8%).

## 4. Larger gene set: what is finished, and when the rest finishes

**Design.** Up screen: the 270 ClinGen haploinsufficiency-3 genes expressed in K562 (TPM ≥ 1) out of 422, one
2 Mb window each, centred on the gene (`demo_window_scan.json`), run in a fixed random order (seed 0). Down screen: the
ClinGen triplosensitivity-3 genes *LMNB1, PLP1, TPSAB1* plus MECP2 and MYC.

**Timeline.**

| step | status |
|---|---|
| Demo set (7 targets, 5 genes, both directions) | **done** |
| Down screen (5 genes) | **done** |
| Up screen, 5 of 270 genes | **done, then stopped on request** |
| Up screen, remaining 265 genes | **not run; estimated ≈ 5 h** |

Why about 5 h: each new window downloads Hi-C (about 90 s, range reads from a 33.8 GB ENCODE file) and CTCF motif hits
(about 100 s from UCSC); alone that is ~4 min per window, and ~7–8 min per window per shard with 6 shards in parallel,
i.e. about 50 windows per hour. The compute after the download takes seconds. To finish, run
`python -m pipeline.screen --set up --expressed --shard I/8` for I in 0..7 (resumable: finished genes are skipped), then
`pipeline.screen_report` and `pipeline.export_web`.

**What the 10 finished records show** (too few to generalise; treat as a smoke test, not a result):

| gene | goal | K562 TPM | outcome |
|---|---|---|---|
| AAGAB | raise | 39 | +20.1% (p99 +7.8%), LOW |
| PQBP1 | raise | 127 | +11.0% (p99 +1.6%), **ELEVATED** 7.5%; 4 other flagged genes nearby |
| ADNP | raise | 16 | +10.5% (p99 +2.8%), LOW |
| LDLR | raise | 44 | +10.1%, same as the demo (same window) |
| MSH6 | raise | 58 | no design reaches 10% (best +2.5%) |
| LMNB1 | lower | 74 | none reaches −10% (best −6.9%, CLEAR) |
| PLP1 | lower | 0 | not expressed in K562: meaningless |
| TPSAB1 | lower | 0 | not expressed in K562: meaningless |

(MECP2 and MYC, the other two down-screen records, are in the table in §1.)

## 5. Caveats that must appear on the site

1. In silico, unvalidated. No wet-lab result supports any number here. (The only external consistency check: the model
   predicts a rise for the LCR → *HBB* staple, matching the sign of forced-looping experiments of Deng et al. 2012 and 2014,
   in other cells.)
2. **One cell line.** Hi-C and expression are K562 (leukaemia). Rett, Sotos, FH and NF1 act in other tissues; a
   gene K562 does not express (PLP1, TPSAB1) says nothing.
3. **Demo windows were chosen after looking** at the 422-gene scan. They are not a random sample.
4. **Risk = ClinGen only.** COSMIC (oncogenes, tumour suppressors) and DepMap (essential genes) could not be obtained.
   CLEAR means "no ClinGen dosage-sensitive gene nearby".
5. **ΔT/T is a regulatory-input change, not an expression change.** Calibration to expression is future work.
6. **Down-regulation is the weakest mechanism** (§2.1); rows from `mode: "down"` should carry that note.
7. Per CLAUDE.md, show tiers rather than false precision: the ΔT/T numbers are for a results table, not a headline.

## 6. File dictionary (`data/derived/website_results.json`)

`demo[]`: `gene, disease, goal (raise|lower), window, target_tpm_k562, flagged_genes_in_window,
n_designs_scored, n_meeting_dose_goal, designs[], rejected_by_safety_gate[]`. A design has `rank,
enhancer_kb, anchor_a_kb, anchor_b_kb, separation_kb, delta_target` (yardstick), `delta_at_p99`, `sweep`
(p50, p75, p90, p95, p99), `sign_stable, confidence (CONFIDENT|TENTATIVE), risk{tier CLEAR|LOW|ELEVATED|BLOCKED, score,
worst_gene, reason}, unflagged_collateral_over_tau[]`, and `q`. `screen.records[]` mirrors `screen/*.json`.
It is **not** `contracts/schema.json` (qualitative tiers, no BLOCKED); mapping it to the contract needs all three
operators' agreement. Suggested map: tier `CLEAR/LOW` → `low`, `ELEVATED` → `medium`, `BLOCKED` → not shown as a design.

## 7. Changes to the pipeline in this round

- **Down-regulation (decoy staple)** added: `search_down` in `pipeline/real_run.py`, documented as step 4d in `method.md`.
- `delta_at_p99` is now stored with every design.
- "Rejected by the gate" now counts only designs that *would have met the dose goal*. This changed the earlier
  claim that NSD1 had a gated design (its target moved −4.8%, below the goal); the real gate example is MECP2 (AVPR2).
- `best_below_floor` now sorts toward the goal for lowering designs (it previously sorted the wrong way).
- A first attempt to require decoys ≥30 kb from every promoter was reverted: it emptied the gene-dense MECP2 window
  and was tuned to a single anecdote. Unflagged collateral is reported instead.
- New: `pipeline/screen.py`, `screen_report.py`, `export_web.py`.

## Sources

ENCODE K562 intact Hi-C ENCFF621AIY, DNase ENCSR000EKS, RNA-seq ENCSR000AEM; UCSC RefSeq and JASPAR CTCF; ClinGen
Dosage Sensitivity Map (2026-10-03); ABC model, Fulco et al. 2019 *Nat Genet*; CRISPRi enhancer mapping in K562, Fulco et al.
2016 *Science* (PMID 27708057); MECP2 duplication syndrome, Van Esch et al. 2005 (PMID 16080119).

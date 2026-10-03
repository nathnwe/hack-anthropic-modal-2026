---
title: "An encyclopedia of human enhancer–gene regulatory interactions"
authors: "Gschwind, A. R.; Mualim, K. S.; Karbalayghareh, A.; et al. (ENCODE Consortium); corresponding: Jesse M. Engreitz"
journal: "Nature"
year: 2026
volume_pages: "657:179–190 (3 September 2026)"
doi: "10.1038/s41586-026-10781-4"
local_pdf: "/Users/nathanewer/nwe_src/hack-anthropic-modal-2026/docs/Gschwind et al. - 2026 - An encyclopedia of human enhancer–gene regulatory interactions.pdf"
date_read: 2026-10-03
note: "The Supplementary Information (Supplementary Notes and Supplementary Tables 1–19) was NOT present in the local PDF; numbers sourced only from those are flagged as unavailable throughout."
---

## ENCODE-rE2G (Gschwind et al., *Nature* 657:179–190, 3 Sept 2026) — engineering summary

Source read: `/Users/nathanewer/nwe_src/hack-anthropic-modal-2026/docs/Gschwind et al. - 2026 - An encyclopedia of human enhancer–gene regulatory interactions.pdf` (31 pages = main text, Methods, Extended Data Figs 1–9, Nature Reporting Summary). **The Supplementary Information / Supplementary Notes / Supplementary Tables are NOT in this PDF** — many precise numbers (per-model score thresholds, the full 12/45-feature lists, Supplementary Tables 1–19) live there and are flagged below as "in Supplementary, not in this PDF".

DOI: https://doi.org/10.1038/s41586-026-10781-4

---

### 1. What ENCODE-rE2G is

- **Model class**: supervised **logistic regression** classifiers (plain, unregularized after variance-stabilization), trained to predict whether an element–gene (E–G) pair is a true regulatory interaction. Deliberately simple/regularized so the output is a **calibrated probability** ("This approach ensures regularization while maintaining that the model outputs calibrated probability estimates"). Features are variance-stabilized (`log(|x| + 0.01)`) before fitting. Training uses **hold-one-chromosome-out cross-validation** (each of chr1–22 and chrX held out in turn); for application to new cell types a model trained on all 23 chromosomes is applied.
- **Candidate elements / pair definition**: elements called with MACS2 (`--shift -75 --extsize 150 --nomodel`) on pooled DNase-seq or ATAC-seq replicates; blacklist removed; **top 150,000 peaks by read count**, resized to **500 bp centred on peak summits**. 500-bp TSS-centred regions added for all genes; overlapping regions merged. Elements classified promoter / genic / intergenic. **All pairs within 5 Mb** are candidate E–G pairs.
- **Feature groups** (12 initial features for the core model; 8 selected for the final model; feature table in Supplementary Table 2, not in this PDF):
  1. **Chromatin state (from chromatin assays)** — quantitative DNase-seq signal at the element and at the promoter; where available H3K27ac ChIP-seq. Read counts over candidate elements, gene bodies and promoters computed with bedtools coverage, then **quantile-normalized to a K562 reference**. Cell-type-specific enhancer activity = **geometric mean of normalized DNase-seq and H3K27ac signals** (if H3K27ac available).
  2. **3D contact frequency (from Hi-C)** — either cell-type-specific ENCODE Hi-C, or a **cell-type-averaged "Megamap" aggregating 65 ENCODE Hi-C datasets** (ENCODE accession **ENCSR492BKE**) to capture cell-type-invariant genome organization. Hi-C processing: SCALE normalization, power-law interpolation of low-coverage bins, doubly stochastic normalization, diagonal replacement, pseudocount adjustment. K562 Hi-C loop calls from the ENCODE portal.
  3. **ABC-derived** — the **ABC score itself** (activity × contact), plus ABC numerator and ABC denominator as separate features. The ABC score is the single most informative individual feature (Fig. 4b,c).
  4. **Genomic position / genomic features** — distance to TSS, **number of intervening TSSs between E and P**, number of intervening candidate elements, local element density within 5 kb, enhancer complexity score.
  5. **Promoter class** — promoter-class annotation (Bergman et al. 2022 sequence-based promoter compatibility classes) and **ubiquitous/housekeeping gene expression** annotation (whether the gene is expressed uniformly across cell types).
  6. **Nearby-enhancer / enhancer–enhancer features** — **sum of activity of all other elements within 5 kb of the perturbed element** (an enhancer-synergy feature). Ablating this feature alone costs 0.01 precision at 70% recall (P_bootstrap = 0.04).
  - Several features also enter as **second-order / squared terms** (e.g. `(contact frequency)²`, `(DNase signal at E)²`) to capture non-linearity.
- **Feature importance ranking** (Fig. 4b,c): most informative groups = **3D contact + distance**, then **enhancer activity**; most informative single feature = **ABC score**. Sequential forward selection order (Fig. 4c): ABC score → ubiquitous expression → #TSSs between E and P → DNase signal at P → #elements between E and P → distance to TSS → (contact frequency)² → DNase signal at E → contact frequency → activity of elements within 5 kb of E → (DNase signal at E)² → #elements within 5 kb of E.

**Training ground truth (CRISPR)**:
- **10,356 tested element–gene pairs in K562**, containing **471 positives**. Assembled by re-analysing and harmonizing three published K562 CRISPR enhancer screens: **Nasser et al. 2021 (mostly CRISPRi-FlowFISH)**, **Gasperini et al. 2019 (Perturb-seq)**, **Schraivogel et al. 2020 (TAP-seq)**.
- Composition of the 471 positives: **462 CRISPRi, 9 CRISPR deletions**.
- The full curated set across all sources: **471 unique positive element–gene pairs** where perturbation gave a significant decrease in expression (**−1 to −93% effects**) and **9,885 negative pairs** where no significant reduction was seen despite ≥80% power to detect **>15–25%** effects.
- Filters: pairs at **1 kb–1 Mb from TSS**, excluding elements overlapping GENCODE v29 promoters or target gene bodies. Negatives require **≥80% power to detect a 15% decrease** (Nasser subset filtered for ≥80% power at 25%). Duplicate pairs across datasets resolved by prioritizing significant results, matching significance status, highest statistical power at 25% effect size.
- Differential expression recomputed with **MAST v1.20.0** (Perturb-seq/TAP-seq), BH FDR < 5%; power by simulation with DESeq2 v1.34.0 negative-binomial parameters, knockdown efficiency a = 0.13 from FlowFISH data, 20 replicates.
- Pairs are weighted/annotated by **P_direct**, the probability the effect is direct (cis) rather than indirect/trans: `P_direct(d) = r_direct(d) / (r_direct(d) + r_indirect)`, where r_indirect is estimated empirically by testing each distal perturbation against 100 random genes on other chromosomes, and r_direct(d) is fitted as a power law of distance to TSS in 50-kb bins within 1 Mb.

**Held-out CRISPR test set**: 8 further screens in **5 cell types** — **K562, WTC11, HCT116, GM12878, Jurkat** — giving **4,378 total E–G pairs and 190 positives** (one passage says 4,405 pairs / 190 positives; both numbers appear, 4,378 is used in the Fig. 2d legend; **157.39 weighted positives**). Processed with SCEPTRE v0.10.0, 10% FDR.

**Genome build: GRCh38 / hg38 throughout** ("All coordinates are reported in GRCh38 unless otherwise specified"). Gene set: **20,678 gene promoters (one per gene symbol)**, 500-bp regions centred on the RefSeq TSS with the largest number of coding isoforms, GENCODE v29, protein_coding / processed_transcript / lincRNA.

---

### 2. Model variants, performance, minimum inputs

| Variant | Inputs | Features |
|---|---|---|
| **ENCODE-rE2G** (core) | **Cell-type-specific DNase-seq only**, plus the cell-type-averaged ENCODE Hi-C Megamap (not cell-type specific) | 8 selected from 12 |
| **ENCODE-rE2G^Extended** | DNase-seq + histone ChIP-seq + cell-type-specific ENCODE Hi-C + ChIA-PET, plus features imported from EpiMap and GraphReg | **45 features** (available in K562 and GM12878 only) |

**Minimum input to run on a new cell type: a single DNase-seq experiment** (the core model is "computable from only cell-type-specific DNase-seq data"), combined with the reference averaged ENCODE Hi-C Megamap. ATAC-seq works as an alternative ("ATAC-only" variant in Fig. 4f,g, but weaker). The authors released **pre-trained ENCODE-rE2G models for eight different combinations of input datasets** so you can match whatever assays a cell type has.

**Performance on the combined K562 CRISPR training benchmark (n = 10,356 pairs, 471 positives; threshold set at 70% recall):**

| Model | AUPRC | Precision at 70% recall | Pearson r vs CRISPRi effect size |
|---|---|---|---|
| ENCODE-rE2G (DNase-only) | **0.66** | **58%** | **−0.45** |
| ABC^(A=DNase, C=Average ENCODE Hi-C) | 0.56 | 48% | −0.40 |
| ENCODE-rE2G^Extended | **0.74** | **73%** | **−0.48** |
| ABC^(A=DNase×H3K27ac, C=ENCODE Hi-C) | 0.61 | 56% | −0.42 |

P_bootstrap = 0.0001 for ENCODE-rE2G vs ABC on AUPRC and precision; P = 1.44×10⁻⁵ and 3.29×10⁻⁷ for the r comparisons.

**Held-out (5 cell types, 190 positives)**: ENCODE-rE2G attains higher weighted AUPRC (P_bootstrap = 0.0003) and precision (P_bootstrap = 0.0001) than ABC and all other models, and **maintains precision 53.2% vs 54.8% on training data** — i.e. generalizes to new cell types with essentially no precision loss. Weighted AUPRC ≈ 0.5 for ENCODE-rE2G on held-out (Fig. 2d), vs ~0.2 for the distance/correlation baselines (P values from 3×10⁻⁴ down to 1×10⁻⁴).

**Assay-swap experiments (Fig. 4d,e,g)** — useful for deciding what to buy for a new cell type:
- Of **523 ENCODE one-dimensional chromatin experiments** tested as the activity term, **DNase-seq was the single best-performing assay**. ChIP-seq for cofactors **EP300, NCOA1, NCOR1** and several TFs achieved similar performance. **H3K27ac was the best histone ChIP-seq assay** (ranked 33rd overall).
- Of **6 methods for estimating 3D contact**, **cell-type-specific ENCODE Hi-C was best** (precision at 70% recall = **53.3%**; Hi-C depth = **23 billion reads, 5-kb resolution**), beating the cross-cell-type average Hi-C (**104 billion reads**) at **48%** (P = 7.62×10⁻⁸) and an inverse-distance function at **44.3%** (P = 6.35×10⁻¹²).
- Marginal utility of adding assays to a DNase-seq base: adding cell-type-specific H3K27ac and/or Hi-C both improve. **Models using DNase-seq beat models using ATAC-seq**, and adding H3K27ac or Hi-C gives a larger improvement for ATAC-based models (i.e. DNase-seq is itself a better enhancer-activity proxy). Pairwise P values in Fig. 4g: 8×10⁻⁴, 2×10⁻³, 6×10⁻¹, 1.19×10⁻³.
- Adding cell-type-specific Hi-C to the DNase-only model recovered **20 regulatory interactions** that the DNase-only model missed, all with **1.56–10.35-fold stronger 3D contact than a power law predicts** and **1.13–4.38-fold stronger contact than the tissue-averaged Hi-C map**. Adding K562 H3K27ac recovered **26**; adding both recovered **27**.

---

### 3. Output per E–G pair

- **Score = the logistic-regression predicted probability of a regulatory effect**, i.e. a calibrated probability in **[0, 1]**. "Each element–gene pair is annotated with a score, corresponding to the probability of a regulatory effect from the logistic regression classifier."
- **Recommended threshold: the score achieving 70% recall on the K562 CRISPR training dataset** (the same convention used by ABC/Nasser 2021). The **numeric score cutoff is reported in Supplementary Table 1/3, which is not in this PDF** — do not hard-code a number; take it from the released model/pipeline. Classification thresholds for *all* benchmarked models were set at 70% recall on the combined CRISPR data.
- **Precision at that threshold**: **58%** for ENCODE-rE2G (DNase-only) and **73%** for ENCODE-rE2G^Extended on the K562 training benchmark; **53.2%** for ENCODE-rE2G on held-out CRISPR data. ABC baselines: 48% / 56%.
- Links are also used continuously — predictor score correlates with CRISPRi effect size at **r = −0.45** (DNase-only) / **−0.48** (Extended).
- For reference, ABC score thresholds in this paper: **0.018** for ABC^(A=DNase, C=Average ENCODE Hi-C) across all 1,458 biosamples; **0.027** for ABC^(A=DNase×H3K27ac, C=ENCODE Hi-C) in K562 and GM12878.

---

### 4. The encyclopedia itself

- **1,458 ENCODE DNase-seq biosamples**, covering **369 unique cell types and tissues**.
- **92,176,227 enhancer–gene regulatory interactions** total (abstract rounds to "more than 92 million").
- **Per biosample**: average **40,210 unique predicted regulatory elements** and **63,221 E–G regulatory interactions** (ENCODE-rE2G^Extended, where available: **108,535** interactions, **38,923** enhancers). On average **40,210 predicted "enhancers"** = distal elements predicted to regulate ≥1 gene.
- Coverage: on average **0.89% of the ~3 billion bp** of the human genome is annotated per biosample.
- **Genome build: GRCh38 (hg38)**. Predictions lifted from hg19 where needed for comparator methods.
- **Formats**: predictions distributed via the **ENCODE portal** (BED/bedpe-style interaction files; accessions listed in Supplementary Table 12, not in this PDF). The combined CRISPR datasets are distributed as **`.tsv.gz`**. Arc/locus plots rendered with Gviz/GenomicInteractions and IGV.

**Exact URLs given in the paper** (Data availability + Code availability + Reporting Summary):

*Data*
- ENCODE portal (all predictions + input epigenomics; file accessions in Supplementary Table 12): `https://www.encodeproject.org`
- Aggregated Hi-C Megamap: ENCODE accession **ENCSR492BKE**
- CTCF ChIA-PET for CTCF features: **ENCSR184YZV**, **ENCSR597AKG**; Hi-C for loop features: **ENCFF134HIZ**, **ENCFF173VDJ**
- Combined K562 CRISPR **training** set: `https://github.com/EngreitzLab/CRISPR_comparison/blob/v1.0.0/resources/crispr_data/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz`
- Same training set including all ENCODE-rE2G^Extended features (Synapse): `https://www.synapse.org/Synapse:syn70814727` — also on the ENCODE portal under **ENCSR998YDI**
- **Held-out** CRISPR set: `https://github.com/EngreitzLab/CRISPR_comparison/blob/v1.0.0/resources/crispr_data/EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz`
- Individual CRISPR screens: **GSE120861** (Gasperini 2019), **GSE135497** (Schraivogel 2020), `https://github.com/EngreitzLab/ABC-GWAS-Paper` (Nasser 2021), **GSE171452** (Morris 2023), **GSE129837** (Xie 2019), **ENCFF904ZDX** (Klann 2021)
- K562 DNase-seq used to define candidate elements for Perturb-seq re-analysis: **ENCSR000EOT**
- GTEx eQTL variants + RNA expression for the eQTL benchmark (Synapse): `https://www.synapse.org/Synapse:syn70198274`
- Fine-mapped UK Biobank GWAS data: `https://www.finucanelab.org/data`
- Baseline predictor files for benchmarking (Synapse): `https://www.synapse.org/Synapse:syn58896208`
- Pairwise CRISPRi MYC perturbation data: ENCODE **ENCSR443VTK**; gRNA library **ENCFF118BUP**
- H3K27ac ChIP-seq (Huang et al. 2018): **GSE107726**; CRISPRi–H3K27ac ChIP-seq at MYC: **GSE225157**
- GraphReg / EPIraction prediction accessions: Supplementary Tables 18 and 19 (not in this PDF)

*Code (v numbers from the Reporting Summary)*
- **ENCODE-rE2G (v1.0.0): `https://github.com/EngreitzLab/ENCODE_rE2G`** ← the model repo; also hosts the Snakemake workflow for applying trained models to new datasets
- ABC (v1.1.2): `https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction`
- GraphReg (v1): `https://github.com/karbalayghareh/GraphReg`
- EPIraction (v1.10.1): `https://github.com/guigolab/EPIraction`
- CIA / CCD features: `https://github.com/wangxi001/CIA`  *(note: the Reporting Summary page prints `wangxi001`; the Code-availability paragraph in the main text prints `https://github.com/wangxi001/CIA` — I read it as `wangxi001`, but the glyphs `i`/`l`/`1` are ambiguous at this render resolution; verify before scripting)*
- EpiMap: `https://github.com/KellisLab/EpiMap_GRCh38_linking`
- CRISPR benchmarking pipeline (v1.0.0): `https://github.com/EngreitzLab/CRISPR_comparison`
- eQTL benchmarking pipeline: `https://github.com/EngreitzLab/eQTLEnrichment`
- GWAS benchmarking pipeline: `https://github.com/Deylab999MSKCC/e2g-benchmarking`
- Combined K562 CRISPR data processing (v1.0.1): `https://github.com/argschwind/ENCODE_CRISPR_data`
- Held-out CRISPR data processing: `https://github.com/EngreitzLab/ENCODE_Test_Dataset_Analysis`
- Manuscript analysis code: `https://github.com/EngreitzLab/ENCODE-Distal-Regulation-Paper`
- No Zenodo DOI is given in this paper.

---

### 5. Distance and degree distributions (Extended Data Fig. 6, over all 92,176,227 pairs)

- **Distance from enhancer to target TSS: median 37,091 kb… — careful, the figure text reads "median: 37,091 kb; mean: 85,143 kb", which is clearly in *bp* despite the kb label** (37 kb median, 85 kb mean). Treat as **median ≈ 37 kb, mean ≈ 85 kb**; the unit label in the figure is internally inconsistent — flagging as uncertain.
- Main-text distance bins: **24.1% of predicted interactions are at less than 10 kb; 74.6% at less than 100 kb**, consistent with CRISPR/eQTL/GWAS distance distributions.
- **Genes per enhancer: median 1, mean 1.57** (range up to 33), over n = 58,625,467 enhancers with ≥1 predicted target gene in one biosample. Elsewhere stated as "each enhancer regulated 1.57 genes (median = 1)".
- **Enhancers per gene: median 5.91… ** — Extended Data Fig. 6c states "median: 5.91" over n = 15,585,465 genes (range to 679); the main text states "Across biosamples, each gene had, on average, **5.91 predicted enhancers (median = 3)**". Use **mean 5.91, median 3**; the figure's "median: 5.91" label appears to be a mislabel of the mean.
- **Cell-type specificity**: for a given E–G pair, the **number of biosamples in which it is predicted has median 99, mean 232** (out of 1,458). Predictions are highly cell-type specific.
- Extremes: **SHOX** 77.57 enhancers per biosample (chrX:277,280–739,464); **DNAJC16** 0.005 enhancers per biosample (chr1:15,337,428–15,760,806).
- **Functional skew by regulatory complexity**: genes expressed in >50% of biosamples with **≥5 enhancers per cell type** are enriched for cell-type-specific processes (angiogenesis, neurogenesis, cell adhesion, ion transport, differentiation); genes with **≤1.5 enhancers per cell type** are enriched for housekeeping processes (rRNA processing, mRNA splicing, mRNA processing, DNA repair, DNA methylation). Enrichment ±~2 fold, BH-adjusted P < 0.05.
- Enhancer activity: ENCODE-rE2G enhancers have elevated DNase-seq signal, but **18.5% fall below the global median** of all DNase peaks — weak elements can have real regulatory effects when close to the promoter.

---

### 6. Stated limitations

- **Indirect / trans effects**: all models, including ENCODE-rE2G, perform worse at predicting interactions more likely to be indirect (**P_direct < 0.9**); these tend to be further from the TSS. Grouped AUPRC in Extended Data Fig. 3a: P_direct ≤ 0.9 group = 7,884 pairs, 53 positives, 0.7% positive rate, markedly lower AUPRC than the >0.9 group (2,472 pairs, 418 positives, 16.9%).
- **Small-effect interactions**: pairs with small CRISPR effect sizes are poorly predicted. ≤5% effect: 9,956 pairs, 71 positives (15.1% of positives); >5%: 400 positives (84.3%). CRISPRi is "sometimes underpowered to detect small effect sizes (for example, <25%)".
- **Non-canonical / H3K27ac-negative enhancers**: elements lacking H3K27ac signal are poorly predicted (H3K27ac-negative group: 8,579 pairs, 464 positives, 5.4% positive rate vs 1,777 pairs / 7 positives, 0.4% for the explicitly negative class — note these numbers are at the edge of legibility in Extended Data Fig. 3a, treat the exact split as uncertain).
- **Ubiquitously expressed genes**: interactions with ubiquitous genes are less enhancer-responsive and have a **lower positive rate in the training data (1.5% vs 7.6%)** — a known weak spot. Ubiquitous-expressed group: 4,804 pairs, 68 positives, 1.4%; other: 5,552 pairs, 403 positives, 7.3%.
- **Very long-range interactions missed**: ENCODE-rE2G and all other models **failed to predict three well-characterized distal MYC enhancers almost 2 Mb from the TSS**.
- **Promoter-proximal elements are excluded by construction**: training pairs restricted to **1 kb–1 Mb**, with elements overlapping annotated promoters (±1 kb of TSS) or the target gene body removed; positives lacking H3K27ac signal (CTCF, H3K27me3 or no-H3K27ac categories) were removed to focus on enhancer-mediated regulation. So the model says nothing reliable about promoter-proximal or CTCF-type elements.
- **Repression is not modelled**: positives are defined solely as *significant decrease* in expression on CRISPRi. Elements that "induce rather than maintain expression" may be missed (ref. 16); CRISPRi "could theoretically affect gene expression by alternative mechanisms, such as spreading of the repressive signal". **No repressor/silencer predictions are made — not stated as supported anywhere in the paper.**
- **Training data is one cell type**: the supervised model is trained entirely on **K562**; cell-type coverage of the encyclopedia is limited to the 369 ENCODE cell types/tissues with DNase-seq. The authors note eQTL data "lack cell type resolution", and call for "larger, unbiased CRISPR perturbation and eQTL datasets with resolution for individual cell types".
- **Quantitative effect size is not predicted**: "To accurately predict quantitative effects on target gene transcription, future models must account for the intrinsic strength of promoters, which can vary by several orders of magnitude."
- **One-dimensional chromatin annotations alone are insufficient**: chromHMM misses **40%** and the cCRE catalogue misses **66%** of ENCODE-rE2G enhancers; **30%** and **12%** respectively of CRISPR-validated enhancers from Nasser 2021 (Extended Data Fig. 7).
- Biases in CRISPR screens (limited cell types, biases in element and gene selection, inability to distinguish direct from indirect) "could therefore bias the predictions of the ENCODE-rE2G model in similar ways".

---

### 7. CRISPRi/CRISPRa effect sizes

- CRISPRa is **not used in this paper** — all perturbations are CRISPRi (KRAB-dCas9) or CRISPR-Cas9 deletion.
- **Single-enhancer effect sizes in the training set span −1% to −93%** change in gene expression. The majority of positives exceed 5% effect (**84.3% of positives are >5%**; 15.1% are ≤5%).
- Negatives are defined with ≥80% power to detect **15%** (Perturb-seq/DC-TAP-seq) or **25%** (FlowFISH) decreases, so effects below ~15% are systematically under-ascertained.
- CRISPRi knockdown efficiency assumed **a = 0.13** (from FlowFISH); FlowFISH effect sizes scaled so the strongest 20-guide window at the TSS has an **85% effect**; "promoter CRISPRi typically shows 80–90% knockdown by qPCR". Background signal corrected by a **1.57-fold** factor.
- **Super-additivity**: across **20 CRISPRi tiling experiments**, **10 genes** had the sum of individual enhancer effects **greater than 100%** of expression. At MYC, **19 of 21 enhancer pairs** (7 enhancers e1–e7, 5,776 gRNA pairs, 10,080-vector library, 76 gRNAs) showed significant super-additive interaction (BH-corrected P < 0.05, two-sided t-test). Interaction strength is distance/contact dependent: **16–39% difference from the additive model for pairs 3–89 kb apart with high contact**, vs **2–10% for pairs 107 kb–1.79 Mb apart**. So enhancers within ~100 kb combine super-additively; distant ones combine approximately additively. No sub-additive/redundant effects were observed.
- **Chromatin crosstalk**: CRISPRi at one enhancer quantitatively reduces H3K27ac at nearby enhancers, significantly more for elements 1–10 kb apart than 10–100 kb (P = 1.79×10⁻⁵) or >100 kb (P = 7.37×10⁻⁶); for CRISPR deletions, 1–10 vs 10–100 kb P = 0.0077, 10–100 vs >100 kb P = 0.0026, 1–10 vs >100 kb P = 0.00025.

---

### 8. Disease variants / GWAS / eQTL

- **Benchmark scale**: **>30,000 fine-mapped eQTLs** and **569 fine-mapped GWAS variants linked to a probable causal gene**; **96 traits × ~2,750 variants per trait (15,000–18,000 total)**; **32 GTEx tissues, 273–2,790 distal variants per tissue (PIP > 0.5)**.
- **eQTL benchmark (GTEx v8, SuSiE PIP > 0.5)**: in GM12878 vs **n = 273** fine-mapped distal noncoding eQTLs from LCLs, ENCODE-rE2G^Extended achieved **31-fold** and ENCODE-rE2G **25-fold** enrichment at **15% recall**. Across **11 additional GTEx tissues**, ENCODE-rE2G averaged **6.5-fold enrichment (range 3.5–11.4)** at a **mean recall of 16% (range 5–24%)**.
- **GWAS benchmark**: **29,245 fine-mapped noncoding variants for 94 UK Biobank traits** (application 31063, SuSiE, up to 361,194 white-British individuals, INFO > 0.8, MAF > 0.01%, HWE P > 1×10⁻¹⁰, up to 10 causal variants per ±1.5 Mb region, credible sets pruned at r² > 0.25, MHC excluded). In K562/GM12878, ENCODE-rE2G enhancers showed **10.6-fold enrichment** for erythroid-/immune-related biomarker traits, significantly above ABC's **9.6-fold** (P = 0.0001, t-test, n = 10).
- **Linking variants to causal genes**: across **473 blood biosamples** and **197 noncoding credible sets** across 11 blood traits with an independently-implicated putative causal gene, ENCODE-rE2G reached **precision 68% vs 67% for ABC**, with **recall 51% vs 41%**. A re-analysis across 31 traits (**n = 560 credible sets**) is in Extended Data Fig. 5b.
- **ENCODE-rE2G ∩ PoPS** (intersecting top-2 ENCODE-rE2G genes with top-2 PoPS polygenic-priority-score genes in a 1-Mb window) is the headline variant-to-gene method: it identified a target gene for **2,173 of 13,012 noncoding credible sets (PIP > 0.1)**, covering **8,234 of 48,305 fine-mapped GWAS variants**, a **3.71-fold enrichment** over control variants. On the 569 silver-standard SNP–gene links it reached **precision 79%, i.e. 1.4-fold and 2.3-fold higher than ENCODE-rE2G alone or PoPS alone**, with recall comparable to ABC-Max, including **new predictions for 1,044 credible sets**. It beats OpenTargets L2G and the L2G ∩ PoPS combination.
- **Pleiotropy of enhancers**: GWAS variants often overlap enhancers predicted to regulate **several genes (average 2.9, maximum 34)**, with **59.8% regulating different genes in different cell types** — the authors frame this as the fundamental complication for variant-to-gene mapping.
- Worked example: **rs875741 (PIP = 0.50)**, mean corpuscular haemoglobin, linked to **CPEB4** only in haematopoietic multipotent progenitors and not in B cells, T cells or K562; MCH variants enriched **16.3-fold** in haematopoietic MPP enhancers. Other examples given for **TFRC, KIT, BCL11A** (details in Supplementary Note 5, not in this PDF).
- Heritability: S-LDSC enrichment of **~20–28×** for fine-mapped SNPs in predicted enhancers across 10 matched trait–biosample pairs (Extended Data Fig. 5a; exact per-model values in Supplementary Table 7, not in this PDF).

**On dosage-sensitive genes specifically: not stated in paper.** The paper discusses housekeeping/ubiquitous vs cell-type-specific genes and polygenic priority scores, but I found no analysis of dosage sensitivity (pLI, haploinsufficiency, triplosensitivity) anywhere in the main text, Methods or Extended Data.

---

### Practical notes for a target-shortlisting pipeline

1. **Pull the encyclopedia rather than re-run the model** unless your target cell type is absent from the 1,458 biosamples. Accessions are in Supplementary Table 12, which you'll need to fetch from the Nature supplementary page (not in this PDF).
2. To run on a new cell type: Snakemake workflow at `https://github.com/EngreitzLab/ENCODE_rE2G`, one DNase-seq (preferred) or ATAC-seq experiment, plus the ENCODE Hi-C Megamap (ENCSR492BKE). Eight pre-trained input-combination models ship with it; pick the one matching your available assays.
3. The score is a **probability**, so it is comparable across cell types and safely thresholdable; but treat the 70%-recall threshold as the authors' convention, and remember precision there is only **~53–58%** (DNase-only) — roughly **1 in 2 shortlisted links is a false positive**. Tiering by score band (rather than a single cut) is well supported by the calibrated output, which fits your `docs/plan.md` rule about tiers over false precision.
4. For engineering, the predicted *effect size* is not available — only link probability. Score correlates with effect size at r ≈ −0.45, so it is a weak proxy at best.
5. Expect **~1.6 genes per enhancer and ~3–6 enhancers per gene**, median distance **~37 kb**, with ~75% of links under 100 kb — useful priors for sizing a design window.
6. Super-additivity means a single-enhancer knockdown at a multi-enhancer locus will under-deliver: pairs within ~100 kb combine **super-additively** (16–39% above additive), so multiplexed targeting of co-located enhancers is the higher-yield strategy.
7. Everything here is K562-trained and in-silico; held-out precision across 5 cell types holds at ~53%, which is the number to quote for honesty in your UI.

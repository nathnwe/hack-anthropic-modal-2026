---
title: "UniversalEPI: robust prediction of cell type-specific and differential chromatin interactions from DNA sequence and chromatin accessibility"
authors: "Grover A, Zhang L, Muser T, Häfliger S, Wang M, Yates J, Indilewitsch M-C, Wang Y, Van Allen EM, Theis FJ, Ibarra IL, Krymova E, Boeva V"
journal: "Nucleic Acids Research"
year: 2026
volume_issue_article: "54(10), gkag485"
doi: "10.1093/nar/gkag485"
url: "https://doi.org/10.1093/nar/gkag485"
local_pdf: "/Users/nathanewer/nwe_src/hack-anthropic-modal-2026/docs/Grover et al. - 2026 - UniversalEPI robust prediction of cell type-specific and differential chromatin interactions from D.pdf"
date_read: 2026-10-03
pages_read: "1-19 (all)"
note: "The Supplementary material (Figs S1-S16, Supplementary Sections S1-S13, Supplementary Tables S1-S2, Supplementary Sheet 1) is NOT present in the local PDF; items sourced only from it are flagged below as unavailable."
---

Source: Grover et al., *UniversalEPI: robust prediction of cell type-specific and differential chromatin interactions from DNA sequence and chromatin accessibility*, **Nucleic Acids Research 2026, 54(10), gkag485**, https://doi.org/10.1093/nar/gkag485. Local copy: `/Users/nathanewer/nwe_src/hack-anthropic-modal-2026/docs/Grover et al. - 2026 - UniversalEPI robust prediction of cell type-specific and differential chromatin interactions from D.pdf` (19 pages). I read all 19 pages of the main article. **The Supplementary (Figs S1–S16, Supplementary Sections S1–S13, Supplementary Tables S1–S2, Supplementary Sheet 1) is NOT in the PDF** — several specifics (exact accessions, ENCODE cell-line list, some protocol details) live there and are marked below as "supplementary only, not available".

---

# 1. Model inputs

- **Two modalities only, at inference: DNA sequence + ATAC-seq.** "To infer the Hi-C signal in a test cell type, UniversalEPI requires DNA sequence and cell-type-specific ATAC-seq profiles as input." No CTCF ChIP-seq, no histone marks, no DNase at inference. (ChIP-seq for CTCF/YY1/SP1 is a *training target* for Stage 1, not an input.)
- **Reference genome: GRCh38 / hg38.** ATAC-seq processed against GRCh38/hg38; Micro-C data hg38-aligned.
- **Input is peak-anchored, not a contiguous window.** The model only sees accessible regions (ATAC-seq peaks). Input to Stage 2 = **401 consecutive accessible regions (peaks)**, each **1 kb long, centered on the position of maximum ATAC-seq signal**, spanning in total **~4 Mb**. Per-token tensor: x_i ∈ R^(1000×5) = one-hot DNA (4 channels) + ATAC-seq signal (1 channel), stacked as X ∈ R^(401×1000×5).
- **All information about the DNA between accessible regions is discarded**; intervening context enters only as the distance between peaks. This is the key architectural simplification (and a stated limitation).
- **Preprocessing (ATAC-seq):** ENCODE "processed signal *P*-value bigwig profiles" plus pseudoreplicated peak files. **Peak de-duplication: within any 500 bp genomic window, keep only the peak with the highest read count, discard the rest** → "~175K peaks for each cell line". **Cross-cell-type bigwig normalization to GM12878 as reference:** identify conserved CTCF sites (peaks common to CTCF ChIP-seq tracks of all five cell lines), extract signal at those sites with `deepTools multiBigwigSummary v3.5.3`, compute **TMM (trimmed mean of M-values, EdgeR)** scaling factor per cell line with GM12878 as reference, multiply the original signal track by that factor. Note: this normalization step as described uses CTCF ChIP-seq peaks to define the conserved site set — relevant if you want to normalize a new cell type.
- **Optional auxiliary input: 36 bp mappability track** from the UCSC Genome Browser, linearly embedded to d_model and concatenated to the token embeddings in Stage 2.
- **Training-only inputs/targets:** CTCF/YY1/SP1 narrowPeak ChIP-seq from ENCODE for GM12878, K562, IMR90, HepG2, A549. CTCF peaks split by motif orientation using **JASPAR consensus motif MA0139.1**, PSSM convolved over forward and reverse strands, cutoff from a **1% false negative rate** against intergenic background. Result: "49% of the peaks (18% of all peaks)" had both strands above cutoff; "Overall, 49% of CTCF peaks contained at least one forward-strand motif hit, and 48% contained at least one reverse-strand motif hit. For 20% of the peaks, no corresponding strand could be determined" — those 20% were **dropped from the dataset during training**. ZNF143 was evaluated as a 4th TF and **excluded** (antibody cross-reactivity with CTCF; no accuracy gain).
- **Training augmentation:** 10% of the accessible regions in the training set were replaced by **inaccessible 1 kb regions chosen randomly from the genome**, so the model also behaves correctly on closed regions (e.g. after mutation).

# 2. Model outputs

- **Predicted quantity: Hi-C interaction value on a quantitative (not distance-normalized) scale** — specifically the **log of a smoothed, ICE-normalized, cross-cell-line distance-stratified z-score-rescaled Hi-C contact value**. The paper deliberately **retains the distance effect** ("we aim to capture the true biological interactions and hence, retain the effect of distance in our Hi-C matrices"), unlike Akita/EPCOT (observed-over-expected) and ChromaFold (z-scores).
- **Target construction pipeline (reproduce exactly if you want comparable numbers):**
  1. Raw Hi-C from 4D Nucleome Data Portal for GM12878, K562, IMR90, HepG2; **5 kb resolution** contact maps.
  2. **ICE normalization** with `HiCExplorer v2.2.1`, filter thresholds **−1.1 and 4.5**, applied **per chromosome, independently per cell line** (intra-chromosomal only).
  3. **Distance-stratified robust z-score rescaling to GM12878** as reference (GM12878 chosen for depth: "~759M interactions as compared to ~275M interactions in the other cell lines"). Per chromosome, store median m_d and MAD mad_d of all interactions at distance d, for **d ∈ [0, 800] bins (each bin = 5 kb)** → i.e. up to **4 Mb** of stored statistics. Smoothed with a **B-spline (SciPy `splrep`), smoothing condition s = 5**. Eq. (1): ŷ^K562_{i,j} = ((y^K562_{i,j} − m^K562_{|j−i|}) / mad^K562_{|j−i|}) · mad^GM12878_{|j−i|} + m^GM12878_{|j−i|}.
  4. **Distance-dependent Gaussian smoothing** with four kernels, variances **0.5, 0.65, 0.8, 1.0** for interactions between regions **<25 kb, 25–100 kb, 100–250 kb, >250 kb** apart respectively.
  5. **Logarithmic transformation** of interaction scores.
- **Output geometry:** Y = {y^(i,j)}, i ≤ j, **i,j ∈ {101,…,300}** — i.e. pairwise interactions among the **central 200 of the 401 input peaks**, arranged in a DeepC-style vertical zigzag pole over the center of the sequence window.
- **Resolution: pairwise between ATAC-seq peaks** (peak-pair level), with Hi-C targets at **5 kb bins**. A separate model was trained on **Micro-C at 1 kb resolution** (GM12878, K562, HCT116 from Hong et al., `cooler v0.10.4 zoomify` from 50 bp → 1 kb .cool), same normalization/smoothing pipeline.
- **Maximum genomic distance: 2 Mb** ("predicts EPIs up to 2 Mb apart"; "regulatory elements that are up to 2 Mb apart"). Note benchmarking against other methods was restricted to **within 1 Mb**, "which was the maximum common receptive field across all existing methods".
- **Also outputs uncertainty per prediction:** total predictive variance [σ_T^(i,j)]² = aleatoric + epistemic, Eq. (6); and a derived **"maximum-confidence fold-change" (FC)** between two conditions, which is the headline differential quantity.

# 3. Architecture and inference cost

Architecture (two stages, trained sequentially):
- **Stage 1 (representation learning), frozen after training:** 1D CNN, **five convolutional layers, {30, 60, 60, 90, 90} channels, kernel sizes {11, 11, 11, 5, 5}**, max-pooling widths **{4, 5, 5, 4, 2}** after each conv layer, **LeakyReLU slope 0.2**, **20% dropout**, fully-connected head → 4-dimensional output = maximal ChIP-seq intensity for **CTCF-forward, CTCF-reverse, SP1, YY1**. Loss: MSE (Eq. 2). Input per peak: 1 kb DNA + matching ATAC-seq. Linear projection heads map each of **N_L = 5** conv-block feature layers to dimension **C = 180**, giving h(X) ∈ R^(401×C×N_L).
- **Stage 2 (Hi-C prediction):** **stochastic gating feature selection** (Yamada et al.) applied at the *feature* level with σ = 0.5, aggregated by average pooling; sparsity regularization **λ = 0.01**; then a **transformer encoder: 4 multi-head attention blocks, 4 heads, dropout 0.1, d_model = 32**, with a **genomic-distance-aware sine-cosine relative positional encoding** (max distance **3 Mbp**, resolution **500 bp**, signed — negative to the left of the center token, positive to the right, +3 Mbp constant offset); token position = position of the peak's maximum ATAC signal; middle token is the reference. Decoder: **2 MLP layers with ReLU** on the concatenation of encoder outputs e_i, e_j. Sequences are **randomly flipped during training** for orientation invariance.
- **Loss:** β-NLL with **β = 0.5** (Eq. 5) for aleatoric uncertainty; **deep ensemble of K = 10** models with different random inits for epistemic uncertainty.
- **Size: 2.5M parameters** (single model; Fig. 3d). Compare: Akita 0.7M, EPCOT 13.1M, C.Origami 13.3M. ⚠️ Fig. 3d labels *both* "UniversalEPI (single model)" and "UniversalEPI (ensemble)" as 2.5M, which looks like a figure labelling artefact (a 10-model ensemble should be ~25M). Flagging as uncertain.
- **Training cost:** Stage 1 Adam, lr 1e−4, weight decay 1e−4, batch 1024, 100 epochs. Stage 2 AdamW, lr 1e−3, batch 32, **20 epochs**; positional-embedding upscaling by √d_model. **"Training the Hi-C and Micro-C prediction model took 12 h on 0.2 of an A100 GPU with a batch size of 32, and 24 h on a single RTX2080Ti with a batch size of 16."**
- **INFERENCE COST (the number you want), Fig. 3d — runtime to predict all test chromosomes of HepG2:**
  - **UniversalEPI single model: 8.4 minutes**
  - **UniversalEPI ensemble (K=10): 84.5 minutes**
  - Akita: 30.6 min; EPCOT: 35.0 min; **C.Origami: 3662.3 min (~61 h)**
  - Graphical abstract restates this as "6 h" (C.Origami) vs "0.1 h" (UniversalEPI single) vs "1.4 h" (ensemble) — ⚠️ the "6 h" in the graphical abstract is inconsistent with 3662.3 min = 61 h in Fig. 3d; trust Fig. 3d's explicit numbers. The 0.1 h / 1.4 h figures match 8.4 min / 84.5 min.
  - **Per-locus runtime is not stated in the paper.** But: whole test-chromosome-set inference at 8.4 min single-model implies a single 4 Mb window (401 peaks → 200×200 upper-triangle pair predictions) is sub-second on GPU. That is an inference by me, not a paper statement.
  - **Hardware for the inference benchmark is not stated** (training used A100 / RTX2080Ti). The paper says the model "can run on standard GPUs". **CPU runtime: not stated in paper.**
  - Practical demo note: for a single locus you only need the ATAC bigwig + peaks around that locus, a 401-peak window, and ideally the 10-model ensemble (the ensemble is what gives the uncertainty and the maximum-confidence FC). Ensemble = 10× the single-model cost but the 10 forward passes are trivially batchable.

# 4. Cell-type specificity and differential prediction

- **Cell-type specificity comes entirely from the ATAC-seq channel** (and the peak set derived from it). DNA sequence is constant; the Stage-1 CNN is trained to predict only three **ubiquitously expressed** TFs (CTCF, YY1, SP1) *by design*, "to prevent the model from capturing DNA motifs of other, potentially cell-type-specific, TFs and thereby enable the model to generalize predictions across different cell types in the second stage." So the model is deliberately made *not* to memorize cell-type-specific motif grammar; specificity is injected via accessibility.
- Demonstrated: "Using cell-line-specific ATAC-seq profiles as inputs, UniversalEPI obtained equally high values of distance-stratified correlation between predicted and ground-truth Hi-C values for all four test cell lines."
- **Differential chromatin interaction, operationally:** run the **same pretrained ensemble twice**, once with condition A's ATAC-seq (and, if modelling a variant, condition A's DNA sequence) and once with condition B's, over the **same peak coordinate set**, then compute the **log₂ fold-change of predicted Hi-C per peak pair**, and refine it to the **"maximum-confidence log₂ FC"** using the ensemble's predicted uncertainty — i.e. the most conservative FC consistent with the uncertainty intervals of the two predictions, which is **0 when the two predictions' uncertainty intervals overlap**. "Using the reported maximum-confidence fold-change (FC), a unique feature of UniversalEPI that removes inherent data noise, the user can focus on EPIs that are significantly altered across conditions." The two conditions can be: two cell types, two time points, two alleles (SNP ref vs alt), or pre/post an in-silico sequence edit.
- Worked examples: macrophage activation (THP-1 + LPS + IFNγ, ATAC-seq + Hi-C at 8 time points, Reed et al.) — **zero-shot** with Model 1 (trained on GM12878 + K562), Pearson's R and Spearman's ρ with ground truth Hi-C **exceeding 90%** at all time points; gain/loss/static loop classification from 0 h vs 24 h ATAC alone (n = 84 loss, 7435 static, 68 gain). Also SNP rs1800734 (G→A) in isogenic COLO320 lines → predicted stronger interaction between the SNP-containing region and the *DCLK3* 3′ UTR in the alternate allele, matching experiment. Also pseudo-bulk scATAC-seq from 8 esophageal adenocarcinoma tumours (differentiated vs undifferentiated states) → differential *LCOR* / *MECOM* promoter interactions.

# 5. Benchmarking

**Compared against:** C.Origami (original, with CTCF ChIP-seq input; and retrained without CTCF for fair input parity), EPCOT (authors' pretrained model), Akita (sequence-only), ChromaFold (scATAC-seq setting), ChINN (anchor-based differential-loop classifier), plus three trivial baselines: **Distance**, **Median**, **Swap** (Supplementary Section S5).

**Metrics:** Spearman's correlation (primary) and Pearson's correlation on smoothed log Hi-C; **distance-stratified Spearman's correlation** computed per distance d in [0, 2 Mb] in **5 kb steps** (removes the trivially strong distance effect); MSE vs uncertainty for calibration; balanced accuracy / macro F1 / macro AUC-ROC for loop-type classification.

**Evaluation protocol:** unseen chromosomes **and** unseen cell types. Chromosome split: **validation = chr5, 12, 13, 21 (~15% of genome); test = chr2, 6, 19 (~15%)**; rest for training. **Model 1** trained on GM12878 + K562, tested on IMR90 + HepG2. **Model 2** trained on IMR90 + HepG2, tested on GM12878 + K562. Bins overlapping the Kelley et al. unmappable (`umap_k36_t10_l32_hg38.bed`) and blacklist (`hg38.blacklist.rep.bed`) regions were flagged and interactions containing a flagged endpoint removed. Benchmarking restricted to **< 1 Mb** (max common receptive field).

**Where it performs well:**
- Only method to reach **Spearman's ρ > 0.90 on all test cell lines** (Fig. 3b); discussion restates ">0.9 in all experimental settings".
- Significantly outperformed C.Origami-with-CTCF (**P < .05**), C.Origami-without-CTCF (**P < .01**), Akita (**P < .001**); Wilcoxon signed-rank.
- Best **distance-stratified** correlation at essentially all distances (Fig. 3c) — i.e. it is not just reproducing the distance decay.
- Per-interaction-class Spearman (Fig. 3e), UniversalEPI: IMR90 PPI 0.93 / dEPI 0.92 / pEPI 0.91 / sEPI 0.89; HepG2 0.95 / 0.94 / 0.93 / 0.93; GM12878 0.95 / 0.93 / 0.93 / 0.93; K562 0.94 / 0.93 / 0.92 / 0.91. (dEPI = distal EPI > 2 kb, pEPI = proximal EPI 200 bp–2 kb, PPI = promoter–promoter, sEPI = super-enhancer–promoter from dbSuper.)
- Micro-C at 1 kb: higher Spearman than EPCOT on held-out HCT116, "This increase was predominantly observed for shorter interactions (<100 kb) while achieving comparable performance for longer-distance interactions (150–200 kb)". Also outperformed Akita and EPCOT after converting to observed-over-expected (divide by distance-stratified mean prediction).
- Uncertainty is calibrated: "higher total uncertainty was generally associated with a greater mean-squared error".
- Ensemble > single model significantly on unseen cell lines (Fig. 2c, Mann–Whitney U, *P* ≤ .05).

**Where it is weaker / failure modes:**
- **On IMR90 and K562, UniversalEPI achieved comparable results to C.Origami** (it only "substantially outperformed C.Origami on HepG2 and GM12878"). In Fig. 3e, C.Origami-with-CTCF is slightly *better* on K562 (0.95/0.92/0.93/0.93 vs 0.94/0.93/0.92/0.91 — mixed). "UniversalEPI outperformed both versions of C.Origami on all but the K562 cell line, where C.Origami with additional CTCF binding information showed slightly higher prediction accuracy."
- Stage-1 TF prediction accuracy on unseen chromosomes/cell types is only moderate: **average Pearson's correlation of 0.78 (CTCF), 0.62 (YY1), 0.47 (SP1)**.
- Distance-stratified correlation falls to roughly 0.4–0.6 in absolute terms (Fig. 3c) once the distance effect is removed — the ρ > 0.9 headline number is inflated by the distance signal.
- Loop-type classification is modest in absolute terms: balanced accuracy **0.57** (vs ChINN 0.34, random 0.33), macro F1 **0.22** (ChINN 0.27 — **ChINN is better on macro F1**), macro AUC-ROC **0.68** (ChINN 0.52, random 0.50). Confusion matrix: of true gains, 37 predicted gain vs 17 static vs 14 loss; of true losses, 9/11/64. Gain/loss classes are tiny (68/84 vs 7435 static).
- **Hi-C at 5 kb contains bins with very few reads**; the authors "removed these and their one-hop neighboring bins from all datasets". So low-coverage bins are out of scope.
- Micro-C "high-resolution data is inherently noisier than the Hi-C data" — the Micro-C model was trained on the smoothed contact matrix.

**Stated generalization limits:**
- Explicitly: "Additional validation experiments may be required to check whether the model can show equally high performance in all cell types and whether highly specialized, cell-type-specific mechanisms that regulate chromatin interactions exist." Validation covered only lung fibroblasts, B lymphocytes, macrophages, CML, and hepatocellular carcinoma cells.
- The model is blind to all non-accessible chromatin and to cell-type-specific TF motifs by construction; intervening insulators/regulatory elements between two anchors are only implicitly captured via the peaks in the 401-peak window.
- Distance: ≤ 2 Mb by design; benchmarks only to 1 Mb.
- Future work suggested: adversarial training to further prevent cell-type-specific information leakage into Stage 1.

# 6. TAD / boundary / CTCF handling

- **No explicit TAD or boundary prediction head.** TADs are only used as an *analysis* construct: 293 TADs across the four cell lines were taken and one endpoint's convergent CTCF site inverted in silico.
- **CTCF motif orientation IS used, and is central.** CTCF ChIP-seq peaks are split into **forward** and **backward** sets using JASPAR **MA0139.1**, and Stage 1 predicts forward and reverse CTCF occupancy as **two separate targets** — explicitly "Since the orientation of CTCF binding to chromatin plays an important role in determining the stability of TADs". DeepLIFT attributions confirmed the Stage-1 CNN learned motifs "closely resembling their respective known consensus motifs".
- **Insulation score: not computed / not stated in paper.** No insulation-score or TAD-calling metric is reported anywhere in the main text.
- The TAD-boundary result: inverting one convergent CTCF site at a TAD boundary produced "predominantly negative values of the log FC across all four cell lines", i.e. reduced interaction at the original boundary. Worked example in IMR90 at **chr6:17,280,000–17,585,000**: substantial reduction between the original TAD boundaries and between the 3′ boundary and several upstream enhancers, plus a **significant increase** in interactions with a **reverse CTCF motif located nearly 200 kb downstream** of the original 3′ boundary.

# 7. Code / data / licence — exact URLs

- **Source code: https://github.com/BoevaLab/UniversalEPI**
- **Zenodo deposit: https://doi.org/10.5281/zenodo.14622040**
- **Tutorial / wiki: https://github.com/BoevaLab/UniversalEPI/wiki**
- **Precomputed Hi-C predictions for 157 ENCODE ATAC-seq datasets (116 cell lines + 41 primary cells), as a UCSC Genome Browser track hub: https://genome.ucsc.edu/cgi-bin/hgHubConnect** — two track types: (i) log-transformed ICE-normalized predictions (quantitative), (ii) z-score-normalized tracks highlighting long-range interactions. The mapping between generated interaction files and ENCODE datasets is in **Supplementary Table S2 (not available in this PDF)**. ⚠️ The hub URL given is just the generic hgHubConnect page, not a direct hub .txt URL — you'll have to find the hub entry via the wiki/Zenodo.
- Primary data sources: ENCODE https://www.encodeproject.org/, 4DNucleome https://data.4dnucleome.org/, GEO https://www.ncbi.nlm.nih.gov/geo/ (accessions in Supplementary Sheet 1 / Supplementary Table S1 — **not available**). EAC scATAC/scRNA from dbGaP https://www.ncbi.nlm.nih.gov/gap/, accession **phs003438.v1**.
- **Article licence: CC BY-NC 4.0** (https://creativecommons.org/licenses/by-nc/4.0/) — non-commercial. **The software licence of the GitHub repo is NOT stated in the paper** — check the repo directly before building anything commercial on it.
- **Model weights: the paper does not explicitly say "weights are released"**; it says source code is available and deposited to Zenodo. Weights are presumably in the Zenodo deposit but this is **not stated in paper** — verify.

# 8. Perturbation prediction: claims and caveats

**Claims made:**
- "This model now makes it possible to carry out *in silico* experiments to assess changes in chromatin folding under diverse biological scenarios."
- **CTCF motif inversion:** 293 TAD-boundary convergent CTCF sites inverted in silico (sequence edit → re-run → log₂ FC). Predominantly negative log₂ FC in all four cell lines; maximum-confidence FCs agreed with the raw FCs.
- **SNP / allele-specific:** rs1800734 (G→A) using isogenic COLO320 ATAC-seq from Liu et al. — correctly predicted stronger interaction between the SNP region and the *DCLK3* 3′ UTR on the alternate allele. "UniversalEPI can detect allele-specific chromatin interaction differences **when paired with corresponding chromatin accessibility data**" — note the inputs changed in *both* channels (sequence and ATAC).
- Future application proposed: "study changes in chromatin organization in the presence of non-coding mutations and structural variations"; "apply the same principles to large cancer datasets such as TCGA".
- Training augmentation with 10% inaccessible regions was explicitly done "ensuring that the model can also predict accurately for various inaccessible regions, which can arise due to mutations".

**Caveats the authors state:**
- **Predicted strong EPIs are not strictly functional:** "As UniversalEPI is trained to predict the Hi-C signal, the strong EPIs predicted by UniversalEPI are not strictly functional. It is known that enhancers and promoters can gain chromatin accessibility and interact even before gene expression is activated. Additional data, such as profiling of histone H3 lysine 27 acetylation (H3K27ac) or enhancer RNA (eRNA) transcription, may be needed to assess the functionality of predicted strong EPIs." **This is the single most important caveat for a target-shortlisting pipeline.**
- The model "operates on the DNA sequence within open chromatin. By ignoring information in closed chromatin, other than assessing the distance between open chromatin regions..." — a perturbation that creates a *new* accessible site cannot be modelled unless you also supply the changed ATAC signal; i.e. **the pipeline cannot predict the accessibility consequences of an edit, only the looping consequences given an assumed accessibility state.**
- Validated only on 5 cell types; generalization to highly specialized cell types unverified.
- The authors do **not** report any validation of predicted *synthetic loop insertion* or large deletion/insertion — only point inversions of CTCF motifs and a SNP. **Deletions, insertions, and synthetic loop design are not evaluated in the paper.**
- Differential calls should be filtered by the maximum-confidence FC; the authors always do this ("retention of only the highly reliable differential interactions").

---

# Practical implications for a target-shortlisting demo

- **Feasible live:** single-model inference for one 4 Mb locus should be well under a second on GPU (extrapolating from 8.4 min for 3 whole test chromosomes); the 10-model ensemble is ~10× and is required for uncertainty/maximum-confidence FC. Precompute the ENCODE track hub predictions where possible instead of running the model.
- **Pipeline prerequisites per cell type:** ATAC-seq P-value bigwig + peak calls on hg38, 500 bp peak de-duplication, TMM normalization to GM12878 using conserved CTCF sites (this last step needs a reference CTCF site set — check the repo for a shipped file, since the paper's recipe requires CTCF ChIP-seq from the five training cell lines).
- **For differential/EPI shortlisting:** output is peak-pair log Hi-C with per-pair uncertainty; the shortlist criterion the authors use is the maximum-confidence log₂ FC (non-zero ⇒ uncertainty intervals don't overlap).
- **Honesty flag for your CLAUDE.md rules:** these are in-silico, unvalidated predictions of *contact*, explicitly not of *function* (authors' own words). Tiered uncertainty is directly available from the model — use its aleatoric+epistemic total uncertainty rather than inventing a confidence score.

**Unsure / not found:** the ensemble parameter count in Fig. 3d (labelled 2.5M, likely a figure error); the inference hardware for the Fig. 3d benchmark; CPU inference feasibility; per-locus runtime; the software licence; whether trained weights are in the Zenodo deposit; the direct URL of the UCSC track hub; all Supplementary content (accessions, Sections S1–S13, Table S1/S2, Sheet 1).

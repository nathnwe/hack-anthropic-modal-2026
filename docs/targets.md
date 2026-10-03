# Deliverable 1 — Key target examples: diseases of gene dosage

**Status:** first pass, built from primary curated sources on 2026-10-03. Every row traces to a
source. Columns marked *(annotation)* are our judgement, not curated fact, and are the ones a
reviewer should attack first.

## Source of truth

ClinGen Dosage Sensitivity Map, downloaded 2026-10-03 (dated `03 Oct, 2026`), GRCh38/hg38:

- Genes: <https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv> — 1,532 genes curated
- Regions: <https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv> — 513 regions curated

ClinGen scores both directions on a 0–3 scale. We only treat **score 3 = "sufficient evidence for
dosage pathogenicity"** as established. Coordinates below are ClinGen's, unmodified.

## The headline asymmetry (and why it shapes the tool)

| | Curated at **gene** level, score 3 | Curated at **region** level, score 3 |
|---|---|---|
| **Haploinsufficiency** ("too little") | **420 genes** | 170 regions |
| **Triplosensitivity** ("too much") | **3 genes** — `LMNB1`, `PLP1`, `TPSAB1` | **20 regions** |

Only one gene, `PLP1`, is score 3 in both directions.

This is not because gain-of-dosage disease is rare — it is because the evidence usually arrives as a
**duplication CNV spanning many genes**, so ClinGen curates the interval rather than assigning blame
to one gene. Three consequences for us:

1. **The "too much" list must be assembled from region curations plus the specific duplication
   syndromes where the driver gene is independently established.** We cannot just filter a column.
2. **Driver attribution inside a duplicated interval is a real, unsolved inference step** — and it is
   exactly the kind of place our tool should report uncertainty rather than a confident single gene.
3. The two directions need different interventions: too-little wants a *new* enhancer–promoter
   contact (staple on), too-much wants an *existing* contact broken or insulated (staple elsewhere,
   or re-anchor a boundary).

## Why dosage disease is the right fit for loop engineering

Three properties make a dosage gene tractable for an E–P staple, and they are the filters the
pipeline applies:

- **Monotonic dose–response.** Phenotype tracks expression level, so a partial, tunable change is
  therapeutic. Contrast a dominant-negative missense allele, where halving total output does nothing
  useful.
- **Distal, cell-type-specific CREs.** The gene is actually enhancer-driven, so there is something to
  staple. A gene run entirely from its own promoter offers no handle.
- **A therapeutic window wide enough to aim at, narrow enough to matter.** If any overexpression is
  tolerated, a blunt tool already works and we add nothing.

`MECP2` satisfies all three and is pathogenic in **both** directions — see below.

## "Too little" — loss of one functional copy (goal: raise output of the remaining allele)

All rows below are ClinGen **haploinsufficiency score 3**. Coordinates GRCh38.

| Gene | Location (GRCh38) | Disease | ClinGen HI PMIDs | Enhancer-tractable? *(annotation)* |
|---|---|---|---|---|
| `MECP2` | chrX:154,021,573-154,097,717 | Rett syndrome | 31206249, 10508514, 16647997 | **Yes — hero case.** Narrow window; both directions pathogenic |
| `SCN1A` | chr2:165,984,641-166,149,161 | Dravet syndrome | 11359211, 11940708, 18930999 | **Yes.** Upregulating the wild-type allele is already an active clinical strategy |
| `BCL11A` | chr2:60,450,520-60,553,924 | — (here the *target* of de-repression, not the disease gene) | 27453576, 25938782, 25979662 | **Yes — precedent anchor.** See below |
| `LDLR` | chr19:11,089,463-11,133,820 | Familial hypercholesterolaemia | 33740630, 15321837 | Yes. Hepatocyte; dose–response to LDL-C is well established |
| `NF1` | chr17:31,094,927-31,377,677 | Neurofibromatosis type 1 | 1757093, 1302608, 34427956, 34750850, 28230061 | Plausible; large gene, many CREs |
| `PTEN` | chr10:87,863,625-87,971,930 | PTEN hamartoma tumour syndrome | 21194675, 32003824 | Plausible; **raising a tumour suppressor is directionally safe** |
| `BMPR2` | chr2:202,376,327-202,567,749 | Pulmonary arterial hypertension | 10903931, 11015450 | Plausible; penetrance is already dose-modified |
| `TBX5` | chr12:114,353,911-114,408,442 | Holt-Oram syndrome | 24664498, 10077612, 12818525 | Developmental window likely closed post-natally — **poor fit, keep as negative control** |
| `FOXG1` | chr14:28,766,787-28,770,277 | FOXG1 syndrome | 28661489, 24836831, 29655203, 24139857, 19578037, 18571142 | Known large regulatory landscape |
| `RUNX1` | chr21:34,787,801-35,049,302 | Familial platelet disorder with AML predisposition | 31723833, 18478040, 19357396, 33616470 | Yes; haematopoietic, accessible cell type |
| `NSD1` | chr5:177,131,798-177,300,213 | Sotos syndrome | 11896389, 17565729, 23190751 | Bidirectional interest — the 5q35 region is also TS3 |
| `HNF1B` | chr17:37,686,431-37,745,059 | Renal cysts and diabetes | 10720943, 15068978, 9398836 | Bidirectional — 17q12 region is also TS3 |

The remaining ~408 HI-3 genes are the long tail the pipeline is meant to triage; these twelve are the
hand-picked demo set.

### The precedent that makes the whole thing credible: `BCL11A`

The one approved therapy that works by this logic is **exatamcel/exagamglogene autotemcel
(Casgevy)**, which disrupts the **erythroid-specific enhancer in intron 2 of `BCL11A`** to lower
BCL11A in red cell precursors and de-repress fetal haemoglobin, treating sickle cell disease and
β-thalassaemia. It is worth stating plainly in the pitch: *regulating a gene by acting on a distal
CRE rather than the coding sequence is already an approved medicine.* Our contribution is not the
concept — it is choosing **which** contact to make or break, for genes where nobody has done the
decade of work that `BCL11A` took.

Note the useful subtlety: that therapy *lowers* BCL11A in order to *raise* HBG. Direction of effect
on the disease gene and on the targeted gene can differ, and the schema should not assume they match.

## "Too much" — extra functional copy or over-activation (goal: lower output, or insulate)

Gene-level TS-3 is only three genes, so the table below mixes gene-level and region-level evidence.
**The evidence level column is the honest part — do not drop it.**

| Driver gene | Location (GRCh38) | Disease | ClinGen evidence level | PMIDs |
|---|---|---|---|---|
| `LMNB1` | chr5:126,776,623-126,837,020 | Autosomal dominant leukodystrophy | **Gene TS score 3** | 16951681, 23649844, 28769756, 30697589 |
| `PLP1` | chrX:103,776,506-103,792,619 | Pelizaeus-Merzbacher disease | **Gene TS score 3** (also HI 3) | 15689360, 10417279, 19376225, 9634530 |
| `TPSAB1` | chr16:1,240,705-1,242,554 | Hereditary alpha-tryptasemia | **Gene TS score 3** | 27749843, 32717252, 32777817, 33465452 |
| `PMP22` | chr17:15,229,779-15,265,326 | Charcot-Marie-Tooth 1A | Region TS 3 (17p12, chr17:14,194,598-15,519,638) | 20301384, 20301532 |
| `MECP2` | chrX:154,021,573-154,097,717 | MECP2 duplication syndrome | Region TS 3 (Xq28, chrX:154,008,529-154,110,279) | 32043567, 29618507, 22679399, 29141583 |
| `APP` | chr21:25,880,550-26,171,128 | Early-onset Alzheimer's / Down syndrome dosage | Region TS 3 (21q21.3); **gene-level TS only 1** | 19684239, 21193246, 37170141 |
| `RAI1` | chr17:17,681,458-17,811,453 | Potocki-Lupski syndrome | Region TS 3 (17p11.2) | 28837307 |
| *ZRS enhancer* | chr7:156,791,102-156,791,874 | Preaxial polydactyly | Region TS 3 — **a CRE, not a gene** | 18178630, 18417549, 19291772 |

### The ZRS row is the most interesting thing in this table

`ZRS` is a ~800 bp limb enhancer roughly **1 Mb** from its target `SHH`. Duplicating the enhancer —
not the gene — causes disease. It is simultaneously:

- the textbook demonstration that **long-range E–P contact is the unit of gene regulation**, which is
  the premise of the whole project;
- a **dosage disease whose causal element is an enhancer**, so the only sane intervention is at the
  contact level; and
- a natural **validation case**: the ground truth of which enhancer talks to which promoter across
  1 Mb is known from human genetics, so if our pipeline cannot recover `ZRS → SHH`, it is wrong.

Proposed: use ZRS→SHH as a **positive control the pipeline must pass** before we believe any novel
prediction. That is a concrete, cheap way to meet the "knows when it's wrong" framing of Track 2.

## The hero gene: `MECP2`

| | |
|---|---|
| Too little | Rett syndrome — ClinGen HI score 3 |
| Too much | MECP2 duplication syndrome — Xq28 region TS score 3 |
| Location | chrX:154,021,573-154,097,717 |

Why it carries the demo:

- **It is the argument for tunability, in one gene.** Gene replacement for Rett risks overshooting
  into duplication syndrome. A therapy that can only go to one fixed level is dangerous here; a
  staple that can be dialled and switched off is the right shape of tool.
- **It exercises both pipeline directions** from one locus, so the demo shows "staple on" and
  "insulate/break" without a second data pull.
- **It shows the curation asymmetry honestly** — one direction gene-level, one region-level. Good
  material for the uncertainty panel rather than something to hide.
- Caution to state out loud: X-linked, so X-inactivation mosaicism in females complicates what
  "dose" even means per cell. This belongs in the limitations panel, not buried.

## Open questions for the team

1. **Demo set size.** Plan says 10-15 genes. The twelve HI-3 rows plus eight too-much rows is twenty.
   Trim to the ones where ENCODE-rE2G actually has a matching biosample — that filter should drive
   the final list, not our taste.
2. **Cell type is the binding constraint, not the gene.** Every row needs a cell type where (a) the
   gene is expressed, (b) rE2G has predictions, (c) a staple could plausibly be delivered. We should
   pick cell types first and let that prune the gene list.
3. **Is `TBX5` worth keeping as a deliberate negative?** A pipeline that rejects targets is more
   credible than one that scores everything. Recommend: yes, keep one or two.
4. **Driver attribution inside TS-3 regions** needs an explicit "we don't know which gene" output
   state for multi-gene intervals.

## Provenance note

ClinGen scores, coordinates, PMIDs and region definitions in this file are copied from the two TSVs
above and are machine-checkable against them. The *(annotation)* column, the tractability judgements,
the ZRS control proposal and the MECP2 argument are **ours**, written 2026-10-03, and are not
externally sourced. The Casgevy/BCL11A mechanism statement should be given a citation before it goes
in the pitch.

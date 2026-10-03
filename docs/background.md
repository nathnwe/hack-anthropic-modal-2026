# Background

## The question
Can we reversibly modulate disease-relevant gene expression by reshaping 3D genome organisation, ideally using existing regulatory machinery?

## Why
- Many promising targets are inaccessible by direct protein targeting.
- Gene and epigenetic therapies usually alter DNA sequence or add new regulatory instructions. Off-target risk and unknown long-term consequences make regulatory approval hard (cf. the 2018 He Jiankui CCR5 edited babies).
- 3D genome organisation is a poorly exploited, low-footprint, reversible level of regulation. Reshaping which enhancers contact which promoters (enhancer-promoter, "E-P", contacts) changes expression without changing the sequence.
- It matters in disease: Dong et al. (2026) linked 3D genome changes during renal cancer progression to activation of a specific oncogene.
- The ingredients exist: cells already use these contacts; programmable CRISPR-based tools can bring chosen DNA regions together ("stapling"); genome maps and AI models are making regulatory architecture predictable.
- **What is missing** is a design system that turns those data into actionable predictions: for a given gene and cell type, which contact should we alter, and what change in expression should we expect?

## What we are building
A **target-shortlisting tool for 3D genome engineering**.

### Deliverables
1. **Key target examples.** Diseases caused by gene dosage, split into too much of a gene vs too little.
2. **Explainer animation (Blender).** A chromatin loop forms in the cell, an enhancer contacts a promoter, expression of the gene rises (increased dosage).
3. **Pipeline.** Based on two papers in this folder: UniversalEPI (Grover 2026) and ENCODE-rE2G (Gschwind 2026).
   1. A human researcher names a dosage-sensitive gene driving a disease.
   2. The model works out: chromosome; TAD; regulating CREs; which two regions to staple to tune expression up, down or off; secondary risks (e.g. does the staple bring a strong enhancer next to an oncogene promoter?).
   3. Considerations are aggregated into a human-readable score/rank (needn't be purely numerical).
   4. Human and model choose the stapling method given cell type, location, reversible vs permanent, etc.
   5. **Dashboard** (Astro, GitHub Pages): wireframe first; shows the workflow, the ranked gene list, and a visualisation of TAD E-P contacts.

## Hackathon: Anthropic x Modal 2026
Tracks: (1) Open Maths Problems (C3); (2) **Originator**: agents that do science and know when they're wrong / reward hacking (benchmarks, lab automation and safety, epistemological agents with calibrated uncertainty); (3) Drug and protein design (Serova: peptide-HLA stability with foundation models); (4) Materials manufacturing (Polaron: SEM batch QC). Sponsor challenges: **Best use of Modal**, Best use of Devin (reproduce a paper, then go further).

**Our track:** Track 3, Drug and protein design: we are developing a new therapeutic approach (reversible 3D-genome "stapling" of a chosen enhancer-promoter contact), with the tool shortlisting and ranking the targets. We also aim for the Modal challenge (UniversalEPI inference on Modal GPUs). Track 2 is no longer the target, but we still surface uncertainty and risks rather than a black-box score, because it is good practice for a therapeutic claim.

### Judging (100 pts, 20 each)
Technicality, Creativity, Usefulness, Demo, Track/sponsor alignment.

### Format
- Submit: 2-minute demo video, GitHub repo link, short description.
- Round 1 (5 min): pitch 1:30 (problem, approach, high-level scheme), live demo 1:30, Q&A 2:00.
- Round 2 (8-10 finalists, on stage, 3 min): pitch 1:30, demo 1:30, no questions.

### Rules to keep in mind
- Build entirely during the event; **no prior commits to the repo**. Our repo already has two setup commits ("uv and git", "Key papers added") made before scaffolding. Ask the organisers whether setup commits and background papers are acceptable, or start a fresh repo if not.
- Teams of 2-5. Open-source libs, APIs and pre-trained models are fine if credited (keep the README credits current).
- Never use data without consent; follow the law, ethical AI practice and Code of Conduct.
- Use partner and sponsor APIs per their terms of service.
- Submit all code and demos before the deadline.

## What the judges should see (demo shape)
One hero gene end to end: researcher picks a gene -> TAD and CRE map -> proposed staple with direction of effect -> risk flags -> ranked tier with plain-English rationale -> recommended stapling method. Plus an explicit "known limitations" panel: predictions are in silico and not experimentally validated.

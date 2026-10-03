# Storyboard v2 (draft for review)

30 s, 24 fps, 720 frames, 4K final. Soft pastel, glossy membrane. Generic hero gene (no real gene named). Minimal on-screen text.

## Decisions so far
- Whole TAD visible throughout, proteins enlarged and not to scale. The TAD is a small illustrative one (about 100-1000 bp). State "not to scale" in the notes.
- CTCF and cohesin already sit at the TAD base (loop pre-formed) as we zoom in; neighbouring loops carry them too.
- All chromatin is one continuous fibre: neighbouring domains are joined to each other and to the hero TAD, ends run off-screen.
- Staple = two dCas9 halves fused directly (no linker). One half binds beside the enhancer first; stapler and enhancer then move together to the promoter, where the other half binds and stays locked.
- Repressor element (silencer) drawn as a generic repressor assembly, marked illustrative.

## Timeline

| # | Time (s) | Frames | Scene | Camera | Action | On-screen text |
|---|---|---|---|---|---|---|
| 1 | 0-3 | 1-72 | Whole cell | Wide, push-in toward the nucleus | Pastel cell, glossy membrane; nucleus (envelope, pores, nucleolus, chromatin), rough ER with ribosomes, Golgi, centrosome and microtubules, mitochondria with cristae, lysosomes, peroxisomes, vesicles drifting | none |
| 2 | 3-6 | 73-144 | Inside the nucleus | Dive through the envelope, settle on one region | Chromatin fibres jiggle; one domain lights up, the rest dims | none |
| 3 | 6-12 | 145-288 | One TAD, baseline | Frames the whole TAD, slow drift | Loop already held by CTCF and cohesin at the base. General transcription factors land on the promoter and recruit Pol II. Several Pol II fire down the gene body; mRNA emerges from each | "Gene" |
| 4 | 12-16 | 289-384 | Enhancer boost | Same framing | Enhancer contacts the promoter complex; Pol II recruitment and firing rate rise; mRNA output rises | "Enhancer elements can increase transcriptional output" |
| 5 | 16-20 | 385-480 | Repressor element | Same framing | Enhancer releases; a repressor element (silencer) moves next to the promoter; fewer Pol II recruited, output falls | "Repressor elements can decrease transcriptional output" |
| 6 | 20-29 | 481-696 | Engineered staple | Pull back slightly to see both sites, then push in | 20-22: the fused dCas9-dCas9 stapler drifts in. 22-23: one half binds next to the enhancer. 23-25.3: stapler and enhancer move together to the promoter. 25.3-26: the other half binds next to the promoter and locks. 26-29: Pol II traffic and mRNA output rise and stay high | "Our stapler molecule enforces co-localisation of promoters and regulatory elements" |
| 7 | 29-30 | 697-720 | Resolve | Slow pull-back, soft glow | Hold on the active gene | none |

## Non-text signal for "more / less" transcription
Use the visible rate of mRNA streaming away from the gene as the output indicator, plus a soft glow around the TAD that brightens with output. Rates are illustrative and kept low enough that Pol II molecules do not overlap (about 0.4/s baseline, 0.8/s enhancer, 0.13/s repressor, 0.9/s stapler).

## Structures (candidates for a future atomistic version; verify each ID on RCSB before use)
| Element | Candidate PDB | Note |
|---|---|---|
| Nucleosome | 1KX5 | repeated and instanced |
| RNA polymerase II | e.g. 5IY6 (human PIC) | verify chain contents |
| TBP / TFIID | e.g. 1CDW (TBP-DNA) | simplified |
| CTCF (zinc fingers on DNA) | to be chosen | verify |
| Cohesin ring | to be chosen | verify |
| dCas9 + guide RNA | e.g. 4OO8 / 5F9R (SpCas9) | dead form is a mutation of this structure |

Illustrative, not from PDB: loop path, DNA between nucleosomes, enhancer and repressor assemblies, the dCas9-dCas9 fusion, mRNA ribbons, cell and organelles.

## Palette
Slate-blue layered background with depth of field. DNA/chromatin pale lilac and white; Pol II pink; general factors lilac; CTCF/cohesin aqua; enhancer peach; repressor blue-grey; dCas9 stapler teal; mRNA warm peach.

## Render notes
- Rendered with `procedural.py` (2.5D NumPy/Pillow renderer, no Blender). Final cut delivered at 1080p.
- Chromatin jiggle: noise-driven control points, not simulation.

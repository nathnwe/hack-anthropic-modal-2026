# Delivery plan

Read [background.md](background.md) first. Rules for agents and humans are in [../CLAUDE.md](../CLAUDE.md).

## 1. Architecture: one contract, three independent branches
The three branches stay independent through a shared data contract: [`contracts/schema.json`](../contracts/schema.json), with fixtures in `contracts/fixtures/`. One record per gene:

| Section | Content | Written by |
|---|---|---|
| `gene` | symbol, chrom, TSS, disease, direction (`too_much`/`too_little`), evidence (sources) | B |
| `locus` | cell type, TAD, CREs with rE2G scores, E-P contacts | A |
| `staples` | candidate anchor pairs, mode (`up`/`down`/`off`), qualitative predicted delta | A |
| `risks` | per-staple flags (e.g. strong enhancer near an oncogene promoter) | B |
| `ranking` | score 0-100, tier A/B/C, rationale, method recommendation | B |

- Fixtures are **placeholders** (`"illustrative": true`), not real coordinates. Replace them with real, sourced data.
- Validate: `uv run python contracts/validate.py` (CI runs it on every PR).
- **Demo scope:** about 10-15 genes, 2-3 cell types. Use precomputed ENCODE-rE2G predictions and UniversalEPI contact maps where possible; no genome-wide runs. Use Modal only if GPU inference is needed for chosen loci, and cache outputs as small files.

## 2. Operators

| Operator | Branch(es) | Owns | Deliverables |
|---|---|---|---|
| **A: Pipeline** | `feat/pipeline` | `pipeline/` | Chromosome, TAD, CREs, candidate staples with predicted direction (workflow steps 1-4); `locus`, `staples` |
| **B: Scoring & science** | `feat/scoring` | `scoring/` | Dosage target list (Deliverable 1), secondary-risk checks, scoring/tiering, staple-method recommender, rationale text; `gene`, `risks`, `ranking` |
| **C: Front-end & media** | `feat/site`, `feat/animation` | `site/`, `animation/` | Wireframe, Astro dashboard on GitHub Pages, Blender animation, demo script, pitch |

- A is the critical path. If A slips, B and C continue on fixtures.
- B verifies every gene/disease claim against ClinGen, gnomAD (pLI), Open Targets or PubMed; the bio-research MCP connectors can help. This is the most hallucination-prone part.
- C starts the animation first (longest, least dependent). Keep the Blender scene a small parametric `bpy` script so storyboard changes are cheap, and time-box it.

## 3. Models and effort

| Task | Model | Effort | Why |
|---|---|---|---|
| Reading the 2 papers; scoring design; risk logic; method rules | Opus 5.5 (Fable 5.1 for a final review pass) | high | Scientific judgment and correctness |
| Pipeline code, Astro components, Blender scripts, tests | Sonnet 5.5 | medium | Fast once the contract exists |
| Bulk tasks: fixtures, table reformatting, README polish, per-gene rationale text | Haiku 4.5 | low | Cheap, parallel; spot-check output |
| Debugging | Sonnet 5.5, escalate to Opus 5.5 | medium -> high | Escalate after one failed attempt |

**Runtime design:** ranking is deterministic (weighted score + auditable flags). Claude only writes the plain-English rationale and method recommendation, from structured JSON, on Sonnet 5.5. This keeps the demo reproducible and fits the "know what you don't know" judging angle: surface uncertainty as tiers and explicit limitations, not false precision.

## 4. Multi-agent use inside each branch
Each operator runs their own Claude Code session on their own branch (ideally a git worktree).
- **Papers (hour 1, A and B):** two parallel research agents, one per paper, each returning a 1-page brief (data formats, download locations, meaning of outputs, caveats). Commit the briefs to `docs/briefs/`.
- **Implementation:** `create-plan` then `implement-plan` for A's dependent phases.
- **Review:** `/code-review` before each merge; B also does a scientific sanity pass on A's outputs (plausible TAD sizes, correct CRE types).
- **Avoid** several agents editing the same files. Divide by file ownership, not task.

## 5. Branch and merge protocol
- `main` is PR-only; one other human approves. Never push straight to `main`.
- Ownership as in section 2; `contracts/` changes need all three operators to agree.
- Rebase on `main` at least every 2 hours.
- GitHub Pages deploys from `main` via `.github/workflows/pages.yml` once `site/package.json` exists. Data lives in `site/public/data/`.
- Large raw data goes in `data/raw/` (git-ignored); commit only small derived JSON.

## 6. Timeline (compress to fit the event)
1. **Hour 0-1, together:** confirm the track; agree contract and demo genes/cell types; paper briefs; C scaffolds Astro and activates the Pages deploy.
2. **Hour 1-5, parallel:** A builds real data for demo genes; B builds target list, scoring and risk rules on fixtures; C builds dashboard on fixtures and storyboards the animation.
3. **Hour 5-7, integration:** swap fixtures for real outputs; tune scoring; wire the TAD/contact visual.
4. **Hour 7-9, polish:** hero-gene end-to-end demo; "known limitations" panel; finish render; README, 2-minute demo video, pitch (1:30 pitch + 1:30 live demo).

## 7. Risks
- **Scope creep on models:** hold the precomputed-data approach and say so in the pitch.
- **Fabricated biology:** every claim has a source link in the JSON.
- **Predicted effect size is the weakest link:** present rank/confidence tiers, not precise numbers.
- **Animation overrun:** a short clean render beats an ambitious unfinished one.
- **Event rules:** "no prior commits" (see background.md). Resolve this with the organisers early.

# Task A: Pipeline (prompt for Opus 5.5, high effort)

Paste everything below the line into a Claude Code session on branch `feat/pipeline`, in a git worktree if you can.

---

You are Operator A on a 3-person hackathon team building a **target-shortlisting tool for 3D genome engineering**: given a dosage-sensitive gene and a cell type, propose which two genomic regions to "staple" (bring together with programmable CRISPR-based tools) or unstaple so that the gene's expression goes up, down or off, using existing regulatory machinery. You own the **pipeline**: the data and modelling layer that produces the `locus` and `staples` sections of each per-gene record. Teammate B consumes your output for risk checks and scoring; teammate C consumes it for a dashboard.

## Read first, in this order
1. `docs/background.md`, `docs/plan.md`, `CLAUDE.md` (rules, ownership, honesty requirements).
2. `contracts/schema.json` and `contracts/fixtures/*.json`. This is the **binding output contract**. The fixtures are placeholders (`"illustrative": true`), not real data.
3. The two papers in `docs/`: UniversalEPI (Grover et al. 2026) and ENCODE-rE2G (Gschwind et al. 2026). Read them properly before designing anything. Extract: what each model takes as input and outputs, genome build, resolution, what cell types/biosamples are covered, how to obtain precomputed predictions or weights, licences, and stated limitations. Write a 1-page brief for each in `docs/briefs/` (`universalepi.md`, `encode-re2g.md`) and commit it. Cite sections or pages. If the PDFs don't state something, say so rather than guessing.

## Your scope (workflow steps 1-4)
For a given gene symbol and cell type, produce:
1. **Chromosome and TSS** of the gene (state the genome build; use the same build as the models).
2. **TAD** containing the gene: coordinates and boundaries, with the source of the TAD call and its resolution.
3. **CREs** regulating the gene: from ENCODE-rE2G predictions, with each element's coordinates, type and rE2G score, restricted to the TAD or the model's window.
4. **E-P contacts** from UniversalEPI (or precomputed contact maps), in the schema's `contacts` format.
5. **Candidate staples**: ranked pairs of anchor regions (`anchor_a`, `anchor_b`) with a `mode` of `up`, `down` or `off`, and a **qualitative** `predicted_delta` tier (e.g. `weak_/moderate_/strong_` + direction). Define explicit, documented heuristics, for example:
   - `up`: strengthen contact between the promoter and a high-scoring enhancer that currently contacts it weakly.
   - `down`/`off`: tether the promoter toward a silencer or insulator, or sequester a dominant enhancer away from the promoter.
   Keep staples within the TAD unless you can justify otherwise. Record the reasoning for each candidate in a human-readable form.

Out of scope: secondary-risk checks, scoring/tiering, method recommendation (B); dashboard and animation (C). Do not write the `gene`, `risks`, or `ranking` sections, beyond what you need to validate the record. You may emit placeholder values there so records validate, clearly marked.

## Constraints and approach
- **Demo scope:** about 10-15 genes and 2-3 cell types. Do not run either model genome-wide. Prefer **precomputed** ENCODE-rE2G predictions and precomputed or lightweight-inference UniversalEPI outputs. Use **Modal** (GPU) only if inference on the chosen loci is genuinely required and cannot be avoided, and cache every result as a small file.
- Start with the hero genes to match the target list B is building. Until B publishes it, use PMP22 (too much), SHANK3 (too little) and MYC (too much, oncogene) as test genes.
- Python 3.12 via `uv`. Add dependencies with `uv add`. Keep raw downloads in `data/raw/` (git-ignored); commit only small derived JSON. Log every data source, version and download URL in `pipeline/DATA_SOURCES.md`, and add credits to the README.
- **Output location:** write per-gene records to `site/public/data/<symbol>.json`. Every record must pass `uv run python contracts/validate.py`. If you need a schema change, stop and propose it in a PR description for B and C. Don't edit the schema unilaterally.
- **Honesty:** no invented coordinates, scores or citations. Every record keeps provenance (which model, which file, which build, which cell type). Where a gene lacks coverage in the chosen cell type, emit an explicit "no prediction available" result instead of a plausible-looking guess. Predictions are in silico and unvalidated; use tiers, not false precision.
- **Structure:** `pipeline/` as a small package with clear stages (`resolve_gene`, `get_tad`, `get_cres`, `get_contacts`, `propose_staples`) plus one CLI entry point (`uv run python -m pipeline.run --gene PMP22 --cell-type <name>`). Deterministic and re-runnable. Add tests for the staple heuristics and a validation test over every output record.

## How to work
1. Plan first (the session is run with the `plan` skill, so don't invoke `create-plan`): after reading the papers, write a short plan covering data availability, the exact files you will use, the staple heuristics, and the risks, and show it to me before implementing. Flag anything that makes the precomputed-data approach infeasible early; that is the biggest risk to the whole project.
2. Implement in phases, one for each stage, committing as you go on `feat/pipeline`. After the first end-to-end record for one gene, open a PR to `main` so B and C get real data early, then continue with the rest.
3. Use parallel sub-agents for read-only research (paper briefs, data-source discovery). Don't have several agents edit the same files.
4. Before each PR, run `/code-review` and `uv run python contracts/validate.py`.

## Definition of done
- Briefs for both papers are committed in `docs/briefs/`.
- Valid records for at least 3 genes, then 10-15, in at least 2 cell types, each with provenance.
- A documented, tested staple-proposal heuristic, and a README section in `pipeline/` explaining how to reproduce a record from scratch.
- A short "limitations and assumptions" note for the pitch (what the predictions can and cannot tell you).

Start by reading the files above and the two papers, then give me the plan.

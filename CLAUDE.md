# Project rules for agents and humans
Target-shortlisting tool for 3D genome engineering. Read `docs/background.md` then `docs/plan.md`.

- Contract first: all data exchanged between streams conforms to `contracts/schema.json`. Run `uv run python contracts/validate.py` before every PR. Changing the schema needs all 3 operators to agree.
- File ownership: A -> `pipeline/`; B -> `scoring/`; C -> `site/`, `animation/`; `contracts/` shared; `docs/` anyone (small PRs).
- Branches: `feat/pipeline`, `feat/scoring`, `feat/site`, `feat/animation`. PR to `main`, one other human approves. Rebase on main at least every 2h. Never push to main directly.
- Honesty: every gene/disease/locus claim carries a source. Fixtures with `"illustrative": true` are placeholders, never present them as real. Predictions are in silico, unvalidated; show uncertainty as tiers, not false precision.
- Ranking stays deterministic and auditable; LLMs only write rationale text from structured input.
- Python 3.12 via uv. Large raw data goes in `data/raw/` (git-ignored); commit only small derived JSON.
- Build during the event, credit all libraries/models/data used (see README credits).

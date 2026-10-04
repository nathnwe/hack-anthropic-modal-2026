"""One JSON for the web team: the seven demo results plus whatever the screen finished.

    PYTHONPATH=. python3 -m pipeline.export_web      # writes data/derived/website_results.json

NOT the contracts/schema.json format (which has qualitative tiers only and needs all three
operators to agree before it changes). This file is a plain results dump; docs/results-report.md
says what each field means and which caveats must travel with it.
"""

from __future__ import annotations

import json

from pipeline.real_run import DERIVED

DEMO = [("ldlr", "LDLR", "Familial hypercholesterolaemia", "raise"),
        ("nf1", "NF1", "Neurofibromatosis type 1", "raise"),
        ("mecp2_up", "MECP2", "Rett syndrome (loss of function)", "raise"),
        ("mecp2_down", "MECP2", "MECP2 duplication syndrome (gain)", "lower"),
        ("myc_down", "MYC", "Oncogene, overexpressed in cancer", "lower"),
        ("runx1", "RUNX1", "Familial platelet disorder with AML predisposition", "raise"),
        ("nsd1", "NSD1", "Sotos syndrome", "raise")]
KEEP = ("rank", "enhancer_kb", "anchor_a_kb", "anchor_b_kb", "separation_kb", "delta_target", "delta_at_p99",
        "sweep", "sign_stable", "confidence", "risk", "unflagged_collateral_over_tau", "q")


def design(r: dict) -> dict:
    return {k: r.get(k) for k in KEEP}


def main() -> None:
    demo = []
    for key, gene, disease, goal in DEMO:
        d = json.loads((DERIVED / f"demo_{key}.json").read_text())
        rows = d["ranked"] + d["best_below_floor"]
        gated = [design(r) for r in rows if r["risk"]["tier"] == "BLOCKED"
                 and (r["delta_target"] if goal == "raise" else -r["delta_target"]) >= 0.10]
        demo.append({
            "gene": gene, "disease": disease, "goal": goal, "window": d["window"],
            "target_tpm_k562": rows[0]["target_tpm"] if rows else None,
            "flagged_genes_in_window": d["flagged_genes_in_window"],
            "n_designs_scored": d["n_staples_scored"], "n_meeting_dose_goal": d["n_meeting_dose_goal"],
            "designs": [design(r) for r in d["ranked"]],
            "rejected_by_safety_gate": gated,
        })
    screen = json.loads((DERIVED / "screen_results.json").read_text())
    out = {
        "cell_line": "K562 (only cell type with Hi-C here)", "assembly": "GRCh38",
        "status": "in silico, unvalidated; real K562 data; windows for the demo chosen post hoc",
        "demo": demo,
        "screen": {"complete": False, "note": "stopped early; random order, seed 0 (see docs/results-report.md)",
                   "records": [{k: r.get(k) for k in ("gene", "mode", "window", "target_tpm", "flagged_in_window",
                                                      "n_scored", "n_meeting_goal", "best", "best_below_floor", "error")}
                               for r in screen]},
    }
    (DERIVED / "website_results.json").write_text(json.dumps(out, indent=1))
    print("wrote", DERIVED / "website_results.json", f"({len(demo)} demo targets, {len(screen)} screened genes)")


if __name__ == "__main__":
    main()

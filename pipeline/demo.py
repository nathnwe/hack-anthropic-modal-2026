"""Pipeline demo on live K562 data: top staple targets, each with a risk score.

    SDKROOT=<older macOS SDK> uv run --with hic-straw --with pybigtools --with numpy \
        python -m pipeline.demo [window ...]

Two kinds of window, because the safety result differs:

  null      no ClinGen dosage-sensitive gene in the window. Reported as a result:
            to the knowledge of the curated list, nothing here is dangerous to
            disturb (risk tier CLEAR).
  flagged   one or more dosage-sensitive genes in the window, so the risk score
            has something to measure and the gate can fire.

Hi-C, DNase, RNA-seq and the gene models are real (ENCODE, UCSC). The risk list
is ClinGen only: COSMIC and DepMap could not be obtained.
"""

from __future__ import annotations

import sys

from pipeline.real_run import DERIVED, Window, run_window

# (window, curated target genes from docs/targets.md; goal: raise dose)
WINDOWS = {
    "ldlr": (Window("LDLR · familial hypercholesterolaemia · FLAGGED neighbours (SMARCA4, PRKCSH)",
                    "chr19", 10_110_000, 12_110_000), ["LDLR"]),
    "nsd1": (Window("NSD1 · Sotos syndrome · FLAGGED neighbour (DDX41)",
                    "chr5", 176_215_000, 178_215_000), ["NSD1"]),
    "runx1": (Window("RUNX1 · familial platelet disorder · NULL (no other flagged gene)",
                     "chr21", 33_920_000, 35_920_000), ["RUNX1"]),
    "nf1": (Window("NF1 · neurofibromatosis type 1 · NULL (no other flagged gene)",
                   "chr17", 30_235_000, 32_235_000), ["NF1"]),
}


def summary(results: list[dict]) -> str:
    """One row per target: the best surviving design, its risk score, and what was rejected."""
    lines = ["| target | window | best design (enhancer -> promoter) | dT/T | risk | score | confidence | rejected by gate |",
             "|---|---|---|---|---|---|---|---|"]
    for out in results:
        for t in out["targets"]:
            mine = [r for r in out["ranked"] if r["target"] == t]
            rej = sum(1 for r in out["best_below_floor"] + out["ranked"] if r["target"] == t and r["risk"]["tier"] == "BLOCKED")
            if mine:
                r = mine[0]
                lines.append(f"| **{t}** | {out['window']} | {r['enhancer_kb']:,} kb -> {r['anchor_b_kb']:,} kb ({r['separation_kb']} kb) "
                             f"| {r['delta_target']:+.1%} | {r['risk']['tier']} | {r['risk']['score']:.1%} | {r['confidence']} | {out['dropped']['blocked_by_risk']} |")
            else:
                lines.append(f"| **{t}** | {out['window']} | none reaches the dose goal | | | | | {out['dropped']['blocked_by_risk']} |")
    return "\n".join(lines)


def main(names: list[str] | None = None) -> None:
    results = []
    for name in names or list(WINDOWS):
        win, targets = WINDOWS[name]
        results.append(run_window(win, f"demo_{name}", targets))
    table = summary(results)
    print(f"\n{'=' * 100}\nSUMMARY: top staple design per target, with risk score\n{'=' * 100}\n{table}\n")
    (DERIVED / "demo_summary.md").write_text(table + "\n")


if __name__ == "__main__":
    main(sys.argv[1:] or None)

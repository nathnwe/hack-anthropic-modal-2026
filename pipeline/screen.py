"""Screen a large gene set on live K562 data: one 2 Mb window per gene, one best design each.

    SDKROOT=<older macOS SDK> PYTHONPATH=. uv run --with hic-straw --with pybigtools --with numpy \
        python -m pipeline.screen [--shard I/N] [--set up|down|all]

Gene sets (all curated elsewhere, none chosen by us after looking at results):
  up    every ClinGen gene with haploinsufficiency score 3 (422 genes; "too little",
        so the goal is to RAISE dose). Windows come from data/derived/demo_window_scan.json.
  down  the ClinGen triplosensitivity-3 genes (LMNB1, PLP1, TPSAB1) plus MECP2
        (Xq28 duplication syndrome) and MYC (oncogene), goal: LOWER dose.

Resumable: a finished gene is skipped. One compact JSON per gene lands in
data/derived/screen/, and `python -m pipeline.screen_report` joins them.
"""

from __future__ import annotations

import json
import random
import sys
import time
import traceback

from pipeline.real_run import DERIVED, Window, run_window

SCREEN = DERIVED / "screen"
SCAN = DERIVED / "demo_window_scan.json"

# Down-regulation targets. Sources: ClinGen dosage map (TS score 3: LMNB1, PLP1, TPSAB1);
# MECP2 duplication syndrome (Xq28; ClinGen region TS); MYC (oncogene; MYC enhancers were mapped in
# K562 by Fulco et al. 2016 Science, doi:10.1126/science.aag2445).
DOWN_WINDOWS = {
    "MECP2": ("chrX", 153_060_000, 155_060_000),
    "MYC": ("chr8", 127_000_000, 129_000_000),
}


def up_windows(expressed_only: bool = False) -> dict[str, tuple[str, int, int]]:
    out = {}
    for r in json.loads(SCAN.read_text()):
        if expressed_only and not (r["tpm"] or 0) >= 1:   # K562 is the only cell with Hi-C: a gene it does not express says nothing
            continue
        chrom, rng = r["win"].split(":")
        a, b = (int(x) for x in rng.split("-"))
        out[r["gene"]] = (chrom, a, b)
    return out


def compact(gene: str, direction: str, out: dict, seconds: float) -> dict:
    """The one-line story of a window: best design, its risk, and what was rejected."""
    mine = [r for r in out["ranked"]]
    below = out["best_below_floor"]
    best = mine[0] if mine else None
    anyrow = best or (below[0] if below else None)
    sp = out["anchor_space"]
    return {
        "gene": gene,
        "mode": direction,
        "window": out["window"],
        "seconds": round(seconds),
        "flagged_in_window": sorted(out["flagged_genes_in_window"]),
        "target_tpm": anyrow["target_tpm"] if anyrow else None,
        "baseline_T": anyrow["baseline_T"] if anyrow else None,
        "permitted_anchors": sp["permitted_anchors"],
        "genes_expressed": sp["genes_expressed"],
        "genes": sp["genes"],
        "excluded_flagged": sp["excluded_flagged"],
        "n_scored": out["n_staples_scored"],
        "n_meeting_goal": out["n_meeting_dose_goal"],
        "dropped": out["dropped"],
        "best": best,
        "best_below_floor": None if best or not below else below[0],
        "ranked_top5": mine[:5],
    }


def run_one(gene: str, direction: str, win: tuple[str, int, int]) -> dict:
    SCREEN.mkdir(parents=True, exist_ok=True)
    f = SCREEN / f"{direction}_{gene}.json"
    if f.exists():
        return json.loads(f.read_text())
    t = time.time()
    chrom, a, b = win
    try:
        out = run_window(
            Window(f"{gene} ({direction})", chrom, a, b),
            f"screen_{direction}_{gene}",
            [gene],
            write=False,
            direction=direction,
            verbose=False,
        )
        rec = compact(gene, direction, out, time.time() - t)
    except (
        Exception
    ) as e:  # a window that cannot be built is recorded, never skipped silently
        rec = {
            "gene": gene,
            "mode": direction,
            "window": f"{chrom}:{a}-{b}",
            "error": f"{type(e).__name__}: {e}",
            "trace": traceback.format_exc(limit=3),
        }
    f.write_text(json.dumps(rec, indent=1))
    return rec


def main(argv: list[str]) -> None:
    shard, n_shards = 0, 1
    which = "all"
    for i, a in enumerate(argv):
        if a == "--shard":
            shard, n_shards = (int(x) for x in argv[i + 1].split("/"))
        if a == "--set":
            which = argv[i + 1]
    expressed_only = "--expressed" in argv
    jobs: list[tuple[str, str, tuple[str, int, int]]] = []
    if which in ("down", "all"):
        wins = dict(DOWN_WINDOWS)
        wins.update({"LMNB1": ("chr5", 125_800_000, 127_800_000),
                     "PLP1": ("chrX", 102_780_000, 104_780_000),
                     "TPSAB1": ("chr16", 240_000, 2_240_000)})
        jobs += [(g, "down", w) for g, w in wins.items()]
    if which in ("up", "all"):
        ups = sorted(up_windows(expressed_only).items())
        random.Random(0).shuffle(ups)  # fixed seed: a run stopped early is a random sample, not the A-C genes
        jobs += [(g, "up", w) for g, w in ups]
    jobs = [j for i, j in enumerate(jobs) if i % n_shards == shard]
    for k, (g, d, w) in enumerate(jobs, 1):
        rec = run_one(g, d, w)
        b = rec.get("best")
        tag = (
            "ERROR " + rec["error"][:80]
            if "error" in rec
            else (
                f"{b['delta_target']:+.1%} {b['risk']['tier']}"
                if b
                else "no design reaches the goal"
            )
        )
        print(
            f"[{shard}/{n_shards}] {k}/{len(jobs)} {d:<4} {g:<10} {rec.get('seconds', '-'):>5}s  {tag}",
            flush=True,
        )


if __name__ == "__main__":
    main(sys.argv[1:])

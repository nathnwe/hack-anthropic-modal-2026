"""Join the per-gene screen records into one table and the headline counts.

    PYTHONPATH=. python3 -m pipeline.screen_report

Writes data/derived/screen_results.json (every gene, compact) and
data/derived/screen_tables.md (the tables the research report quotes).
Pure bookkeeping: no new modelling happens here.
"""

from __future__ import annotations

import json
from collections import Counter

from pipeline.screen import SCREEN
from pipeline.real_run import DERIVED


def load() -> list[dict]:
    return [json.loads(f.read_text()) for f in sorted(SCREEN.glob("*.json"))]


def status(r: dict) -> str:
    if "error" in r:
        return "error"
    if r["best"]:
        return "design"
    return "none_reaches_goal"


def row(r: dict) -> str:
    b = r["best"]
    k = b["risk"]
    other = f"{k['n_flagged_in_window'] - (1 if b['target_flag'] else 0)}"
    return (f"| {r['gene']} | {r['target_tpm']:.0f} | {b['enhancer_kb']:,} -> {b['anchor_b_kb']:,} ({b['separation_kb']} kb) "
            f"| {b['delta_target']:+.0%} | {k['tier']} ({k['score']:.1%}) | {b['confidence'].lower()} | {other} |")


HEAD = ("| gene | TPM | enhancer -> promoter-side anchor (kb) | dT/T | risk | confidence | other flagged genes in window |\n"
        "|---|---|---|---|---|---|---|")


def main() -> None:
    recs = load()
    (DERIVED / "screen_results.json").write_text(json.dumps(recs, indent=1))
    lines = []
    for mode in ("up", "down"):
        rs = [r for r in recs if r["mode"] == mode]
        c = Counter(status(r) for r in rs)
        expr = [r for r in rs if "error" not in r and (r["target_tpm"] or 0) >= 1] \
            + [r for r in rs if "error" not in r and r["target_tpm"] is None]
        lines.append(f"## {mode.upper()} screen: {len(rs)} genes · {dict(c)}\n")
        designs = [r for r in rs if status(r) == "design"]
        tiers = Counter(r["best"]["risk"]["tier"] for r in designs)
        conf = Counter(r["best"]["confidence"] for r in designs)
        lines.append(f"best-design risk tiers: {dict(tiers)}; confidence: {dict(conf)}\n")
        ex = [r for r in designs if (r["target_tpm"] or 0) >= 1]
        lines.append(f"genes expressed in K562 (TPM>=1) with a design: {len(ex)} of "
                     f"{sum(1 for r in rs if (r['target_tpm'] or 0) >= 1 or r['target_tpm'] is None)} expressed\n")
        for title, sel in (("CLEAR and CONFIDENT, expressed", [r for r in ex if r["best"]["risk"]["tier"] == "CLEAR" and r["best"]["confidence"] == "CONFIDENT"]),
                           ("LOW / ELEVATED, expressed", [r for r in ex if r["best"]["risk"]["tier"] in ("LOW", "ELEVATED")])):
            sel = sorted(sel, key=lambda r: -abs(r["best"]["delta_target"]))
            lines.append(f"### {mode}: {title} ({len(sel)})\n\n{HEAD}")
            lines += [row(r) for r in sel[:40]]
            lines.append("")
        errs = [r for r in rs if "error" in r]
        if errs:
            lines.append(f"errors ({len(errs)}): " + "; ".join(f"{r['gene']}: {r['error'][:60]}" for r in errs[:10]) + "\n")
    (DERIVED / "screen_tables.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:60]))


if __name__ == "__main__":
    main()

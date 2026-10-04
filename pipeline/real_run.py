"""The method, run on real K562 data at the beta-globin locus.

    SDKROOT=<older macOS SDK> uv run --with hic-straw --with pybigtools \
        --with numpy python -m pipeline.real_run

Prints a report and writes data/derived/globin_k562.json (small, committable).

Two things are checked before any staple is simulated:

  sanity   Does baseline ABC recover known globin biology? The locus control
           region (LCR, ~chr11:5.28 Mb) should dominate the regulatory input
           of the embryonic/fetal globin genes. If it does not, nothing
           downstream is worth reading.
  gate     Which genes in the window are flagged at all.

Then every (target gene, strong enhancer) staple is simulated and ranked.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pipeline import abc as abc_mod
from pipeline import anchors as anchors_mod
from pipeline import sweep as sweep_mod
from pipeline.locus import BIN_KB
from pipeline.real_data import load_k562_locus
from scoring.rank import Candidate, rank
from scoring.risk import TAU, concerns

CHROM, START, END = "chr11", 4_700_000, 6_700_000
LCR_KB = (5_270, 5_310)  # LCR hypersensitive sites, hg38
GLOBIN = ["HBE1", "HBG2", "HBG1", "HBD", "HBB"]
TOP_ENHANCERS = 8
EXPRESSED_TPM = 1.0
MIN_SEPARATION_BINS = 6  # 30 kb: a staple shorter than this is not a staple
ANCHOR_REACH = 8  # search +/-40 kb for a permitted anchor
OUT = Path(__file__).resolve().parents[1] / "data" / "derived" / "globin_k562.json"


def kb(bin_: int) -> int:
    return START // 1000 + bin_ * BIN_KB


def nearest_permitted(accessible, forbidden, target_bin):
    ok = anchors_mod.candidate_anchors(accessible, forbidden, target_bin, ANCHOR_REACH)
    return min(ok, key=lambda b: abs(b - target_bin)) if ok else None


def staple_effect(loc, order, gidx, t0, a, b, target):
    """Simulate one staple; return dT/T for every gene, the headline q and sweep stability."""
    sep = abs(a - b)
    q = sweep_mod.yardstick(loc.contact, sep)
    totals = lambda m: dict(zip(order, abc_mod.total_input(loc.activities, m, loc.elements, gidx)))
    delta = {x: (v - t0[x]) / t0[x] for x, v in totals(sweep_mod.apply_staple(loc.contact, a, b, q)).items()}
    sweep = []
    for qi in sweep_mod.q_grid(loc.contact, sep):
        ti = totals(sweep_mod.apply_staple(loc.contact, a, b, qi))
        sweep.append((ti[target] - t0[target]) / t0[target])
    return delta, q, bool(all(d > 0 for d in sweep)), sweep


def directed(loc, order, gidx, t0, forbidden, accessible, label):
    """LCR -> each globin promoter. Deng et al. 2012 Cell / 2014 Nature forced
    exactly this loop (ZF-Ldb1) and activated HBB, so the sign is checkable."""
    lcr_bin = int(np.argmax(np.where((loc.elements * BIN_KB + START // 1000 >= LCR_KB[0]) &
                                     (loc.elements * BIN_KB + START // 1000 <= LCR_KB[1]), loc.activities, -1)))
    lcr_bin = int(loc.elements[lcr_bin])
    a = nearest_permitted(accessible, forbidden, lcr_bin)
    print(f"\nDIRECTED — staple the LCR ({kb(lcr_bin)} kb) to each globin promoter  [{label}]")
    print(f"  anchor A at {kb(a) if a is not None else None} kb")
    rows = []
    for g in GLOBIN:
        b = nearest_permitted(accessible, forbidden, loc.genes[g])
        if a is None or b is None or abs(a - b) < 2:
            print(f"  {g:<5} no permitted anchor"); continue
        delta, q, stable, sweep = staple_effect(loc, order, gidx, t0, a, b, g)
        coll = sorted(((x, d) for x, d in delta.items() if x != g and abs(d) > TAU), key=lambda kv: -abs(kv[1]))
        print(f"  {g:<5} anchor B {kb(b)} kb  sep {abs(a-b)*BIN_KB:>4} kb  dT/T {delta[g]:+8.1%}  "
              f"sweep p50->p99 {sweep[0]:+.1%}..{sweep[-1]:+.1%}  moved>{TAU:.0%}: "
              + (", ".join(f"{x} {d:+.0%}" for x, d in coll[:4]) or "—"))
        rows.append({"target": g, "anchor_a_kb": kb(a), "anchor_b_kb": kb(b), "separation_kb": abs(a - b) * BIN_KB,
                     "delta_target": round(delta[g], 4), "sweep": [round(x, 4) for x in sweep],
                     "collateral": [{"gene": x, "delta": round(d, 4)} for x, d in coll[:8]],
                     "all_globin": {x: round(delta[x], 4) for x in GLOBIN if x in delta}})
    return rows


def main() -> dict:
    loc, prov = load_k562_locus(CHROM, START, END)
    order = loc.gene_order
    gidx = np.array([loc.genes[g] for g in order])
    n = len(loc.contact)

    # ---- baseline ---------------------------------------------------------
    t0 = dict(
        zip(order, abc_mod.total_input(loc.activities, loc.contact, loc.elements, gidx))
    )
    share = abc_mod.abc_share(loc.activities, loc.contact, loc.elements, gidx)
    gpos = {g: i for i, g in enumerate(order)}

    print(f"\n{prov['window']}  ·  {prov['counts']}\n")
    print("SANITY CHECK — top ABC element per globin gene (expect the LCR, ~5,280 kb)")
    sanity = {}
    for g in GLOBIN:
        if g not in gpos:
            continue
        col = share[:, gpos[g]]
        top = int(np.argmax(col))
        e_kb = kb(int(loc.elements[top]))
        in_lcr = LCR_KB[0] <= e_kb <= LCR_KB[1]
        sanity[g] = {
            "top_element_kb": e_kb,
            "share": round(float(col[top]), 3),
            "is_lcr": in_lcr,
        }
        print(
            f"  {g:<5} top element at {e_kb} kb, ABC share {col[top]:.2f}   {'LCR ✓' if in_lcr else 'not the LCR ✗'}"
        )

    # An enhancer is protected only if it drives a gene that is actually on.
    # Without this, 56 of 57 elements are 'functional' for some unexpressed
    # olfactory receptor and no anchor survives.
    on = np.array([loc.expression.get(g, 0.0) >= EXPRESSED_TPM for g in order])
    functional = np.where((share[:, on] >= abc_mod.ABC_THRESHOLD).any(axis=1))[0]
    functional_bins = loc.elements[functional]
    print(f"{int(on.sum())} of {len(order)} genes expressed in K562 (TPM >= {EXPRESSED_TPM:g})")
    print(f"{len(functional_bins)} of {len(loc.elements)} candidate elements pass ABC >= {abc_mod.ABC_THRESHOLD}"
          " for an expressed gene -> excluded as anchors")

    forbidden = anchors_mod.forbidden_mask(
        n, functional_bins, loc.genes.values(), loc.gene_bodies.values(), loc.ctcf
    )
    permitted = int((loc.accessible & ~forbidden).sum())
    print(f"ATAC peaks {int(loc.accessible.sum())} -> permitted anchors {permitted}")
    print(
        f"flagged (ClinGen) genes in window: {len(loc.flags)}  -> the safety gate cannot fire here\n"
    )

    # ---- every staple: (strong enhancer) x (target gene) ------------------
    strongest = np.argsort(-loc.activities)[:TOP_ENHANCERS]
    candidates, rows = [], []
    for ei in strongest:
        e_bin = int(loc.elements[ei])
        a = nearest_permitted(loc.accessible, forbidden, e_bin)
        if a is None:
            continue
        for g in order:
            b = nearest_permitted(loc.accessible, forbidden, loc.genes[g])
            if b is None or abs(a - b) < MIN_SEPARATION_BINS:
                continue
            sep = abs(a - b)
            q = sweep_mod.yardstick(loc.contact, sep)
            after = sweep_mod.apply_staple(loc.contact, a, b, q)
            t1 = dict(
                zip(
                    order,
                    abc_mod.total_input(loc.activities, after, loc.elements, gidx),
                )
            )
            delta = {x: (t1[x] - t0[x]) / t0[x] for x in order}

            grid = sweep_mod.q_grid(loc.contact, sep)
            tops = []
            for qi in grid:
                ti = dict(
                    zip(
                        order,
                        abc_mod.total_input(
                            loc.activities,
                            sweep_mod.apply_staple(loc.contact, a, b, qi),
                            loc.elements,
                            gidx,
                        ),
                    )
                )
                tops.append(g if (ti[g] - t0[g]) / t0[g] > 0 else None)
            stable = all(t == g for t in tops)

            found = concerns(delta, loc.flags, TAU, target=g)
            collateral = sorted(
                ((x, d) for x, d in delta.items() if x != g and abs(d) > TAU),
                key=lambda kv: -abs(kv[1]),
            )
            cand = Candidate(
                name=f"{g} <- {kb(e_bin)}kb",
                delta_t=delta[g],
                provenance="hic",
                distance_kb=sep * BIN_KB,
                abstained=False,
                rank_stable=stable,
                safe=not found,
            )
            candidates.append(cand)
            rows.append(
                {
                    "target": g,
                    "enhancer_kb": kb(e_bin),
                    "anchor_a_kb": kb(a),
                    "anchor_b_kb": kb(b),
                    "separation_kb": sep * BIN_KB,
                    "q": round(q, 3),
                    "delta_target": round(delta[g], 4),
                    "tier": cand.tier,
                    "safe": cand.safe,
                    "rank_stable": stable,
                    "collateral_over_tau": [
                        {"gene": x, "delta": round(d, 4)} for x, d in collateral
                    ],
                }
            )

    relaxed_ok = np.ones(n, dtype=bool)          # SPEC DEVIATION: ignore accessibility, keep every other exclusion
    d_strict = directed(loc, order, gidx, t0, forbidden, loc.accessible, "strict step 3: accessible and not forbidden")
    d_relaxed = directed(loc, order, gidx, t0, forbidden, relaxed_ok, "RELAXED: any non-forbidden bin (spec deviation)")

    ranked = rank(candidates)
    print(
        f"{len(candidates)} staples simulated, {len(ranked)} survive the gate. Top 12 by dT/T:\n"
    )
    print(
        f"  {'target':<9}{'enhancer':>9}{'sep kb':>8}{'dT/T':>9}  {'tier':<10}collateral genes moved > {TAU:.0%}"
    )
    by_name = {c.name: r for c, r in zip(candidates, rows)}
    for c in ranked[:12]:
        r = by_name[c.name]
        coll = (
            ", ".join(
                f"{x['gene']} {x['delta']:+.0%}" for x in r["collateral_over_tau"][:3]
            )
            or "—"
        )
        print(
            f"  {r['target']:<9}{r['enhancer_kb']:>8}k{r['separation_kb']:>8}{r['delta_target']:>+9.1%}  {c.tier:<10}{coll}"
        )

    globin_rows = sorted(
        (
            r
            for r in rows
            if r["target"] in GLOBIN and r["enhancer_kb"] in range(*LCR_KB)
        ),
        key=lambda r: -r["delta_target"],
    )
    result = {
        "provenance": prov,
        "sanity_abc_globin": sanity,
        "n_functional_elements": int(len(functional_bins)),
        "n_permitted_anchors": permitted,
        "n_flagged_genes": len(loc.flags),
        "n_staples": len(rows),
        "staples": sorted(rows, key=lambda r: -r["delta_target"])[:40],
        "lcr_to_globin": globin_rows,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=1))
    print(f"\nwrote {OUT.relative_to(OUT.parents[2])}")
    return result


if __name__ == "__main__":
    main()

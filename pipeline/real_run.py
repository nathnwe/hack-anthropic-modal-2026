"""The method, run on real K562 data at the beta-globin locus.

    SDKROOT=<older macOS SDK> uv run --with hic-straw --with pybigtools \
        --with numpy python -m pipeline.real_run

Writes data/derived/globin_k562.json and data/derived/lcr_validation.json
(small, committable). Four things happen, in this order:

  1 sanity     Does baseline ABC put the locus control region (LCR, ~5.28 Mb)
               on the globin genes? A consistency check, not an independent
               confirmation: the LCR wins on activity (5.09 vs a window median
               of 0.72), not on contact.
  2 search     Every (strong enhancer x target gene) staple under the strict
               step-3 anchor rule, ranked.
  3 validation LCR -> each globin promoter with the exclusions WAIVED, the way
               the published forced-looping experiments did it (Deng et al.
               2012 Cell; 2014 Nature). A sign-and-competition check only:
               those experiments were not in K562.
  4 sensitivity The same validation across sigma, because the collateral
               readout depends on it far more than the target does.

What this run does NOT exercise, and says so in the output:
  * the safety gate: no ClinGen-flagged gene lies in this window, and COSMIC
    and DepMap could not be obtained, so both gate tests are vacuous here;
  * UniversalEPI: not run, so the abstention test is a constant and
    provenance is constant 'hic'. Two of the four confidence tests are live.
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
from scoring.rank import MAX_CONFIDENT_DISTANCE_KB, MIN_EFFECT, Candidate, rank
from scoring.risk import TAU, concerns

CHROM, START, END = "chr11", 4_700_000, 6_700_000
LCR_KB = (5_270, 5_310)  # LCR hypersensitive sites, hg38
GLOBIN = ["HBE1", "HBG2", "HBG1", "HBD", "HBB"]
TOP_ENHANCERS = 8
EXPRESSED_TPM = 1.0
MIN_SEPARATION_BINS = 6  # 30 kb: a staple shorter than this is not a staple
ANCHOR_REACH = 8  # search +/-40 kb for a permitted anchor
SIGMAS = (1.5, 3.0, 4.5)
DERIVED = Path(__file__).resolve().parents[1] / "data" / "derived"

CONFIG = {
    "TAU": TAU,
    "MIN_EFFECT": MIN_EFFECT,
    "TOP_ENHANCERS": TOP_ENHANCERS,
    "EXPRESSED_TPM": EXPRESSED_TPM,
    "MIN_SEPARATION_BINS": MIN_SEPARATION_BINS,
    "ANCHOR_REACH": ANCHOR_REACH,
    "MAX_CONFIDENT_DISTANCE_KB": MAX_CONFIDENT_DISTANCE_KB,
    "SIGMA_BINS": sweep_mod.SIGMA_BINS,
    "RENORM_ITERS": sweep_mod.RENORM_ITERS,
}
# What the tier actually rests on in this run. Printed beside every row.
TESTS_EVALUATED = {
    "abstention": "NOT RUN: UniversalEPI was not used, constant False",
    "provenance": "constant: measured Hi-C",
    "distance": "live",
    "sign_stability": "live",
}


def kb(bin_: int) -> int:
    return START // 1000 + bin_ * BIN_KB


def nearest_permitted(accessible, forbidden, target_bin):
    ok = anchors_mod.candidate_anchors(accessible, forbidden, target_bin, ANCHOR_REACH)
    return min(ok, key=lambda b: abs(b - target_bin)) if ok else None


class Scorer:
    """Baseline T for one locus, and the effect of any staple on every gene."""

    def __init__(self, loc):
        self.loc = loc
        self.order = loc.gene_order
        self.gidx = np.array([loc.genes[g] for g in self.order])
        self.t0 = self.totals(loc.contact)
        bad = [g for g, v in self.t0.items() if not v > 0]
        if bad:  # dT/T is undefined at T = 0; refuse rather than emit inf
            raise ValueError(f"zero baseline regulatory input for {bad}")

    def totals(self, m):
        t = abc_mod.total_input(self.loc.activities, m, self.loc.elements, self.gidx)
        return dict(zip(self.order, t))

    def delta(self, a, b, q, sigma=sweep_mod.SIGMA_BINS):
        m = sweep_mod.apply_staple(self.loc.contact, a, b, q, sigma)
        return {g: (v - self.t0[g]) / self.t0[g] for g, v in self.totals(m).items()}

    def effect(self, a, b, target, sigma=sweep_mod.SIGMA_BINS):
        q, noop = sweep_mod.staple_strength(self.loc.contact, a, b)
        delta = self.delta(a, b, q, sigma)
        sweep = [
            self.delta(a, b, qi, sigma)[target]
            for qi in sweep_mod.q_grid(self.loc.contact, abs(a - b))
        ]
        return {
            "delta": delta,
            "q": q,
            "noop": noop,
            "sweep": sweep,
            "sign_stable": bool(all(s > 0 for s in sweep)),
        }


def sanity(loc, share, order):
    gpos = {g: i for i, g in enumerate(order)}
    print(
        "SANITY — top ABC element per globin gene (the LCR wins on activity, not contact)"
    )
    out = {}
    for g in GLOBIN:
        col = share[:, gpos[g]]
        top = int(np.argmax(col))
        second = float(np.sort(col)[-2])
        e_kb = kb(int(loc.elements[top]))
        out[g] = {
            "top_element_kb": e_kb,
            "share": round(float(col[top]), 3),
            "margin_over_second": round(float(col[top]) / second, 2),
            "is_lcr": LCR_KB[0] <= e_kb <= LCR_KB[1],
        }
        print(
            f"  {g:<5} {e_kb} kb  share {col[top]:.2f}  margin {col[top] / second:.2f}x  "
            f"{'LCR' if out[g]['is_lcr'] else 'not the LCR'}"
        )
    return out


def search(sc, forbidden):
    loc = sc.loc
    cands, rows, dropped = [], [], {"noop": 0, "no_anchor": 0}
    for ei in np.argsort(-loc.activities)[:TOP_ENHANCERS]:
        e_bin = int(loc.elements[ei])
        a = nearest_permitted(loc.accessible, forbidden, e_bin)
        for g in sc.order:
            b = nearest_permitted(loc.accessible, forbidden, loc.genes[g])
            if a is None or b is None or abs(a - b) < MIN_SEPARATION_BINS:
                dropped["no_anchor"] += 1
                continue
            r = sc.effect(a, b, g)
            if r["noop"]:
                dropped["noop"] += 1
                continue
            found = concerns(r["delta"], loc.flags, TAU, target=g)
            cand = Candidate(
                name=f"{g} <- {kb(e_bin)}kb",
                delta_t=r["delta"][g],
                provenance="hic",
                distance_kb=abs(a - b) * BIN_KB,
                abstained=False,
                sign_stable=r["sign_stable"],
                safe=not found,
            )
            coll = sorted(
                ((x, d) for x, d in r["delta"].items() if x != g and abs(d) > TAU),
                key=lambda kv: -abs(kv[1]),
            )
            cands.append(cand)
            rows.append(
                {
                    "name": cand.name,
                    "target": g,
                    "target_tpm": round(loc.expression.get(g, 0.0), 2),
                    "baseline_T": round(sc.t0[g], 2),
                    "enhancer_kb": kb(e_bin),
                    "anchor_a_kb": kb(a),
                    "anchor_b_kb": kb(b),
                    "separation_kb": abs(a - b) * BIN_KB,
                    "q": round(r["q"], 3),
                    "delta_target": round(r["delta"][g], 4),
                    "sweep": [round(x, 4) for x in r["sweep"]],
                    "sign_stable": r["sign_stable"],
                    "safe": cand.safe,
                    "tier": cand.tier,
                    "collateral_over_tau": [
                        {"gene": x, "delta": round(d, 4)} for x, d in coll
                    ],
                }
            )
    ranked = rank(cands)
    order = {c.name: i + 1 for i, c in enumerate(ranked)}
    for r in rows:
        r["rank"] = order.get(r["name"])
    dropped["unsafe"] = int(sum(not c.safe for c in cands))
    dropped["below_min_effect"] = int(sum(c.safe and c.delta_t < MIN_EFFECT for c in cands))
    return ranked, rows, dropped


def validation(sc, forbidden_unused=None):
    """LCR -> each globin promoter with every exclusion waived."""
    loc = sc.loc
    lcr = int(
        loc.elements[
            np.argmax(
                np.where(
                    (loc.elements * BIN_KB + START // 1000 >= LCR_KB[0])
                    & (loc.elements * BIN_KB + START // 1000 <= LCR_KB[1]),
                    loc.activities,
                    -1,
                )
            )
        ]
    )
    out = {"lcr_kb": kb(lcr), "per_gene": {}, "sigma_sensitivity": {}}
    print(
        f"\nVALIDATION — LCR ({kb(lcr)} kb) -> each globin promoter, exclusions waived (published design)"
    )
    print(
        f"  {'gene':<6}{'TSS kb':>7}{'sep':>5}{'C[a,b] pct':>11}{'TPM':>8}{'dT/T':>9}  note"
    )
    for g in GLOBIN:
        b = loc.genes[g]
        r = sc.effect(lcr, b, g)
        band = np.diagonal(loc.contact, offset=abs(lcr - b))
        pct = float((band < loc.contact[lcr, b]).mean() * 100)
        note = (
            "already the strongest contact at this separation: staple is a NO-OP"
            if r["noop"]
            else ""
        )
        coll = sorted(
            ((x, d) for x, d in r["delta"].items() if x != g and abs(d) > 0.05),
            key=lambda kv: -abs(kv[1]),
        )
        print(
            f"  {g:<6}{kb(b):>7}{abs(lcr - b) * BIN_KB:>5}{pct:>10.0f}%{loc.expression.get(g, 0):>8.0f}"
            f"{r['delta'][g]:>+9.1%}  {note or ', '.join(f'{x} {d:+.0%}' for x, d in coll[:3])}"
        )
        out["per_gene"][g] = {
            "tss_kb": kb(b),
            "sep_kb": abs(lcr - b) * BIN_KB,
            "contact_percentile": round(pct, 1),
            "tpm": round(loc.expression.get(g, 0.0), 2),
            "baseline_T": round(sc.t0[g], 2),
            "q": round(r["q"], 3),
            "noop": r["noop"],
            "delta": round(r["delta"][g], 4),
            "sweep": [round(x, 4) for x in r["sweep"]],
            "collateral": [{"gene": x, "delta": round(d, 4)} for x, d in coll[:6]],
        }
    print("\nSIGMA SENSITIVITY — LCR -> HBB")
    b = loc.genes["HBB"]
    for s in SIGMAS:
        r = sc.effect(lcr, b, "HBB", sigma=s)
        others = sorted(
            ((x, d) for x, d in r["delta"].items() if x != "HBB"),
            key=lambda kv: -abs(kv[1]),
        )
        n2 = int(sum(abs(d) > 0.02 for x, d in r["delta"].items() if x != "HBB"))
        print(
            f"  sigma {s:>3}: HBB {r['delta']['HBB']:+.1%}   largest collateral {others[0][0]} {others[0][1]:+.1%}"
            f"   other genes > 2%: {n2}"
        )
        out["sigma_sensitivity"][str(s)] = {
            "hbb": round(r["delta"]["HBB"], 4),
            "largest_collateral": others[0][0],
            "largest_collateral_delta": round(others[0][1], 4),
            "genes_over_2pct": n2,
        }
    return out


def main() -> dict:
    loc, prov = load_k562_locus(CHROM, START, END)
    sc = Scorer(loc)
    n = len(loc.contact)
    share = abc_mod.abc_share(loc.activities, loc.contact, loc.elements, sc.gidx)
    print(f"\n{prov['window']}  ·  {prov['counts']}\n")
    sanity_out = sanity(loc, share, sc.order)

    on = np.array([loc.expression.get(g, 0.0) >= EXPRESSED_TPM for g in sc.order])
    functional_bins = loc.elements[
        np.where((share[:, on] >= abc_mod.ABC_THRESHOLD).any(axis=1))[0]
    ]
    forbidden = anchors_mod.forbidden_mask(
        n, functional_bins, loc.genes.values(), loc.gene_bodies.values(), loc.ctcf
    )
    permitted = int((loc.accessible & ~forbidden).sum())
    print(
        f"\n{int(on.sum())}/{len(sc.order)} genes expressed; {len(functional_bins)}/{len(loc.elements)} elements "
        f"protected; ATAC peaks {int(loc.accessible.sum())} -> permitted anchors {permitted}"
    )
    print(
        f"ClinGen-flagged genes in window: {len(loc.flags)}  -> BOTH safety-gate tests are vacuous here"
    )
    print(f"Confidence tests: {TESTS_EVALUATED}")

    ranked, rows, dropped = search(sc, forbidden)
    print(
        f"\nSEARCH — {len(rows)} staples scored, {len(ranked)} meet the dose goal (dT/T >= {MIN_EFFECT:.0%}); "
        f"dropped: {dropped}"
    )
    for c in ranked[:12]:
        r = next(x for x in rows if x["name"] == c.name)
        print(
            f"  #{r['rank']:<3}{r['target']:<8}{r['enhancer_kb']:>6}k {r['delta_target']:>+8.1%}  {c.tier}  (TPM {r['target_tpm']})"
        )
    if not ranked:
        print(
            "  no staple in this window meets the dose goal under the strict anchor rule"
        )
    best = max(rows, key=lambda r: r["delta_target"], default=None)
    if best:
        print(
            f"  best available, below the floor: {best['name']} {best['delta_target']:+.1%} (TPM {best['target_tpm']})"
        )

    val = validation(sc)
    DERIVED.mkdir(parents=True, exist_ok=True)
    (DERIVED / "lcr_validation.json").write_text(json.dumps(val, indent=1))
    result = {
        "provenance": prov,
        "config": CONFIG,
        "tests_evaluated": TESTS_EVALUATED,
        "sanity_abc_globin": sanity_out,
        "n_functional_elements": int(len(functional_bins)),
        "n_permitted_anchors": permitted,
        "n_flagged_genes": len(loc.flags),
        "safety_gate": "vacuous: no flagged genes in window; COSMIC/DepMap unavailable",
        "n_staples_scored": len(rows),
        "n_meeting_dose_goal": len(ranked),
        "dropped": dropped,
        "staples": sorted(rows, key=lambda r: -r["delta_target"])[:40],
    }
    (DERIVED / "globin_k562.json").write_text(json.dumps(result, indent=1))
    print("\nwrote data/derived/globin_k562.json, data/derived/lcr_validation.json")
    return result


if __name__ == "__main__":
    main()

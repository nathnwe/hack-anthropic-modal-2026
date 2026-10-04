"""The method on real K562 data, for any 2 Mb window.

    SDKROOT=<older macOS SDK> uv run --with hic-straw --with pybigtools \
        --with numpy python -m pipeline.demo

(`pipeline/demo.py` is the entry point; this module is the engine.)

For one window it: scores every gene's baseline regulatory input T; finds the
permitted staple anchors; simulates every (strong enhancer -> target gene) staple;
and ranks the ones that raise the target by at least MIN_EFFECT, each with a
risk score and a confidence tier.

Anchor rules (all relaxed from the first draft, see docs/method.md step 3):
  * an anchor must be within FLANK_BINS of an ATAC peak, not on one;
  * the promoter exclusion is the TSS bin only;
  * gene bodies are exons, not whole transcripts;
  * enhancers of EXPRESSED genes, and every CTCF motif hit, are excluded.

What a run does NOT exercise, and says so in its output:
  * UniversalEPI is not run, so abstention is a constant and provenance is
    constant 'hic'. Two of the four confidence tests are live;
  * COSMIC and DepMap are unobtainable, so risk rests on ClinGen alone.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from pipeline import abc as abc_mod
from pipeline import anchors as anchors_mod
from pipeline import sweep as sweep_mod
from pipeline.locus import BIN_KB
from pipeline.real_data import load_k562_locus
from scoring.rank import MAX_CONFIDENT_DISTANCE_KB, MIN_EFFECT, Candidate, rank
from scoring.risk import ELEVATED_AT, TAU, risk_score

TOP_ENHANCERS = 10
TARGET_ENHANCERS = 15  # when ranking designs for one named target, consider this many strongest enhancers
EXPRESSED_TPM = 1.0
MIN_SEPARATION_BINS = 6  # 30 kb: a staple shorter than this is not a staple
ANCHOR_REACH = 8  # search +/-40 kb for a permitted anchor
SIGMAS = (1.5, 3.0, 4.5)
DERIVED = Path(__file__).resolve().parents[1] / "data" / "derived"

CONFIG = {
    "TAU": TAU,
    "ELEVATED_AT": ELEVATED_AT,
    "MIN_EFFECT": MIN_EFFECT,
    "TOP_ENHANCERS": TOP_ENHANCERS,
    "EXPRESSED_TPM": EXPRESSED_TPM,
    "MIN_SEPARATION_BINS": MIN_SEPARATION_BINS,
    "ANCHOR_REACH": ANCHOR_REACH,
    "FLANK_BINS": anchors_mod.FLANK_BINS,
    "TSS_PAD_BINS": anchors_mod.TSS_PAD_BINS,
    "MAX_CONFIDENT_DISTANCE_KB": MAX_CONFIDENT_DISTANCE_KB,
    "SIGMA_BINS": sweep_mod.SIGMA_BINS,
    "RENORM_ITERS": sweep_mod.RENORM_ITERS,
    "DOWN": {"DECOY_REACH_BINS": 100, "DECOY_KEEP_OUT_BINS": 12,
             "N_DECOYS": 40, "DOWN_ENHANCERS": 5},
}
TESTS_EVALUATED = {
    "abstention": "NOT RUN: UniversalEPI not used, constant False",
    "provenance": "constant: measured Hi-C",
    "distance": "live",
    "sign_stability": "live",
}


@dataclass(frozen=True)
class Window:
    label: str
    chrom: str
    start: int
    end: int

    @property
    def start_kb(self) -> int:
        return self.start // 1000

    def kb(self, bin_: int) -> int:
        return self.start_kb + int(bin_) * BIN_KB


class Scorer:
    """Baseline T for one locus, and the effect of any staple on every gene."""

    def __init__(self, loc):
        self.loc = loc
        full = loc.gene_order
        self.gidx = np.array([loc.genes[g] for g in full])
        self.order = full
        t = self.totals(loc.contact)
        # A gene with no contact data (T = 0, e.g. an unmappable stretch) has an undefined dT/T.
        # Exclude it, but never silently: it is recorded, and a flagged one is a hole in the risk score.
        self.excluded = [g for g in full if not t[g] > 0]
        self.order = [g for g in full if t[g] > 0]
        self.gidx = np.array([loc.genes[g] for g in self.order])
        self.t0 = {g: t[g] for g in self.order}
        self.excluded_flagged = [g for g in self.excluded if g in loc.flags]

    def totals(self, m):
        t = abc_mod.total_input(self.loc.activities, m, self.loc.elements, self.gidx)
        return dict(zip(self.order, t))

    def delta(self, a, b, q, sigma=sweep_mod.SIGMA_BINS):
        m = sweep_mod.apply_staple(self.loc.contact, a, b, q, sigma)
        return {g: (v - self.t0[g]) / self.t0[g] for g, v in self.totals(m).items()}

    def effect(self, a, b, target, sigma=sweep_mod.SIGMA_BINS, sweep_if_over=None, direction=1):
        """Headline effect, and the q-sweep unless the headline is below `sweep_if_over`.
        `direction` is +1 to raise the target's dose, -1 to lower it."""
        q, noop = sweep_mod.staple_strength(self.loc.contact, a, b)
        delta = self.delta(a, b, q, sigma)
        sweep = None
        if sweep_if_over is None or direction * delta[target] >= sweep_if_over:
            sweep = [
                self.delta(a, b, qi, sigma)[target]
                for qi in sweep_mod.q_grid(self.loc.contact, abs(a - b))
            ]
        return {
            "delta": delta,
            "q": q,
            "noop": noop,
            "sweep": sweep,
            "sign_stable": bool(sweep is not None and all(direction * s > 0 for s in sweep)),
        }


def nearest_permitted(open_mask, forbidden, target_bin):
    ok = anchors_mod.candidate_anchors(open_mask, forbidden, target_bin, ANCHOR_REACH)
    return min(ok, key=lambda b: abs(b - target_bin)) if ok else None


def anchor_space(loc, share, order):
    """(forbidden mask, open-chromatin mask, per-category counts) under the relaxed rules."""
    n = len(loc.contact)
    on = np.array([loc.expression.get(g, 0.0) >= EXPRESSED_TPM for g in order])
    functional = loc.elements[
        np.where((share[:, on] >= abc_mod.ABC_THRESHOLD).any(axis=1))[0]
    ]
    spans = loc.exons or list(loc.gene_bodies.values())
    forbidden = anchors_mod.forbidden_mask(
        n, functional, loc.genes.values(), spans, loc.ctcf
    )
    open_mask = anchors_mod.near_open_chromatin(loc.accessible)
    ok = open_mask & ~forbidden
    return (
        forbidden,
        open_mask,
        {
            "genes_expressed": int(on.sum()),
            "genes": len(order),
            "elements_protected": int(len(functional)),
            "elements": int(len(loc.elements)),
            "atac_peaks": int(loc.accessible.sum()),
            "bins_near_open_chromatin": int(open_mask.sum()),
            "forbidden_bins": int(forbidden.sum()),
            "permitted_anchors": int(ok.sum()),
            "bins": n,
        },
    )


def search(win, sc, forbidden, open_mask, targets=None):
    """Rank staple designs. With `targets`, rank designs for those named genes (the
    shortlisting use); without, sweep every gene (exploratory, noisy)."""
    loc = sc.loc
    cands, rows, dropped = [], [], {"noop": 0, "no_anchor": 0, "duplicate": 0}
    seen = set()
    n_enh = TARGET_ENHANCERS if targets else TOP_ENHANCERS
    for ei in np.argsort(-loc.activities)[:n_enh]:
        e_bin = int(loc.elements[ei])
        a = nearest_permitted(open_mask, forbidden, e_bin)
        for g in (targets or sc.order):
            b = nearest_permitted(open_mask, forbidden, loc.genes[g])
            if a is None or b is None or abs(a - b) < MIN_SEPARATION_BINS:
                dropped["no_anchor"] += 1
                continue
            if (a, b, g) in seen:                 # two enhancers resolving to one anchor pair is one design
                dropped["duplicate"] += 1
                continue
            seen.add((a, b, g))
            r = sc.effect(a, b, g, sweep_if_over=MIN_EFFECT)
            if r["noop"]:
                dropped["noop"] += 1
                continue
            risk = risk_score(r["delta"], loc.flags, TAU, target=g)
            cand = Candidate(
                name=f"{g} <- {win.kb(e_bin)}kb",
                delta_t=r["delta"][g],
                provenance="hic",
                distance_kb=abs(a - b) * BIN_KB,
                abstained=False,
                sign_stable=r["sign_stable"],
                safe=risk.tier != "BLOCKED",
            )
            unflagged = sorted(
                (
                    (x, d)
                    for x, d in r["delta"].items()
                    if x != g and x not in loc.flags and abs(d) > TAU
                ),
                key=lambda kv: -abs(kv[1]),
            )
            cands.append(cand)
            rows.append(
                {
                    "name": cand.name,
                    "mode": "up",
                    "target": g,
                    "target_tpm": round(loc.expression.get(g, 0.0), 2),
                    "target_flag": loc.flags.get(g),
                    "baseline_T": round(sc.t0[g], 2),
                    "enhancer_kb": win.kb(e_bin),
                    "anchor_a_kb": win.kb(a),
                    "anchor_b_kb": win.kb(b),
                    "separation_kb": abs(a - b) * BIN_KB,
                    "q": round(r["q"], 3),
                    "delta_target": round(r["delta"][g], 4),
                    "sweep": None
                    if r["sweep"] is None
                    else [round(x, 4) for x in r["sweep"]],
                    "delta_at_p99": None if r["sweep"] is None else round(r["sweep"][-1], 4),
                    "sign_stable": r["sign_stable"],
                    "confidence": cand.tier,
                    "risk": {
                        "tier": risk.tier,
                        "score": round(risk.score, 4),
                        "worst_gene": risk.worst_gene,
                        "n_flagged_in_window": risk.n_flagged_in_window,
                        "reason": risk.reason,
                    },
                    "unflagged_collateral_over_tau": [
                        {"gene": x, "delta": round(d, 4)} for x, d in unflagged[:6]
                    ],
                }
            )
    ranked = rank(cands)
    place = {c.name: i + 1 for i, c in enumerate(ranked)}
    for r in rows:
        r["rank"] = place.get(r["name"])
    dropped["blocked_by_risk"] = int(sum(not c.safe and c.delta_t >= MIN_EFFECT for c in cands))   # would have met the dose goal
    dropped["below_min_effect"] = int(
        sum(c.safe and c.delta_t < MIN_EFFECT for c in cands)
    )
    return ranked, rows, dropped


DECOY_REACH = 100          # decoy anchor at most +/-500 kb from the enhancer anchor (the confident-distance limit)
DECOY_KEEP_OUT = 12        # ... and at least 60 kb from the target's own promoter, or it would raise the target
N_DECOYS = 40              # evenly subsample the permitted decoys per enhancer
DOWN_ENHANCERS = 5         # the target's five strongest enhancers by ABC share


def search_down(win, sc, forbidden, open_mask, targets, share):
    """Rank staple designs that LOWER a target's dose (sequestration).

    The mechanism is the same one that lets a point edit move nothing: the contact
    matrix is conserved, so pulling an enhancer into a new contact elsewhere (a
    "decoy" anchor) takes contact away from its existing partners, the target among
    them. Everything downstream (risk score, gate, floor, sign stability) is the
    up-regulation machinery with the sign flipped. This rests on the row-sum
    conservation assumption, the least tested step in the model (docs/method.md).
    """
    loc = sc.loc
    n = len(loc.contact)
    cands, rows, dropped = [], [], {"noop": 0, "no_anchor": 0, "duplicate": 0}
    seen = set()
    permitted = np.where(open_mask & ~forbidden)[0]
    for g in targets:
        gi = sc.order.index(g)
        t_bin = loc.genes[g]
        top = np.argsort(-share[:, gi])[:DOWN_ENHANCERS]
        for ei in top:
            e_bin = int(loc.elements[ei])
            a = nearest_permitted(open_mask, forbidden, e_bin)
            if a is None:
                dropped["no_anchor"] += 1
                continue
            pool = [int(b) for b in permitted
                    if MIN_SEPARATION_BINS <= abs(b - a) <= DECOY_REACH and abs(b - t_bin) >= DECOY_KEEP_OUT]
            if len(pool) > N_DECOYS:
                pool = [pool[i] for i in np.linspace(0, len(pool) - 1, N_DECOYS).astype(int)]
            for b in pool:
                if (a, b, g) in seen:
                    dropped["duplicate"] += 1
                    continue
                seen.add((a, b, g))
                r = sc.effect(a, b, g, sweep_if_over=MIN_EFFECT, direction=-1)
                if r["noop"]:
                    dropped["noop"] += 1
                    continue
                risk = risk_score(r["delta"], loc.flags, TAU, target=g)
                cand = Candidate(
                    name=f"{g} v {win.kb(e_bin)}kb->{win.kb(b)}kb",
                    delta_t=-r["delta"][g],                      # dose REDUCTION, so rank() can be reused
                    provenance="hic",
                    distance_kb=abs(a - b) * BIN_KB,
                    abstained=False,
                    sign_stable=r["sign_stable"],
                    safe=risk.tier != "BLOCKED",
                )
                unflagged = sorted(
                    ((x, d) for x, d in r["delta"].items()
                     if x != g and x not in loc.flags and abs(d) > TAU),
                    key=lambda kv: -abs(kv[1]),
                )
                cands.append(cand)
                rows.append({
                    "name": cand.name, "mode": "down", "target": g,
                    "target_tpm": round(loc.expression.get(g, 0.0), 2),
                    "target_flag": loc.flags.get(g),
                    "baseline_T": round(sc.t0[g], 2),
                    "enhancer_kb": win.kb(e_bin),
                    "enhancer_abc_share": round(float(share[ei, gi]), 3),
                    "anchor_a_kb": win.kb(a), "anchor_b_kb": win.kb(b),
                    "decoy_kb": win.kb(b),
                    "separation_kb": abs(a - b) * BIN_KB,
                    "q": round(r["q"], 3),
                    "delta_target": round(r["delta"][g], 4),
                    "sweep": None if r["sweep"] is None else [round(x, 4) for x in r["sweep"]],
                    "delta_at_p99": None if r["sweep"] is None else round(r["sweep"][-1], 4),
                    "sign_stable": r["sign_stable"],
                    "confidence": cand.tier,
                    "risk": {"tier": risk.tier, "score": round(risk.score, 4), "worst_gene": risk.worst_gene,
                             "n_flagged_in_window": risk.n_flagged_in_window, "reason": risk.reason},
                    "unflagged_collateral_over_tau": [{"gene": x, "delta": round(d, 4)} for x, d in unflagged[:6]],
                })
    ranked = rank(cands)
    place = {c.name: i + 1 for i, c in enumerate(ranked)}
    for r in rows:
        r["rank"] = place.get(r["name"])
    dropped["blocked_by_risk"] = int(sum(not c.safe and c.delta_t >= MIN_EFFECT for c in cands))   # would have met the dose goal
    dropped["below_min_effect"] = int(sum(c.safe and c.delta_t < MIN_EFFECT for c in cands))
    return ranked, rows, dropped


def show(win, loc, ranked, rows, dropped, space, limit=10, direction="up"):
    print(
        f"\n{'=' * 100}\n{win.label}  ·  {win.chrom}:{win.start:,}-{win.end:,}  ·  K562\n{'=' * 100}"
    )
    print(
        f"{space['genes_expressed']}/{space['genes']} genes expressed · {space['elements_protected']}/{space['elements']} "
        f"elements protected · {space['atac_peaks']} ATAC peaks · {space['permitted_anchors']} permitted anchors"
    )
    flagged = (
        ", ".join(
            f"{g} ({f})"
            for g, f in sorted(loc.flags.items(), key=lambda kv: loc.genes[kv[0]])
        )
        or "none"
    )
    print(f"ClinGen dosage-sensitive genes in window ({len(loc.flags)}): {flagged}")
    print(
        f"{len(rows)} staples scored · {len(ranked)} reach dose {'rise' if direction == 'up' else 'fall'} >= {MIN_EFFECT:.0%} · dropped {dropped}\n"
    )
    sign = 1 if direction == "up" else -1
    blocked = sorted((r for r in rows if r["risk"]["tier"] == "BLOCKED" and sign * r["delta_target"] >= MIN_EFFECT),
                     key=lambda r: -r["risk"]["score"])
    if blocked:
        print(f"  REJECTED BY THE SAFETY GATE ({len(blocked)}):")
        for r in blocked[:4]:
            print(f"    x {r['target']:<8}{r['separation_kb']:>5} kb  target dT/T {r['delta_target']:>+7.1%}   {r['risk']['reason']}")
        print()
    if not ranked:
        print("  no staple reaches the dose goal in this window")
        pick = max if direction == "up" else min
        best = pick(rows, key=lambda r: r["delta_target"], default=None)
        if best:
            print(f"  best available: {best['name']} {best['delta_target']:+.1%}")
        return
    print(
        f"  {'#':<3}{'target':<9}{'TPM':>7}{'sep kb':>8}{'dT/T':>9}   {'risk':<9}{'score':>7}  {'confidence':<11}worst flagged gene / collateral"
    )
    by = {r["name"]: r for r in rows}
    for c in ranked[:limit]:
        r = by[c.name]
        k = r["risk"]
        note = k["reason"] if k["tier"] != "CLEAR" else "no flagged gene in window"
        extra = ", ".join(
            f"{x['gene']} {x['delta']:+.0%}"
            for x in r["unflagged_collateral_over_tau"][:2]
        )
        print(
            f"  {r['rank']:<3}{r['target']:<9}{r['target_tpm']:>7.0f}{r['separation_kb']:>8}{r['delta_target']:>+9.1%}   "
            f"{k['tier']:<9}{k['score']:>7.1%}  {r['confidence']:<11}{note}"
            + (f"  | also moves {extra}" if extra else "")
        )


def run_window(win: Window, name: str, targets=None, write: bool = True, direction: str = "up",
               verbose: bool = True) -> dict:
    loc, prov = load_k562_locus(win.chrom, win.start, win.end)
    sc = Scorer(loc)
    share = abc_mod.abc_share(loc.activities, loc.contact, loc.elements, sc.gidx)
    forbidden, open_mask, space = anchor_space(loc, share, sc.order)
    space["excluded_genes"], space["excluded_flagged"] = sc.excluded, sc.excluded_flagged
    if direction == "down":
        ranked, rows, dropped = search_down(win, sc, forbidden, open_mask, targets, share)
    else:
        ranked, rows, dropped = search(win, sc, forbidden, open_mask, targets)
    if verbose:
        show(win, loc, ranked, rows, dropped, space, direction=direction)
    out = {
        "label": win.label,
        "window": prov["window"],
        "targets": targets,
        "direction": direction,
        "provenance": prov,
        "config": CONFIG,
        "tests_evaluated": TESTS_EVALUATED,
        "anchor_space": space,
        "flagged_genes_in_window": {g: f for g, f in loc.flags.items()},
        "n_staples_scored": len(rows),
        "n_meeting_dose_goal": len(ranked),
        "dropped": dropped,
        "ranked": [
            r for r in sorted((r for r in rows if r["rank"]), key=lambda r: r["rank"])
        ],
        "best_below_floor": sorted(
            (r for r in rows if not r["rank"]),
            key=lambda r: -r["delta_target"] if direction == "up" else r["delta_target"],   # closest to the goal first
        )[:5],
    }
    if write:
        DERIVED.mkdir(parents=True, exist_ok=True)
        (DERIVED / f"{name}.json").write_text(json.dumps(out, indent=1))
    return out

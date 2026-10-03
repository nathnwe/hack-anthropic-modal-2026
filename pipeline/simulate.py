"""Run the whole method end to end on one locus and return the numbers.

The document and every figure read from this function, so a figure cannot
disagree with the prose: there is one source of arithmetic.

    uv run --with numpy python -m pipeline.simulate
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pipeline import abc as abc_mod
from pipeline import anchors as anchors_mod
from pipeline import sweep as sweep_mod
from pipeline.locus import BIN_KB, DEMO_ANCHOR_A, DEMO_ANCHOR_B, Locus, demo_locus
from scoring.rank import Candidate
from scoring.risk import TAU, concerns

SWEEP_PERCENTILES = (50, 75, 90, 95, 99)


@dataclass
class Result:
    locus: Locus
    before: np.ndarray
    after: np.ndarray
    q: float
    t_before: dict[str, float]
    t_after: dict[str, float]
    delta: dict[str, float]           # dT/T per gene
    share_before: np.ndarray
    share_after: np.ndarray
    forbidden: np.ndarray
    sweep: dict[str, np.ndarray]      # q grid and dT/T per gene across it
    concerns: list
    candidate: Candidate


def _totals(loc: Locus, contact: np.ndarray) -> dict[str, float]:
    order = loc.gene_order
    idx = np.array([loc.genes[g] for g in order])
    t = abc_mod.total_input(loc.activities, contact, loc.elements, idx)
    return dict(zip(order, t))


def run(loc: Locus | None = None, anchor_a: int = DEMO_ANCHOR_A,
        anchor_b: int = DEMO_ANCHOR_B, target: str = "TARGET") -> Result:
    loc = loc or demo_locus()
    before = loc.contact
    sep = abs(anchor_b - anchor_a)
    q = sweep_mod.yardstick(before, sep)
    after = sweep_mod.apply_staple(before, anchor_a, anchor_b, q)

    t0, t1 = _totals(loc, before), _totals(loc, after)
    delta = {g: (t1[g] - t0[g]) / t0[g] for g in t0}

    order = loc.gene_order
    gidx = np.array([loc.genes[g] for g in order])
    share0 = abc_mod.abc_share(loc.activities, before, loc.elements, gidx)
    share1 = abc_mod.abc_share(loc.activities, after, loc.elements, gidx)

    forbidden = anchors_mod.forbidden_mask(
        len(before), loc.elements, loc.genes.values(),
        loc.gene_bodies.values(), loc.ctcf)

    grid = sweep_mod.q_grid(before, sep, SWEEP_PERCENTILES)
    curves = {g: [] for g in order}
    ranks = []
    for qi in grid:
        ti = _totals(loc, sweep_mod.apply_staple(before, anchor_a, anchor_b, qi))
        d = {g: (ti[g] - t0[g]) / t0[g] for g in order}
        for g in order:
            curves[g].append(d[g])
        ranks.append(max(order, key=lambda g: d[g]))
    curves = {g: np.array(v) for g, v in curves.items()}

    found = concerns(delta, loc.flags, TAU, target=target)
    cand = Candidate(
        name=f"{target} <- anchor@{anchor_a}",
        delta_t=delta[target],
        provenance="hic",
        distance_kb=sep * BIN_KB,
        abstained=False,
        rank_stable=len(set(ranks)) == 1,
        safe=not found,
    )
    return Result(loc, before, after, q, t0, t1, delta, share0, share1,
                  forbidden, {"q": grid, **curves}, found, cand)


if __name__ == "__main__":
    r = run()
    print(f"q = {r.q:.4f}   separation = {r.candidate.distance_kb:.0f} kb")
    print(f"anchors permitted: A={not r.forbidden[DEMO_ANCHOR_A]} "
          f"B={not r.forbidden[DEMO_ANCHOR_B]}")
    for g in r.locus.gene_order:
        flag = r.locus.flags.get(g, "")
        print(f"  {g:<7} {r.delta[g]:+7.1%}  {flag}")
    print("concerns:", [str(c) for c in r.concerns] or "none")
    print("tier:", r.candidate.tier, " safe:", r.candidate.safe,
          " stable:", r.candidate.rank_stable)

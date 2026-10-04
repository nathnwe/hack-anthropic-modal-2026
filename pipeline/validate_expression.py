"""Does baseline T track measured K562 expression across the window?

    uv run --with hic-straw --with pybigtools --with numpy python -m pipeline.validate_expression

Spearman rank correlation of T_G against log(1 + TPM) over every coding gene
in the window, against a null that keeps the same matrix, elements and
expression but scatters the promoter positions at random. Writes
data/derived/abc_vs_expression.json.

Read the result for what it is: ABC knows nothing about expression, so this
measures whether genes sit in active regions, not whether a given gene is on
(an unexpressed olfactory receptor, OR51V1, scores T = 746).
"""

from __future__ import annotations

import json

import numpy as np

from pipeline import abc as abc_mod
from pipeline.real_data import load_k562_locus
from pipeline.real_run import CHROM, DERIVED, END, START


def ranks(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x))
    r[order] = np.arange(len(x))
    for v in np.unique(x):                      # average ties
        m = x == v
        if m.sum() > 1:
            r[m] = r[m].mean()
    return r


def spearman(a, b) -> float:
    return float(np.corrcoef(ranks(np.asarray(a)), ranks(np.asarray(b)))[0, 1])


def main(n_null: int = 2000, seed: int = 0) -> dict:
    loc, _ = load_k562_locus(CHROM, START, END)
    order = loc.gene_order
    gidx = np.array([loc.genes[g] for g in order])
    t = abc_mod.total_input(loc.activities, loc.contact, loc.elements, gidx)
    y = np.log1p([loc.expression[g] for g in order])
    rho = spearman(t, y)

    rng = np.random.default_rng(seed)
    n = len(loc.contact)
    null = np.array([spearman(abc_mod.total_input(loc.activities, loc.contact, loc.elements,
                                                  rng.choice(n, len(order), replace=False)), y)
                     for _ in range(n_null)])
    out = {"n_genes": len(order), "spearman_T_vs_log1p_tpm": round(rho, 3),
           "null_mean": round(float(null.mean()), 3), "null_95th_pct": round(float(np.percentile(null, 95)), 3),
           "fraction_of_null_below_observed": float((null < rho).mean()), "n_null": n_null,
           "top_T_genes": [[order[i], round(float(t[i])), round(float(np.expm1(y[i])))] for i in np.argsort(-t)[:6]]}
    DERIVED.mkdir(parents=True, exist_ok=True)
    (DERIVED / "abc_vs_expression.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    main()

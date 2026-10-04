"""The intervention: what a staple does to a contact matrix.

Two design choices here carry the whole safety argument.

1. The edit is a KERNEL, not a point.

   ABC normalises per gene. If you set a single cell C[A, B] you can only ever
   move gene B -- every other gene's denominator is untouched. Enhancer
   hijacking becomes arithmetically impossible and the off-target analysis
   returns zeros by construction. A model that cannot express the failure mode
   cannot rule it out. So the staple deposits a Gaussian blob:

       C'[i, j] = max( C[i, j], q * K(i - A, j - B) )

   sigma ~= 2-3 bins (10-15 kb). Empirical support: ENCODE-rE2G found
   significant CRISPRi crosstalk between elements 1-10 kb apart, P = 1.8e-5.

2. ROW SUMS ARE PRESERVED.

   A locus has a finite contact budget. Without this, a staple adds closeness
   for free, nothing can go down, and enhancer theft is unrepresentable.

Both are assumptions, not measurements. They are listed as such in docs/method.md.
"""

from __future__ import annotations

import numpy as np

SIGMA_BINS = 3.0  # ~15 kb at 5 kb resolution
RENORM_ITERS = 100   # 4 left up to 0.9% row-sum error; 100 converges to ~1e-16


def kernel(n: int, anchor_a: int, anchor_b: int, sigma: float = SIGMA_BINS) -> np.ndarray:
    """Symmetric Gaussian blob centred on the stapled pair."""
    i, j = np.indices((n, n))
    k = np.exp(-(((i - anchor_a) ** 2 + (j - anchor_b) ** 2) / (2 * sigma ** 2)))
    return np.maximum(k, k.T)


def apply_staple(contact: np.ndarray, anchor_a: int, anchor_b: int, q: float,
                 sigma: float = SIGMA_BINS, conserve: bool = True) -> np.ndarray:
    """C' = max(C, q*K), then restore every row's original contact budget."""
    n = contact.shape[0]
    out = np.maximum(contact, q * kernel(n, anchor_a, anchor_b, sigma))
    if not conserve:
        return out
    budget = contact.sum(axis=1)
    for _ in range(RENORM_ITERS):
        scale = np.divide(budget, out.sum(axis=1), out=np.ones(n), where=out.sum(axis=1) > 0)
        out = out * scale[:, None]
        out = (out + out.T) / 2.0
    return out


def staple_strength(contact: np.ndarray, a: int, b: int, tol: float = 0.95) -> tuple[float, bool]:
    """(q, noop). The yardstick, plus a flag for a staple that cannot help.

    The yardstick is the strongest real contact at this separation ANYWHERE in
    the window. If the pair being stapled already is that contact, then
    max(C, q*K) leaves the anchor pixel unchanged and renormalisation then
    slightly weakens it (verified: LCR->HBG2, 49.13 -> 47.62). Such a staple is
    reported as already-at-yardstick, not scored as if it did something.
    """
    q = yardstick(contact, abs(a - b))
    return q, bool(contact[a, b] >= tol * q)


def yardstick(contact: np.ndarray, separation_bins: int) -> float:
    """A describable q instead of a free parameter.

    The strongest real contact anywhere in this window at the given genomic
    separation. 'As tight as the tightest natural loop at 100 kb' is a claim a
    reviewer can check; '0.4' is not.
    """
    return float(np.diagonal(contact, offset=separation_bins).max())


def q_grid(contact: np.ndarray, separation_bins: int,
           percentiles=(50, 75, 90, 95, 99)) -> np.ndarray:
    """Sweep the staple strength across real contacts at this separation.

    Nobody has measured the contact frequency a dCas staple produces, so we
    refuse to pick one and report the dose-response instead.
    """
    band = np.diagonal(contact, offset=separation_bins)
    return np.percentile(band, percentiles)

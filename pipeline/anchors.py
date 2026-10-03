"""Where a staple is allowed to land.

An anchor is a stretch of accessible chromatin the dCas machinery can bind
without destroying the thing we are trying to use. We subtract, from the ATAC
peaks, everything whose function we would disrupt by parking a protein on it.
This is a filter, not a score: it removes the impossible rather than ranking
the possible.
"""

from __future__ import annotations

import numpy as np

TSS_PAD_BINS = 1  # +/- 5 kb around a TSS


def forbidden_mask(n: int, elements, gene_tss, gene_bodies, ctcf_sites,
                   tss_pad: int = TSS_PAD_BINS) -> np.ndarray:
    """True where a staple must not be placed."""
    mask = np.zeros(n, dtype=bool)
    for e in elements:
        mask[e] = True                       # the enhancer itself
    for t in gene_tss:
        mask[max(0, t - tss_pad):t + tss_pad + 1] = True   # promoter
    for start, end in gene_bodies:
        mask[start:end + 1] = True           # transcribed region
    for c in ctcf_sites:
        mask[c] = True                       # boundary element
    return mask


def candidate_anchors(accessible: np.ndarray, forbidden: np.ndarray,
                      near: int, window: int = 8) -> list[int]:
    """Accessible, permitted bins within `window` bins of a point of interest."""
    lo, hi = max(0, near - window), min(len(accessible), near + window + 1)
    return [b for b in range(lo, hi) if accessible[b] and not forbidden[b]]

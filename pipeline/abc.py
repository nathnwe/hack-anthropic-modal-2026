"""Activity-by-Contact: the baseline arithmetic.

Fulco et al. 2019, Nat Genet 51:1664-1669.

Two numbers come out of this module and they are NOT the same thing:

    share[e, g]  the ABC score -- what fraction of gene g's regulatory input
                 element e is responsible for. Normalised per gene, so it
                 sums to 1 down each column. This is ATTRIBUTION.

    total[g]     T_G = sum_e A_e * C(e, g). Not normalised. This is DOSAGE,
                 and it is the quantity a staple is supposed to move.

Reading the share after an intervention tells you about redistribution, not
about increase. Report dT/T; use the share to say which element did it.
"""

from __future__ import annotations

import numpy as np

ABC_THRESHOLD = 0.02  # Fulco 2019, calibrated against CRISPRi-FlowFISH


def activity(atac: np.ndarray, h3k27ac: np.ndarray | None = None) -> np.ndarray:
    """Element activity A_e.

    The published definition is the geometric mean of DNase and H3K27ac read
    counts. We default to accessibility alone, which is what ENCODE-rE2G's own
    core model uses, because H3K27ac is unavailable for most cell types. Pass
    h3k27ac to recover the published form.
    """
    atac = np.asarray(atac, dtype=float)
    if h3k27ac is None:
        return atac
    return np.sqrt(atac * np.asarray(h3k27ac, dtype=float))


def total_input(activities: np.ndarray, contact: np.ndarray,
                elements: np.ndarray, genes: np.ndarray) -> np.ndarray:
    """T_G for each gene: the unnormalised sum of activity-weighted contact."""
    return activities @ contact[np.ix_(elements, genes)]


def abc_share(activities: np.ndarray, contact: np.ndarray,
              elements: np.ndarray, genes: np.ndarray) -> np.ndarray:
    """ABC score for every (element, gene) pair. Columns sum to 1."""
    weighted = activities[:, None] * contact[np.ix_(elements, genes)]
    denom = weighted.sum(axis=0)
    return np.divide(weighted, denom, out=np.zeros_like(weighted), where=denom > 0)

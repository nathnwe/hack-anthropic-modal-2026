"""A synthetic locus with realistic statistics.

ILLUSTRATIVE. Element positions, activities and gene positions are invented.
The contact matrix is generated from a distance-decay power law with TAD
blocks and a few CTCF loops, which is the shape real Hi-C has, but it is not
real Hi-C. Every figure built on this is labelled accordingly.

The point of the demo is that the ARITHMETIC downstream of here is the real
pipeline arithmetic. Only the inputs are made up.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

BIN_KB = 5
N_BINS = 400          # 2 Mb window
DECAY = 1.0           # contact ~ distance^-1, the usual polymer scaling


@dataclass(frozen=True)
class Locus:
    contact: np.ndarray
    elements: np.ndarray          # bin index of each candidate enhancer
    activities: np.ndarray        # A_e
    genes: dict[str, int]         # symbol -> TSS bin
    gene_bodies: dict[str, tuple[int, int]]
    ctcf: list[int]
    tad_edges: list[int]
    accessible: np.ndarray        # bool per bin
    flags: dict[str, str] = field(default_factory=dict)
    illustrative: bool = True

    @property
    def gene_order(self) -> list[str]:
        return sorted(self.genes, key=self.genes.get)

    def kb(self, b: int) -> float:
        return b * BIN_KB


def _contact_matrix(n: int, tad_edges: list[int], loops, rng) -> np.ndarray:
    i, j = np.indices((n, n))
    d = np.abs(i - j)
    c = 1.0 / np.power(d + 1.0, DECAY)
    tad = np.zeros(n, dtype=int)
    for k in range(len(tad_edges) - 1):
        tad[tad_edges[k]:tad_edges[k + 1]] = k
    c *= np.where(tad[:, None] == tad[None, :], 1.0, 0.30)   # boundary insulation
    for a, b, strength in loops:                              # CTCF corner peaks
        blob = np.exp(-(((i - a) ** 2 + (j - b) ** 2) / (2 * 2.5 ** 2)))
        c += strength * np.maximum(blob, blob.T)
    c *= rng.lognormal(0.0, 0.18, size=(n, n))                # experimental noise
    c = (c + c.T) / 2.0
    np.fill_diagonal(c, 0.0)
    return c


def demo_locus(seed: int = 7) -> Locus:
    """One 2 Mb window: a target gene, five bystanders, six enhancers.

    Laid out so the demo exercises every branch -- an intended gain, a
    collateral gain on a near neighbour, a loss on a flagged gene whose
    enhancer we steal contact from, and distal genes that must not move.
    """
    rng = np.random.default_rng(seed)
    tad_edges = [0, 150, 280, N_BINS]
    ctcf = [150, 280, 90, 340]
    loops = [(150, 280, 0.22), (90, 150, 0.16), (280, 340, 0.14)]
    contact = _contact_matrix(N_BINS, tad_edges, loops, rng)

    elements = np.array([60, 108, 192, 232, 262, 350])
    activities = np.array([1.5, 1.1, 1.3, 2.8, 1.6, 1.2])

    genes = {"FAR": 40, "MIDDLE": 120, "ONCO": 176, "PARTNER": 180,
             "NEIGHB": 302, "TARGET": 308, "DISTAL": 372}
    gene_bodies = {"FAR": (40, 46), "MIDDLE": (120, 128), "ONCO": (170, 176),
                   "PARTNER": (180, 186), "NEIGHB": (296, 302),
                   "TARGET": (308, 316), "DISTAL": (372, 378)}
    flags = {"ONCO": "oncogene", "MIDDLE": "haploinsufficient"}

    accessible = np.zeros(N_BINS, dtype=bool)
    accessible[rng.choice(N_BINS, size=110, replace=False)] = True
    accessible[elements] = True
    accessible[[234, 306, 178, 110]] = True   # the anchor pairs the demo staples
    for t in genes.values():
        accessible[t] = True
    for c in ctcf:
        accessible[c] = True

    return Locus(contact=contact, elements=elements, activities=activities,
                 genes=genes, gene_bodies=gene_bodies, ctcf=ctcf,
                 tad_edges=tad_edges, accessible=accessible, flags=flags)


#: The staple the demo applies: a strong enhancer in TAD 1 brought to TARGET.
DEMO_ANCHOR_A, DEMO_ANCHOR_B = 234, 306
DEMO_ELEMENT, DEMO_TARGET = 232, "TARGET"

#: A second staple that fails the safety gate: the intended gene PARTNER sits
#: 40 kb from an oncogene, and the kernel cannot tell them apart.
UNSAFE_ANCHOR_A, UNSAFE_ANCHOR_B = 110, 178
UNSAFE_TARGET = "PARTNER"

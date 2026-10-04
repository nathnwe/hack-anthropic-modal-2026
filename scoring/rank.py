"""Confidence, then the ordering.

Confidence is not a vibe. It is one question -- does the answer survive the
things we do not know? -- asked four ways, all deterministic:

  abstained       UniversalEPI's K=10 deep ensemble returns a maximum-confidence
                  fold change of exactly 0 when the prediction intervals
                  overlap. That is a free, principled abstain.
  provenance      measured Hi-C 53.3% precision > averaged megamap 48.0% >
                  inverse-distance 44.3% (ENCODE-rE2G benchmark).
  distance        every model tested failed three real MYC enhancers at ~2 Mb.
  sign stability  is the predicted change in the same direction at every q from
                  p50 to p99 of real contacts? (Full gene-ordering stability is
                  too strict: it was never satisfied by any of 592 real staples.)

Output order: SAFE is a gate; among survivors, sort by dT/T. Confidence rides
along as a label rather than a sort key, so a confident small effect never
outranks an uncertain large one by accident -- the researcher sees both.
"""

from __future__ import annotations

from dataclasses import dataclass

MAX_CONFIDENT_DISTANCE_KB = 500
MIN_EFFECT = 0.10   # assumption: a staple must move the target's dose by at least this much to be a candidate
TRUSTED_PROVENANCE = {"hic", "universalepi"}
PROVENANCE_PRECISION = {"hic": 0.533, "megamap": 0.480, "power_law": 0.443}


@dataclass(frozen=True)
class Candidate:
    name: str
    delta_t: float            # dT/T on the target gene
    provenance: str           # hic | megamap | power_law
    distance_kb: float
    abstained: bool
    sign_stable: bool
    safe: bool

    @property
    def confident(self) -> bool:
        return (not self.abstained
                and self.provenance in TRUSTED_PROVENANCE
                and self.distance_kb <= MAX_CONFIDENT_DISTANCE_KB
                and self.sign_stable)

    @property
    def tier(self) -> str:
        if self.abstained:
            return "ABSTAIN"
        return "CONFIDENT" if self.confident else "TENTATIVE"


def rank(candidates: list[Candidate], min_effect: float = MIN_EFFECT) -> list[Candidate]:
    """Gate on safety, drop anything that does not raise the target's dose by
    at least min_effect, then order by predicted dose change. Deterministic:
    ties on the rounded value break by name."""
    survivors = [c for c in candidates
                 if c.safe and not c.abstained and c.delta_t >= min_effect]
    return sorted(survivors, key=lambda c: (-round(c.delta_t, 4), c.name))

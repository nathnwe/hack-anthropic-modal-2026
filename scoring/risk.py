"""Which genes must not move, and the gate that enforces it.

Three fixed lists, no model in the loop:

  oncogene / tumour suppressor   COSMIC Cancer Gene Census
  dosage-sensitive               ClinGen haploinsufficiency + triplosensitivity
  essential                      DepMap common-essential

A candidate fails SAFE if any flagged gene in the window changes its total
regulatory input by more than tau. tau is an assumption; we start at 10%.

Asymmetric trust, per CLAUDE.md: a language model may RAISE a concern with a
citation. It may never clear one. Silence from a model is not evidence of
safety.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

TAU = 0.10  # relative change in T that counts as "moved"

FLAG_LABELS = {
    "oncogene": "COSMIC Cancer Gene Census",
    "tumour_suppressor": "COSMIC Cancer Gene Census",
    "triplosensitive": "ClinGen dosage sensitivity",
    "haploinsufficient": "ClinGen dosage sensitivity",
    "essential": "DepMap common-essential",
}


@dataclass(frozen=True)
class Concern:
    gene: str
    flag: str
    delta: float
    source: str

    def __str__(self) -> str:
        return f"{self.gene} ({self.flag}) {self.delta:+.1%} -- {self.source}"


def concerns(deltas: dict[str, float], flags: dict[str, str],
             tau: float = TAU, target: str | None = None) -> list[Concern]:
    """Every flagged gene that moved too much, worst first.

    Two independent tests, because they fail for different reasons:

      absolute   the flagged gene moved more than tau
      relative   the flagged gene moved at least as much as the gene we were
                 aiming at. However small the numbers, a staple that does more
                 to an oncogene than to its target is not a candidate.

    Fail closed. A flagged gene whose change is undefined (NaN, e.g. from a
    zero baseline) is reported as a concern rather than waved through, and a
    target missing from `deltas` is an error rather than a silently disabled
    relative test.
    """
    if target is not None and target not in deltas:
        raise ValueError(f"target {target!r} is not in deltas; the relative test would be skipped")
    intended = abs(deltas[target]) if target is not None else None
    if intended is not None and not math.isfinite(intended):
        raise ValueError(f"target {target!r} has an undefined change ({deltas[target]})")
    found = []
    for g, d in deltas.items():
        if g not in flags or g == target:
            continue
        undefined = not math.isfinite(d)
        over_tau = abs(d) > tau
        # a null staple (intended == 0) must not make every untouched flagged gene "dominate"
        dominates = intended is not None and abs(d) > 0 and abs(d) >= intended
        if undefined or over_tau or dominates:
            found.append(Concern(g, flags[g], d, FLAG_LABELS.get(flags[g], "unsourced")))
    return sorted(found, key=lambda c: (not math.isfinite(c.delta), abs(c.delta)), reverse=True)


def is_safe(deltas: dict[str, float], flags: dict[str, str], tau: float = TAU,
            target: str | None = None) -> bool:
    """A hard gate, not a term in a sum.

    A weighted score would let a large predicted effect outrank a safety
    concern. A gate cannot. Pass `target` so the relative test runs too.
    """
    return not concerns(deltas, flags, tau, target)

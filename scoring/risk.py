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
    """
    intended = abs(deltas[target]) if target in deltas else None
    found = []
    for g, d in deltas.items():
        if g not in flags or g == target:
            continue
        over_tau = abs(d) > tau
        dominates = intended is not None and abs(d) >= intended
        if over_tau or dominates:
            found.append(Concern(g, flags[g], d, FLAG_LABELS.get(flags[g], "unsourced")))
    return sorted(found, key=lambda c: -abs(c.delta))


def is_safe(deltas: dict[str, float], flags: dict[str, str], tau: float = TAU) -> bool:
    """A hard gate, not a term in a sum.

    A weighted score would let a large predicted effect outrank a safety
    concern. A gate cannot.
    """
    return not concerns(deltas, flags, tau)

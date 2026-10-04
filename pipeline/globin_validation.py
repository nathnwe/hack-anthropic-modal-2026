"""Checks against known beta-globin biology. Not part of the shortlist.

  sanity       baseline ABC puts the LCR on the globin genes (it wins on activity,
               not contact, so this is a consistency check);
  validation   LCR -> each globin promoter with every exclusion waived, the way
               forced-looping experiments did it (Deng et al. 2012 Cell; 2014
               Nature; neither in K562), so a sign-and-competition check only;
  sensitivity  the same across sigma, because the collateral readout depends on
               it far more than the target does.

    SDKROOT=<older macOS SDK> uv run --with hic-straw --with pybigtools --with numpy \
        python -m pipeline.globin_validation
"""

from __future__ import annotations

import json

import numpy as np

from pipeline import abc as abc_mod
from pipeline.locus import BIN_KB
from pipeline.real_data import load_k562_locus
from pipeline.real_run import DERIVED, SIGMAS, Scorer, Window

WIN = Window("beta-globin", "chr11", 4_700_000, 6_700_000)
LCR_KB = (5_270, 5_310)
GLOBIN = ["HBE1", "HBG2", "HBG1", "HBD", "HBB"]


def main() -> dict:
    loc, _ = load_k562_locus(WIN.chrom, WIN.start, WIN.end)
    sc = Scorer(loc)
    share = abc_mod.abc_share(loc.activities, loc.contact, loc.elements, sc.gidx)
    gpos = {g: i for i, g in enumerate(sc.order)}
    out = {"sanity": {}, "per_gene": {}, "sigma_sensitivity": {}}

    print("SANITY — top ABC element per globin gene (the LCR wins on activity, not contact)")
    for g in GLOBIN:
        col = share[:, gpos[g]]
        top = int(np.argmax(col))
        e_kb = WIN.kb(loc.elements[top])
        out["sanity"][g] = {"top_element_kb": e_kb, "share": round(float(col[top]), 3),
                            "margin_over_second": round(float(col[top] / np.sort(col)[-2]), 2),
                            "is_lcr": LCR_KB[0] <= e_kb <= LCR_KB[1]}
        print(f"  {g:<5} {e_kb} kb share {col[top]:.2f} margin {out['sanity'][g]['margin_over_second']}x")

    kbs = loc.elements * BIN_KB + WIN.start_kb
    lcr = int(loc.elements[np.argmax(np.where((kbs >= LCR_KB[0]) & (kbs <= LCR_KB[1]), loc.activities, -1))])
    print(f"\nVALIDATION — LCR ({WIN.kb(lcr)} kb) -> each globin promoter, exclusions waived")
    print(f"  {'gene':<6}{'sep':>5}{'C pct':>7}{'TPM':>7}{'dT/T':>9}  note")
    for g in GLOBIN:
        b = loc.genes[g]
        r = sc.effect(lcr, b, g)
        pct = float((np.diagonal(loc.contact, offset=abs(lcr - b)) < loc.contact[lcr, b]).mean() * 100)
        coll = sorted(((x, d) for x, d in r["delta"].items() if x != g and abs(d) > 0.05), key=lambda kv: -abs(kv[1]))
        note = "NO-OP: already the strongest contact at this separation" if r["noop"] else \
            ", ".join(f"{x} {d:+.0%}" for x, d in coll[:3])
        print(f"  {g:<6}{abs(lcr - b) * BIN_KB:>5}{pct:>6.0f}%{loc.expression.get(g, 0):>7.0f}{r['delta'][g]:>+9.1%}  {note}")
        out["per_gene"][g] = {"sep_kb": abs(lcr - b) * BIN_KB, "contact_percentile": round(pct, 1),
                              "tpm": round(loc.expression.get(g, 0.0), 2), "q": round(r["q"], 3), "noop": r["noop"],
                              "delta": round(r["delta"][g], 4), "sweep": [round(x, 4) for x in r["sweep"]],
                              "collateral": [{"gene": x, "delta": round(d, 4)} for x, d in coll[:6]]}

    print("\nSIGMA SENSITIVITY — LCR -> HBB")
    for s in SIGMAS:
        r = sc.effect(lcr, loc.genes["HBB"], "HBB", sigma=s)
        others = sorted(((x, d) for x, d in r["delta"].items() if x != "HBB"), key=lambda kv: -abs(kv[1]))
        n2 = int(sum(abs(d) > 0.02 for _, d in others))
        print(f"  sigma {s}: HBB {r['delta']['HBB']:+.1%}  largest collateral {others[0][0]} {others[0][1]:+.1%}  genes >2%: {n2}")
        out["sigma_sensitivity"][str(s)] = {"hbb": round(r["delta"]["HBB"], 4), "largest_collateral": others[0][0],
                                            "largest_collateral_delta": round(others[0][1], 4), "genes_over_2pct": n2}
    DERIVED.mkdir(parents=True, exist_ok=True)
    (DERIVED / "lcr_validation.json").write_text(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    main()

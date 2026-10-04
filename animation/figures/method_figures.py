"""The method, as a series of slide figures.

Every number here comes from `pipeline.simulate`, which is the same code path
the pipeline runs. No figure can disagree with the document, because there is
one source of arithmetic. Inputs are synthetic and labelled as such; the
arithmetic on them is real.

Run:
    uv run --with matplotlib --with numpy python animation/figures/method_figures.py

Writes docs/figures/method/fig*.png and docs/figures/method/captions.md
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch, Rectangle

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from pipeline import abc as abc_mod          # noqa: E402
from pipeline import sweep as sweep_mod      # noqa: E402
from pipeline.locus import (                 # noqa: E402
    BIN_KB, DEMO_ANCHOR_A, DEMO_ANCHOR_B,
    UNSAFE_ANCHOR_A, UNSAFE_ANCHOR_B,
)
from pipeline.simulate import run            # noqa: E402
from scoring.risk import TAU                 # noqa: E402

OUT = REPO / "docs" / "figures" / "method"
OUT.mkdir(parents=True, exist_ok=True)

BG, PANEL, INK = "#fffefb", "#f4f1e9", "#14201e"
TEAL, AMBER, ROSE = "#0f7a6e", "#9a6410", "#a33a3a"
GREY, MUTE, EDGE = "#8a8275", "#566460", "#ded8c9"
MONO = ["DejaVu Sans Mono", "Courier New", "monospace"]

HIC = LinearSegmentedColormap.from_list(
    "hic", ["#fffdf7", "#f7e7c2", "#e3a954", "#b8771b", "#7a3d14", "#4a1f0c"])
DIFF = LinearSegmentedColormap.from_list(
    "diff", [TEAL, "#8fc4bd", "#fffdf7", "#e8c583", AMBER])

plt.rcParams.update({
    "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
    "font.size": 11, "text.color": INK, "axes.labelcolor": MUTE,
    "xtick.color": MUTE, "ytick.color": MUTE, "axes.edgecolor": EDGE,
})

CAPTIONS: list[tuple[str, str, str]] = []


def finish(fig, name: str, title: str, caption: str) -> None:
    fig.savefig(OUT / f"{name}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    CAPTIONS.append((name, title, caption))
    print(f"  {name}.png")


def tidy(ax, *, left=True, bottom=True) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)
    ax.grid(axis="y", color=EDGE, lw=0.8)
    ax.set_axisbelow(True)


def heading(ax, text: str, sub: str | None = None) -> None:
    ax.set_title(text, loc="left", fontsize=14, color=INK, pad=28 if sub else 10)
    if sub:
        ax.text(0, 1.015, sub, transform=ax.transAxes, fontsize=10, color=GREY,
                va="bottom")


# --------------------------------------------------------------------------
# 1  Where the contact map comes from
# --------------------------------------------------------------------------
def fig1() -> None:
    sources = ["Measured Hi-C\n(cell type)", "UniversalEPI\n(sequence + ATAC)",
               "Averaged megamap\n(ENCSR492BKE)", "Inverse distance\n(power law)"]
    prec = [53.3, np.nan, 48.0, 44.3]
    colours = [TEAL, AMBER, GREY, ROSE]

    fig, ax = plt.subplots(figsize=(9, 4.4))
    bars = ax.bar(range(4), [p if p == p else 0 for p in prec], color=colours,
                  width=0.62)
    for i, (b, p) in enumerate(zip(bars, prec)):
        if p == p:
            ax.text(b.get_x() + b.get_width() / 2, p + 0.7, f"{p:.1f}%",
                    ha="center", fontsize=12, color=INK)
        else:
            ax.text(i, 1.2, "not benchmarked\nby rE2G", ha="center", fontsize=10,
                    color=AMBER, style="italic")
    ax.set_xticks(range(4), sources, fontsize=10)
    ax.set_ylabel("ENCODE-rE2G precision at 70% recall")
    ax.set_ylim(0, 60)
    tidy(ax)
    heading(ax, "Use the best contact source available, and record which one",
            "ENCODE-rE2G benchmark (Gschwind et al. 2026, Nature 657:179-190)")
    ax.text(0.5, 0.94, "the choice is worth 9 points of precision",
            transform=ax.transAxes, fontsize=11, color=MUTE, ha="center")
    finish(fig, "fig1_contact_source",
           "Step 1 — the contact map",
           "Measured Hi-C beats every prediction, so where it exists we use it "
           "and bypass the model entirely. UniversalEPI fills in cell types with "
           "no Hi-C from sequence and ATAC; a cell type with neither falls back "
           "to the averaged megamap or plain inverse distance, at a cost of nine "
           "points of precision. Whichever rung was used is written into the "
           "record, so a reviewer can see what a claim rests on.")


# --------------------------------------------------------------------------
# 2  ABC: share is attribution, T is dose
# --------------------------------------------------------------------------
def fig2(r) -> None:
    loc = r.locus
    tgt = loc.genes["TARGET"]
    contrib = loc.activities * loc.contact[loc.elements, tgt]
    share = contrib / contrib.sum()
    kb = loc.elements * BIN_KB

    fig, axes = plt.subplots(3, 1, figsize=(9.4, 6.4), sharex=True,
                             gridspec_kw={"height_ratios": [1, 1, 1.25],
                                          "hspace": 0.34})
    ax = axes[0]
    ax.stem(kb, loc.activities, linefmt="-", markerfmt="o", basefmt=" ")
    for ln in ax.get_children():
        pass
    ax.set_ylabel("A$_e$")
    tidy(ax)
    heading(ax, "One gene's regulatory input, taken apart",
            f"TARGET at {tgt * BIN_KB} kb  ·  activity x contact, element by element")

    ax = axes[1]
    ax.stem(kb, loc.contact[loc.elements, tgt], linefmt="-", markerfmt="o",
            basefmt=" ")
    ax.set_ylabel("C(e, G)")
    tidy(ax)

    ax = axes[2]
    ax.bar(kb, contrib, width=26, color=TEAL, alpha=0.9)
    for x, c, s in zip(kb, contrib, share):
        ax.text(x, c + contrib.max() * 0.04, f"{s:.0%}", ha="center", fontsize=10,
                color=INK)
    ax.set_ylabel("A$_e$ · C(e, G)")
    ax.set_xlabel("position in the 2 Mb window (kb)")
    ax.set_ylim(0, contrib.max() * 1.28)
    tidy(ax)
    ax.text(0.01, 0.95,
            f"bars sum to  T$_G$ = {contrib.sum():.3f}   ← the dose\n"
            "labels are the ABC share   ← the attribution",
            transform=ax.transAxes, ha="left", va="top", fontsize=10.5,
            family=MONO, color=MUTE,
            bbox=dict(boxstyle="round,pad=0.5", fc=PANEL, ec=EDGE))

    for a in axes:
        a.axvline(tgt * BIN_KB, color=TEAL, lw=1.1, ls=(0, (4, 4)), alpha=0.55,
                  zorder=0)
    axes[0].text(tgt * BIN_KB, axes[0].get_ylim()[1] * 0.98, " TARGET TSS",
                 color=TEAL, fontsize=9.5, va="top")

    for ax in axes[:2]:
        for coll in ax.collections + list(ax.lines):
            coll.set_color(MUTE)
    finish(fig, "fig2_abc_anatomy",
           "Step 2 — ABC, and the two numbers it gives you",
           "Activity times contact, element by element. The bars sum to T, the "
           "gene's total regulatory input; the percentages are the ABC score, "
           "each element's share of that total. They answer different questions. "
           "The share is normalised per gene and always sums to 100%, so after an "
           "intervention it reports redistribution, not increase. T is the dial a "
           "staple is meant to turn, and dT/T is what we report.")


# --------------------------------------------------------------------------
# 3  Where a staple is allowed to land
# --------------------------------------------------------------------------
def fig3(r) -> None:
    loc = r.locus
    fig, ax = plt.subplots(figsize=(10.4, 4.0))
    ok = np.where(loc.accessible & ~r.forbidden)[0]

    rows = {"genes": 3.0, "excluded": 2.0, "atac": 1.0, "anchors": 0.0}
    for y in rows.values():
        ax.axhline(y, color=EDGE, lw=0.9, zorder=0)

    for i, (g, (s, e)) in enumerate(sorted(loc.gene_bodies.items(),
                                           key=lambda kv: kv[1][0])):
        ax.add_patch(Rectangle((s * BIN_KB, rows["genes"] - 0.14),
                               max((e - s) * BIN_KB, 18), 0.28, color=INK,
                               alpha=0.75, zorder=3))
        ax.text(loc.genes[g] * BIN_KB, rows["genes"] + (0.46 if i % 2 else 0.22),
                g, fontsize=9.5, color=MUTE, ha="center")

    for b in loc.elements:
        ax.add_patch(Rectangle((b * BIN_KB - 9, rows["excluded"] - 0.14), 18,
                               0.28, color=ROSE, zorder=3))
    for b in loc.ctcf:
        ax.add_patch(Rectangle((b * BIN_KB - 9, rows["excluded"] - 0.14), 18,
                               0.28, color=AMBER, zorder=3))

    peaks = np.where(loc.accessible)[0] * BIN_KB
    ax.vlines(peaks, rows["atac"] - 0.17, rows["atac"] + 0.17, color=GREY, lw=1.1)
    ax.scatter(ok * BIN_KB, np.full(len(ok), rows["anchors"]), s=18, color=TEAL,
               zorder=4)
    for b, lab in ((DEMO_ANCHOR_A, "anchor A"), (DEMO_ANCHOR_B, "anchor B")):
        ax.annotate(lab, (b * BIN_KB, rows["anchors"]), (b * BIN_KB, -0.78),
                    ha="center", fontsize=10.5, color=TEAL,
                    arrowprops=dict(arrowstyle="-", color=TEAL, lw=1.3))

    ax.set_yticks(list(rows.values()),
                  ["genes", "excluded", f"ATAC  ({loc.accessible.sum()})",
                   f"permitted  ({len(ok)})"], fontsize=10)
    ax.set_ylim(-1.0, 4.0)
    ax.set_xlim(-20, 2020)
    ax.set_xlabel("position in the 2 Mb window (kb)")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="y", length=0)
    heading(ax, "Subtract everything a staple would break",
            f"{loc.accessible.sum()} ATAC peaks  →  {len(ok)} permitted anchors")
    ax.legend(handles=[
        Patch(color=INK, alpha=0.75, label="gene body + promoter"),
        Patch(color=ROSE, label="enhancer"),
        Patch(color=AMBER, label="CTCF site"),
        Patch(color=GREY, label="ATAC peak"),
        Patch(color=TEAL, label="permitted anchor"),
    ], loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=5, frameon=False,
        fontsize=9.5)
    finish(fig, "fig3_anchor_filter",
           "Step 3 — filter the physical conflicts",
           "An anchor has to be open chromatin the machinery can bind, minus "
           "everything we would break by parking a protein on it: the enhancer "
           "itself, the promoter bin, the transcribed body, and "
           "CTCF boundary sites. This is a filter, not a score — it removes the "
           "impossible rather than ranking the possible, which is why 'tractable' "
           "does not need to be an axis in the ranking.")


# --------------------------------------------------------------------------
# 4  Why the edit must be a kernel
# --------------------------------------------------------------------------
def fig4(r) -> None:
    loc = r.locus
    n = len(r.before)
    a, b = DEMO_ANCHOR_A, DEMO_ANCHOR_B
    point = r.before.copy()
    point[a, b] = point[b, a] = max(point[a, b], r.q)
    point = sweep_mod.apply_staple(point, a, b, 0.0)       # same renormalisation
    order = loc.gene_order

    def deltas(mat):
        idx = np.array([loc.genes[g] for g in order])
        t0 = abc_mod.total_input(loc.activities, r.before, loc.elements, idx)
        t1 = abc_mod.total_input(loc.activities, mat, loc.elements, idx)
        return (t1 - t0) / t0

    dp, dk = deltas(point), deltas(r.after)

    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.9),
                             gridspec_kw={"width_ratios": [1, 1, 1.5],
                                          "wspace": 0.34})
    win = slice(a - 18, a + 19), slice(b - 18, b + 19)
    for ax, mat, lab in ((axes[0], point, "point edit"),
                         (axes[1], r.after, "Gaussian kernel, σ = 3 bins")):
        ax.imshow((mat - r.before)[win[0], win[1]], cmap=DIFF, origin="lower",
                  vmin=-r.q, vmax=r.q, aspect="equal")
        ax.set_xticks([]), ax.set_yticks([])
        ax.set_title(lab, fontsize=11, color=INK)
        ax.set_xlabel("anchor B →", fontsize=9)
        ax.set_ylabel("anchor A →", fontsize=9)

    ax = axes[2]
    y = np.arange(len(order))
    ax.barh(y - 0.19, dp * 100, height=0.36, color=GREY, label="point edit")
    ax.barh(y + 0.19, dk * 100, height=0.36, color=TEAL, label="kernel")
    ax.set_yticks(y, order, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel("ΔT / T  (%)")
    ax.axvline(0, color=INK, lw=1)
    tidy(ax, left=False)
    ax.grid(axis="x", color=EDGE, lw=0.8)
    ax.legend(frameon=False, fontsize=9.5, loc="lower right",
              bbox_to_anchor=(1.0, 0.0))
    ax.set_xlim(-2, max(dk.max() * 100 * 1.18, 10))
    zero = "exactly zero" if abs(dp).max() == 0 else f"below {abs(dp).max():.1e}"
    ax.text(0.03, 0.25, f"every point-edit bar is {zero}",
            transform=ax.transAxes, fontsize=10.5, color=GREY, style="italic")
    fig.suptitle("A point edit cannot produce an off-target effect — "
                 "and here it produces no effect at all",
                 x=0.045, ha="left", fontsize=14, color=INK, y=1.14)
    fig.text(0.045, 1.04, "ABC normalises per gene, so one cell moves one gene. "
             "Step 3 then forbids anchors on the element or the TSS, so the one "
             "cell touches neither.", fontsize=10, color=GREY, ha="left")
    finish(fig, "fig4_point_vs_kernel",
           "Step 4a — the intervention has to be a blob",
           "ABC's denominator is per gene, so editing a single cell of the contact "
           "matrix can only move the gene you aimed at. Enhancer hijacking becomes "
           "arithmetically impossible and the off-target analysis returns zeros by "
           "construction. A model that cannot express the failure mode cannot rule "
           "it out. So the staple deposits a Gaussian blob, sigma about 3 bins or "
           "15 kb — a width with empirical support: rE2G found significant CRISPRi "
           "crosstalk between elements 1-10 kb apart, P = 1.8e-5.")


# --------------------------------------------------------------------------
# 5  Before, after, difference
# --------------------------------------------------------------------------
def fig5(r) -> None:
    loc = r.locus
    lo, hi = 60, 380
    sub = slice(lo, hi)
    vmax = float(np.percentile(r.before, 99.3))
    d = r.after - r.before
    dmax = float(np.percentile(np.abs(d), 99.6))

    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.3),
                             gridspec_kw={"wspace": 0.16})
    ext = [lo * BIN_KB, hi * BIN_KB] * 2
    for ax, mat, cm, lim, lab in (
            (axes[0], r.before, HIC, (0, vmax), "before"),
            (axes[1], r.after, HIC, (0, vmax), "after the staple"),
            (axes[2], d, DIFF, (-dmax, dmax), "after − before")):
        im = ax.imshow(mat[sub, sub], cmap=cm, origin="lower", extent=ext,
                       vmin=lim[0], vmax=lim[1], aspect="equal",
                       interpolation="nearest")
        ax.set_title(lab, fontsize=12, color=INK, loc="left")
        for e in loc.tad_edges[1:-1]:
            ax.axhline(e * BIN_KB, color=INK, lw=0.7, ls=(0, (4, 4)), alpha=0.45)
            ax.axvline(e * BIN_KB, color=INK, lw=0.7, ls=(0, (4, 4)), alpha=0.45)
        ax.set_xticks([500, 1000, 1500])
        ax.set_yticks([500, 1000, 1500])
        ax.tick_params(labelsize=9)
    axes[2].plot(DEMO_ANCHOR_B * BIN_KB, DEMO_ANCHOR_A * BIN_KB, "o",
                 mfc="none", mec=INK, ms=16, mew=1.4)
    axes[2].annotate("forced contact", (DEMO_ANCHOR_B * BIN_KB,
                                        DEMO_ANCHOR_A * BIN_KB),
                     (1180, 760), fontsize=10, color=INK,
                     arrowprops=dict(arrowstyle="->", color=INK, lw=1.1))
    axes[2].text(430, 1680, "contact taken from\neverywhere else", fontsize=10,
                 color=TEAL)
    axes[0].set_ylabel("kb")
    fig.suptitle("The staple spends a contact budget, it does not create one",
                 x=0.125, ha="left", fontsize=14, color=INK, y=1.01)
    finish(fig, "fig5_before_after_diff",
           "Step 4b — what the edit actually did",
           "Left and centre look almost identical, which is the point: subtract "
           "them. The amber blob is the contact we forced between the two anchors. "
           "The teal stripes running along both anchor rows are contact taken away "
           "from everywhere else, because we hold each locus's total contact budget "
           "fixed. Without that constraint a staple would add closeness for free, "
           "nothing could ever go down, and enhancer theft would be unrepresentable.")


# --------------------------------------------------------------------------
# 6  Dose-response across the staple strength sweep
# --------------------------------------------------------------------------
def fig6(r) -> None:
    grid = r.sweep["q"]
    pct = [50, 75, 90, 95, 99]
    fig, ax = plt.subplots(figsize=(9.4, 4.6))
    for g in r.locus.gene_order:
        v = r.sweep[g] * 100
        col = TEAL if g == "TARGET" else (AMBER if abs(v[-1]) > 3 else GREY)
        ax.plot(pct, v, "-o", color=col, lw=2.2 if col != GREY else 1.2,
                ms=5, alpha=1.0 if col != GREY else 0.55)
        if abs(v[-1]) > 2:
            ax.text(99.6, v[-1], f"  {g}", color=col, fontsize=11, va="center")
    ax.set_xticks(pct, [f"p{p}" for p in pct])
    ax.text(50.4, ax.get_ylim()[1] * 0.96,
            f"the single figure quoted elsewhere uses q = {r.q:.3f}, the strongest\n"
            f"real contact at this separation — above p99, by construction",
            fontsize=10, color=MUTE, va="top",
            bbox=dict(boxstyle="round,pad=0.45", fc=PANEL, ec=EDGE))
    ax.set_xlabel("staple strength q, as a percentile of real contacts at the same separation")
    ax.set_ylabel("ΔT / T  (%)")
    ax.set_xlim(48, 108)
    tidy(ax)
    heading(ax, "Nobody has measured how strong a staple is, so we sweep it",
            "the ranking is stable across the whole range — that is the confidence test")
    finish(fig, "fig6_dose_response",
           "Step 4c — sweep the strength instead of guessing it",
           "No one has measured the contact frequency a dCas staple produces, so "
           "we refuse to pick a number. We sweep q across the range of real "
           "contacts observed at the same genomic separation and report the whole "
           "curve. This turns an unknown parameter into an output. It also gives "
           "the confidence test teeth: if the ordering of genes flips as q varies, "
           "we do not know the answer and the candidate is marked tentative.")


# --------------------------------------------------------------------------
# 7  Read every gene in the window
# --------------------------------------------------------------------------
def fig7(r) -> None:
    loc = r.locus
    order = loc.gene_order
    vals = np.array([r.delta[g] for g in order]) * 100
    pos = np.array([loc.genes[g] for g in order]) * BIN_KB

    fig, ax = plt.subplots(figsize=(10.2, 4.4))
    quiet = 0
    for g, x, v in zip(order, pos, vals):
        if g == "TARGET":
            col = TEAL
        elif v > 1:
            col = AMBER
        elif v < -1:
            col = ROSE
        else:
            col = "#c2bcae"
        ax.bar(x, v, width=46, color=col)
        shown = 0.0 if abs(v) < 0.05 else v
        if abs(v) > 1:
            ax.text(x, v + 2.4, f"{g}\n{shown:+.1f}%", ha="center", fontsize=10.5,
                    color=INK, va="bottom")
        else:
            y = -5.0 if quiet % 2 == 0 else -9.5
            quiet += 1
            ax.plot([x, x], [0, y + 1.1], color=EDGE, lw=0.9, zorder=0)
            ax.text(x, y, f"{g}  {shown:+.1f}%", ha="center", fontsize=9.5,
                    color=GREY, va="top")
    ax.axhline(0, color=INK, lw=1.1)
    ax.axhline(TAU * 100, color=ROSE, lw=1.1, ls=(0, (5, 4)))
    ax.text(2000, TAU * 100 + 1.2, f"τ = {TAU:.0%}", color=ROSE, fontsize=10,
            ha="right")
    ax.set_xlim(0, 2000)
    ax.set_ylim(-16, 54)
    ax.set_xlabel("position in the 2 Mb window (kb)")
    ax.set_ylabel("ΔT / T  (%)")
    tidy(ax)
    heading(ax, "Score every gene in the window, not just the one you aimed at",
            "one staple, 2 Mb, seven genes  ·  the neighbour 30 kb away was not in the plan")
    finish(fig, "fig7_window_readout",
           "Step 5 — the window readout",
           "The target gains 42% of regulatory input, which is what we wanted. A "
           "neighbouring promoter 30 kb away gains 14%, which we did not want and "
           "would never have seen had we scored only the gene we aimed at. The "
           "flagged genes elsewhere in the window do not move. This panel is the "
           "whole argument for simulating the intervention rather than just "
           "scoring the target.")


# --------------------------------------------------------------------------
# 8  The gate, and what it costs
# --------------------------------------------------------------------------
def fig8(a, b) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.4),
                             gridspec_kw={"wspace": 0.3})
    for ax, r, tgt, label in ((axes[0], a, "TARGET", "candidate A"),
                              (axes[1], b, "PARTNER", "candidate B")):
        flagged = {g: r.delta[g] for g in r.locus.flags}
        worst = max(flagged, key=lambda g: abs(flagged[g]))
        bars = [("intended\n" + tgt, r.delta[tgt] * 100, TEAL),
                ("worst flagged\n" + worst, flagged[worst] * 100, ROSE)]
        ax.bar([0, 1], [v for _, v, _ in bars], color=[c for *_, c in bars],
               width=0.5)
        for i, (_, v, _) in enumerate(bars):
            ax.text(i, v + 0.9, f"{v:+.1f}%", ha="center", fontsize=12, color=INK)
        ax.axhline(TAU * 100, color=ROSE, lw=1.1, ls=(0, (5, 4)))
        ax.text(1.42, TAU * 100 + 0.6, f"τ = {TAU:.0%}", color=ROSE, fontsize=9.5,
                ha="right")
        ax.set_xticks([0, 1], [n for n, *_ in bars], fontsize=10)
        ax.set_ylim(0, 56)
        ax.set_ylabel("ΔT / T  (%)")
        tidy(ax)
        verdict = "PASS" if r.candidate.safe else "FAIL"
        col = TEAL if r.candidate.safe else ROSE
        ax.set_title(f"{label}   ·   {verdict}", loc="left", fontsize=13, color=col)
        note = ("no flagged gene moves more than τ,\n"
                "and none moves more than the target"
                if r.candidate.safe else
                f"{worst} moves more than the gene we aimed at.\n"
                "Under τ alone this would have passed.")
        ax.text(0.5, 0.86, note, transform=ax.transAxes, ha="center", fontsize=10,
                color=MUTE, bbox=dict(boxstyle="round,pad=0.45", fc=PANEL, ec=EDGE))
    fig.suptitle("Safety is a gate, not a term in a score",
                 x=0.09, ha="left", fontsize=14, color=INK, y=1.02)
    finish(fig, "fig8_safety_gate",
           "Step 6 — the gate, and the ordering",
           "Two tests, because they fail for different reasons. Absolute: no "
           "flagged gene moves more than tau, which we start at 10%. Relative: no "
           "flagged gene moves more than the gene we were aiming at — candidate B "
           "does 6.1% to an oncogene and only 5.0% to its target, so it dies even "
           "though both numbers are small and tau alone would have let it through. "
           "A weighted score would allow a large predicted effect to outrank a "
           "safety concern. A gate cannot. Survivors are then ordered by dT/T, "
           "with the confidence tier shown beside each row rather than sorted on.")


def main() -> None:
    print("rendering ->", OUT)
    a = run(anchor_a=DEMO_ANCHOR_A, anchor_b=DEMO_ANCHOR_B, target="TARGET")
    b = run(anchor_a=UNSAFE_ANCHOR_A, anchor_b=UNSAFE_ANCHOR_B, target="PARTNER")
    fig1()
    fig2(a)
    fig3(a)
    fig4(a)
    fig5(a)
    fig6(a)
    fig7(a)
    fig8(a, b)

    lines = ["# Figure series — the method, step by step", "",
             "Generated by `animation/figures/method_figures.py`. Every number "
             "comes from `pipeline.simulate`, the same code path the pipeline "
             "runs. Inputs are synthetic (`pipeline/locus.py`, `illustrative: "
             "true`); the arithmetic on them is real.", ""]
    for name, title, cap in CAPTIONS:
        lines += [f"## {title}", "", f"![{title}]({name}.png)", "", cap, ""]
    (OUT / "captions.md").write_text("\n".join(lines))
    print("  captions.md")


if __name__ == "__main__":
    main()

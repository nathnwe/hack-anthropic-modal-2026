"""Render the off-target simulation diagram.

What the pipeline actually simulates: the staple lands exactly where we aimed,
and we ask which *other* genes in the window change as a result. The figure is
computed, not drawn -- the contact edit, the row-sum renormalisation and the
per-gene totals are the real arithmetic from the spec, run on a synthetic
contact matrix with realistic statistics. Only the inputs are invented.

Three reads, top to bottom:
  1  before and after, with the staple arc
  2  the difference map -- what the edit actually did to the window
  3  change in total regulatory input per gene, aligned to the map above

Palette and box helpers are shared with pipeline_diagram so the two figures
stay visually consistent; that module's own figure only renders under main.

Run:
    uv run --with matplotlib --with numpy \
        python animation/figures/offtarget_diagram.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Rectangle

from pipeline_diagram import (  # noqa: E402  (shared palette and box helpers)
    PANEL,
    INK,
    TEAL,
    AMBER,
    ROSE,
    GREY,
    MUTE,
    EDGE,
    MONO,
    HIC_CMAP,
    band_label,
    T_SUB,
    T_TINY,
)

# ------------------------------------------------------------------ canvas
W_IN, H_IN = 8.0, 4.3
XMAX = 100.0
YMAX = XMAX * H_IN / W_IN
MAP_X0, MAP_X1 = 1.0, 99.0

DIFF_CMAP = LinearSegmentedColormap.from_list(
    "diff",
    ["#0b5f96", "#4aa8de", "#c4e3f6", "#fffdf7", "#f0cf93", "#c98d2a", "#8a5410"],
)

RNG = np.random.default_rng(7)

# -------------------------------------------------------------- the window
N = 180  # 5 kb bins across 900 kb
TAD_EDGES = [0, 62, 118, N]
ANCHOR_A, ANCHOR_B = 34, 150  # where the two halves of the staple sit
SIGMA = 3.0  # kernel width in bins, ~15 kb

ELEMENTS = [(28, 1.5), (34, 3.1), (60, 1.1), (100, 1.4), (130, 0.9), (160, 1.8)]
GENES = [
    ("neighbour", 40),
    ("distal", 96),
    ("TARGET", 150),
    ("neighbour", 146),
    ("distal", 172),
]


def linear_contacts() -> np.ndarray:
    """Polymer decay + TAD blocks + CTCF corner dots, in linear space."""
    i, j = np.indices((N, N))
    m = 1.0 / (1.0 + np.abs(i - j)) ** 0.95
    for a, b in zip(TAD_EDGES[:-1], TAD_EDGES[1:]):
        m = np.where((i >= a) & (i < b) & (j >= a) & (j < b), m * 2.4, m)
    for a, b in zip(TAD_EDGES[:-1], TAD_EDGES[1:]):
        for ci, cj in ((a + 2, b - 3), (b - 3, a + 2)):
            m += 0.5 * np.exp(-(((i - ci) ** 2 + (j - cj) ** 2) / 8.0))
    m = m * RNG.lognormal(0.0, 0.20, size=m.shape)
    return (m + m.T) / 2.0


def staple(m: np.ndarray, q: float) -> np.ndarray:
    """C' = max(C, q*K), then preserve every row's total contact budget."""
    i, j = np.indices((N, N))
    k = np.exp(-(((i - ANCHOR_A) ** 2 + (j - ANCHOR_B) ** 2) / (2 * SIGMA**2)))
    k = q * np.maximum(k, k.T)
    out = np.maximum(m, k)
    for _ in range(4):  # symmetric row-budget restore
        out = out * (m.sum(1) / out.sum(1))[:, None]
        out = (out + out.T) / 2.0
    return out


def total_input(m: np.ndarray) -> dict[int, float]:
    """T_G = sum over elements of activity x contact to that gene's TSS."""
    return {g: sum(a * m[e, g] for e, a in ELEMENTS) for _, g in GENES}


BEFORE = linear_contacts()
# Yardstick for the forced contact: as tight as the strongest real loop at
# 100 kb anywhere in this window. Describable, and not a free parameter.
Q = float(np.diagonal(BEFORE, offset=20).max())
AFTER = staple(BEFORE, Q)
DIFF = AFTER - BEFORE

T0, T1 = total_input(BEFORE), total_input(AFTER)
DELTA = {g: (T1[g] - T0[g]) / T0[g] for _, g in GENES}

VMAX = float(np.percentile(np.log1p(BEFORE * 60), 99.2))
DMAX = float(np.percentile(np.abs(DIFF), 99.0))


def gx(b: float) -> float:
    """Genomic bin -> x on the shared axis, so every panel lines up."""
    return MAP_X0 + (MAP_X1 - MAP_X0) * b / N


def tri(ax, mat, x0, x1, ybase, ytop, *, cmap, vmin, vmax, boundaries=False):
    bw, bh = x1 - x0, ytop - ybase
    tr = (
        mtransforms.Affine2D()
        .translate(-N / 2, -N / 2)
        .rotate_deg(-45)
        .scale(1 / np.sqrt(2))
        .translate(N / 2, 0)
        .scale(bw / N, bh / (N / 2.0))
        .translate(x0, ybase)
        + ax.transData
    )
    im = ax.imshow(
        mat,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        origin="lower",
        extent=(0, N, 0, N),
        interpolation="bilinear",
        zorder=2,
    )
    im.set_transform(tr)
    im.set_clip_path(Rectangle((x0, ybase), bw, bh, transform=ax.transData))
    if boundaries:
        for b in TAD_EDGES[1:-1]:
            bx = x0 + bw * b / N
            ax.plot(
                [bx, bx],
                [ybase, ybase + bh * 0.55],
                lw=1.1,
                ls=(0, (3, 2.5)),
                color=AMBER,
                zorder=4,
            )




def falloff(offsets_kb):
    """How much a probe promoter that far from the target TSS picks up."""
    out = []
    for kb in offsets_kb:
        g = ANCHOR_B - int(round(kb / 5))
        t0 = sum(a * BEFORE[e, g] for e, a in ELEMENTS)
        t1 = sum(a * AFTER[e, g] for e, a in ELEMENTS)
        out.append((t1 - t0) / t0)
    return np.array(out)


OFFSETS = np.arange(0, 61, 5)
FALLOFF = falloff(OFFSETS)


def build(ax):
    # ===================================================== 1  before / after
    band_label(ax, 1, 52.9, "1  THE STAPLE LANDS EXACTLY WHERE WE AIMED")

    for x0, x1, title, mat, arc in [
        (1, 26, "before", BEFORE, False),
        (30, 55, "after", AFTER, True),
    ]:
        ax.text(x0, 50.9, title, fontsize=T_SUB, color=INK, weight="bold",
                va="center")
        tri(ax, np.log1p(mat * 60), x0, x1, 42.0, 49.6,
            cmap=HIC_CMAP, vmin=0, vmax=VMAX, boundaries=True)
        ax.plot([x0, x1], [42.0, 42.0], color=INK, lw=1.1, alpha=0.5,
                zorder=5)
        if arc:
            pa = x0 + (x1 - x0) * ANCHOR_A / N
            pb = x0 + (x1 - x0) * ANCHOR_B / N
            th = np.linspace(0, np.pi, 60)
            mid, half = (pa + pb) / 2, (pb - pa) / 2
            ax.plot(mid - half * np.cos(th), 41.8 - 1.8 * np.sin(th),
                    color=TEAL, lw=1.8, zorder=6)
            for p in (pa, pb):
                ax.add_patch(plt.Circle((p, 41.8), 0.5, color=TEAL, zorder=7))
            ax.text(mid, 38.8, "staple", fontsize=T_TINY, color=TEAL,
                    ha="center", va="center", weight="bold")

    # --- how far the accident reaches ------------------------------------
    ax.add_patch(FancyBboxPatch(
        (62, 38.6), 37, 12.2, boxstyle="round,pad=0,rounding_size=0.9",
        facecolor=PANEL, edgecolor=EDGE, lw=1.2, zorder=2))
    ax.text(64, 49.6, "How far the accident reaches", fontsize=T_SUB,
            color=INK, weight="bold", va="center", zorder=5)

    px0, px1, py0, py1 = 67.0, 96.0, 41.5, 47.6
    fmax = max(FALLOFF.max(), 1e-9)
    fx = px0 + (px1 - px0) * OFFSETS / 60.0
    fy = py0 + (py1 - py0) * FALLOFF / fmax
    ax.plot([px0, px1], [py0, py0], color=INK, lw=1.0, alpha=0.5, zorder=4)
    ax.plot(fx, fy, color=AMBER, lw=1.8, zorder=5)
    ax.fill_between(fx, py0, fy, color=AMBER, alpha=0.16, lw=0, zorder=4)
    for kb in (0, 20, 40, 60):
        x = px0 + (px1 - px0) * kb / 60.0
        ax.text(x, 40.5, str(kb), fontsize=T_TINY, color=GREY, ha="center",
                va="center", zorder=5)
    ax.text(px1, 39.4, "kb from the target TSS", fontsize=T_TINY, color=GREY,
            ha="right", va="center", zorder=5)
    j = int(np.argmin(np.abs(OFFSETS - 20)))
    ax.add_patch(plt.Circle((fx[j], fy[j]), 0.55, color=AMBER, zorder=6))
    ax.text(fx[j] + 1.4, fy[j] + 1.4,
            f"a promoter 20 kb away\npicks up {FALLOFF[j]:+.0%}",
            fontsize=T_TINY, color=AMBER, va="center", linespacing=1.4,
            zorder=6)
    ax.text(px0 - 0.8, fy[0], f"{FALLOFF[0]:+.0%}", fontsize=T_TINY,
            color=TEAL, ha="right", va="center", weight="bold", zorder=6)
    ax.text(px0 - 0.8, py0, "0", fontsize=T_TINY, color=GREY, ha="right",
            va="center", zorder=6)
    ax.text(px0 + 0.8, 48.7, "the target", fontsize=T_TINY, color=TEAL,
            va="center", zorder=6)

    # ===================================================== 2  difference map
    band_label(ax, 1, 37.0, "2  SUBTRACT \u2014 WHAT THE EDIT ACTUALLY DID")
    tri(ax, np.clip(np.sign(DIFF) * (np.abs(DIFF) / DMAX) ** 0.45, -1, 1),
        MAP_X0, MAP_X1, 23.0, 35.2, cmap=DIFF_CMAP, vmin=-1, vmax=1)
    ax.plot([MAP_X0, MAP_X1], [23.0, 23.0], color=INK, lw=1.2, alpha=0.55,
            zorder=5)
    ax.text(gx((ANCHOR_A + ANCHOR_B) / 2), 31.6, "forced contact",
            fontsize=T_TINY, color="#8a5410", ha="center", va="center",
            weight="bold")
    ax.text(gx(96), 25.4, "contact taken from everywhere else",
            fontsize=T_TINY, color="#0b5f96", ha="center", va="center")
    for lbl, col, x in [("gained", "#8a5410", 74), ("lost", "#0b5f96", 88)]:
        ax.add_patch(Rectangle((x, 36.3), 1.6, 1.2, color=col, zorder=5))
        ax.text(x + 2.4, 36.9, lbl, fontsize=T_TINY, color=MUTE,
                va="center", zorder=5)

    # ===================================================== 3  per-gene read
    band_label(ax, 1, 20.0, "3  READ EVERY GENE IN THE WINDOW")

    zero = 12.4
    ax.plot([MAP_X0, MAP_X1], [zero, zero], color=INK, lw=1.2, zorder=4)
    for e, _a in ELEMENTS:
        ax.add_patch(Rectangle((gx(e) - 0.5, 21.4), 1.0, 0.8, color=GREY,
                               zorder=5))
    ax.text(MAP_X1, 21.8, "elements", fontsize=T_TINY, color=GREY,
            ha="right", va="center")
    ax.text(MAP_X0, 18.6, "\u0394 total regulatory input per gene",
            fontsize=T_TINY, color=MUTE, va="center")

    worst = max(abs(v) for v in DELTA.values())
    order = [g for _n, g in sorted(GENES, key=lambda t: t[1])]
    xs = [gx(g) for g in order]
    for _ in range(60):                       # relax until labels stop touching
        for i in range(len(xs) - 1):
            gap = xs[i + 1] - xs[i]
            if gap < 11.0:
                push = (11.0 - gap) / 2
                xs[i] -= push
                xs[i + 1] += push
        xs = [min(max(x, 7.0), 93.0) for x in xs]
    label_x = dict(zip(order, xs))
    for name, g in GENES:
        d = DELTA[g]
        if g == ANCHOR_B:
            col, tag = TEAL, "intended"
        elif d > 0.05:
            col, tag = AMBER, "off-target gain"
        elif d < -0.015:
            col, tag = ROSE, "off-target loss"
        else:
            col, tag = "#c2bcae", "unchanged"
        h = (d / worst) * 5.0
        ax.add_patch(Rectangle((gx(g) - 0.8, zero if h >= 0 else zero + h),
                               1.6, max(abs(h), 0.18), color=col, zorder=5))
        ax.text(gx(g), zero + h + (0.8 if h >= 0 else -0.8), f"{d:+.0%}",
                fontsize=T_TINY, color=col, ha="center",
                va="bottom" if h >= 0 else "top", family=MONO, zorder=6)
        lx = label_x[g]
        if abs(lx - gx(g)) > 0.4:
            ax.plot([gx(g), lx], [7.5, 7.2], color="#c2bcae", lw=0.9,
                    zorder=4)
        ax.text(lx, 6.8, name, fontsize=T_TINY, color=INK, ha="center",
                va="center", weight="bold" if name == "TARGET" else "normal")
        ax.text(lx, 5.6, tag, fontsize=T_TINY, color=col, ha="center",
                va="center")

    # ===================================================== 4  the numbers
    ax.plot([MAP_X0, MAP_X1], [4.4, 4.4], color=EDGE, lw=1.0, zorder=3)
    for k, (big, cap) in enumerate([
        ("1.57", "genes per enhancer,\nand up to 34"),
        ("P = 1.8×10⁻⁵", "CRISPRi crosstalk\nat 1–10 kb"),
        ("16–39%", "super-additive\nat 3–89 kb"),
        ("53.2%", "rE2G precision out\nof cell type"),
        ("r = −0.45", "score vs effect size:\nwe rank, not dose"),
    ]):
        x = 1 + k * 19.8
        ax.plot([x, x], [0.8, 3.8], color=EDGE, lw=1.2, zorder=3)
        ax.text(x + 1.4, 3.3, big, fontsize=T_SUB, color=AMBER, family=MONO,
                weight="bold", va="center", zorder=5)
        ax.text(x + 1.4, 1.7, cap, fontsize=T_TINY, color=MUTE, va="center",
                linespacing=1.45, zorder=5)


REPO = Path(__file__).resolve().parents[2]


def main(out=None):
    out = str(out or REPO / "docs/figures/offtarget.png")
    fig = plt.figure(figsize=(W_IN, H_IN), dpi=320)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, XMAX)
    ax.set_ylim(0, YMAX)
    ax.set_aspect("equal")
    ax.axis("off")
    build(ax)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=320)
    print("wrote", out)
    for name, g in GENES:
        print(f"  {name:>10} @{g:4d}  {DELTA[g]:+.1%}")


if __name__ == "__main__":
    main()

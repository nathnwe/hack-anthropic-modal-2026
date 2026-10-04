"""Render the technical flow diagram: inputs -> models -> contact -> ranking.

Two bands. The top band builds one contact map from whatever evidence the
cell type has. The bottom band uses that map: place anchors, force a contact,
re-score the window, grade it on four axes, rank.

Node fill marks who owns the module: blue for `pipeline/`, amber for
`scoring/`. Thumbnails (sequence, ATAC, Hi-C, the dose-response curve, the
window bars) are drawn from synthetic data with realistic statistics --
polymer distance decay, TAD blocks, CTCF-anchored corner dots, peaky ATAC
signal. Nothing here is pipeline output and the slide says so.

Output is sized for a slide: the figure is deliberately small in inches so
that type baked into it lands above 24 px once the image is displayed at
about 1,550 px wide.

Run:
    uv run --with matplotlib --with numpy \
        python animation/figures/pipeline_diagram.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

# ------------------------------------------------------------------ palette
BG = "#fffefb"
PANEL = "#f4f1e9"
INK = "#14201e"
TEAL = "#1e9ae0"  # light blue (deck accent)
AMBER = "#9a6410"
ROSE = "#a33a3a"
GREY = "#8a8275"
MUTE = "#566460"
EDGE = "#ded8c9"

HIC_CMAP = LinearSegmentedColormap.from_list(
    "hic", ["#fffdf7", "#f7e7c2", "#e3a954", "#b8771b", "#7a3d14", "#4a1f0c"]
)
MONO = ["DejaVu Sans Mono", "Courier New", "monospace"]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial"],
    "figure.facecolor": BG,
    "savefig.facecolor": BG,
})

RNG = np.random.default_rng(11)

# -------------------------------------------------------- synthetic contacts
N = 180                       # bins; read as 5 kb bins across a 900 kb window
TAD_EDGES = [0, 62, 118, N]
E_BIN, P_BIN = 34, 150


def contact_matrix(loops: list[tuple[int, int, float]]) -> np.ndarray:
    i, j = np.indices((N, N))
    m = 1.0 / (1.0 + np.abs(i - j)) ** 0.95
    for a, b in zip(TAD_EDGES[:-1], TAD_EDGES[1:]):
        m = np.where((i >= a) & (i < b) & (j >= a) & (j < b), m * 2.4, m)
    for a, b in zip(TAD_EDGES[:-1], TAD_EDGES[1:]):
        for ci, cj in ((a + 2, b - 3), (b - 3, a + 2)):
            m += 0.5 * np.exp(-(((i - ci) ** 2 + (j - cj) ** 2) / 8.0))
    for a, b, amp in loops:
        m += amp * np.exp(-(((i - a) ** 2 + (j - b) ** 2) / 6.0))
        m += amp * np.exp(-(((i - b) ** 2 + (j - a) ** 2) / 6.0))
    m = m * RNG.lognormal(0.0, 0.20, size=m.shape)
    return np.log1p((m + m.T) / 2.0 * 60.0)


BASE = contact_matrix([])
STAPLED = contact_matrix([(E_BIN, P_BIN, 3.2)])
VMAX = float(np.percentile(BASE, 99.2))

# ------------------------------------------------------------------ canvas
# Displayed ~1,488 px wide, so one point of type is ~2.6 px on the slide and
# one data unit is 0.08 in. A node's usable width is its width less 3 units;
# copy is written to fit that rather than trimmed afterwards. Three bands,
# read as a snake: build the map, compute and intervene, score and rank.
W_IN, H_IN = 8.0, 4.18
XMAX = 100.0
YMAX = XMAX * H_IN / W_IN

T_TITLE, T_PATH, T_SUB, T_TAG, T_TINY = 9.4, 8.0, 8.3, 8.0, 7.7
LINE = 1.35


def hic(ax, mat, x0, x1, ybase, ytop, *, boundaries=False):
    """Place a contact matrix as the conventional triangle inside a rect."""
    bw, bh = x1 - x0, ytop - ybase
    tr = (
        mtransforms.Affine2D()
        .translate(-N / 2, -N / 2).rotate_deg(-45).scale(1 / np.sqrt(2))
        .translate(N / 2, 0).scale(bw / N, bh / (N / 2.0)).translate(x0, ybase)
        + ax.transData
    )
    im = ax.imshow(mat, cmap=HIC_CMAP, vmin=0.0, vmax=VMAX, origin="lower",
                   extent=(0, N, 0, N), interpolation="bilinear", zorder=2)
    im.set_transform(tr)
    im.set_clip_path(Rectangle((x0, ybase), bw, bh, transform=ax.transData))
    if boundaries:
        for b in TAD_EDGES[1:-1]:
            bx = x0 + bw * b / N
            ax.plot([bx, bx], [ybase, ybase + bh * 0.58], lw=1.2,
                    ls=(0, (3, 2.5)), color=AMBER, zorder=4)


def node(ax, x0, y0, x1, y1, *, title, path=None, owner=TEAL, subs=(),
         emits=None, fill=BG, number=None):
    ax.add_patch(FancyBboxPatch(
        (x0, y0), x1 - x0, y1 - y0,
        boxstyle="round,pad=0,rounding_size=0.9",
        facecolor=fill, edgecolor=EDGE, lw=1.3, zorder=3))
    ax.add_patch(Rectangle((x0, y0), 0.55, y1 - y0, facecolor=owner,
                           edgecolor="none", zorder=4))
    tx, y = x0 + 1.6, y1 - 1.9
    if number is not None:
        ax.text(tx, y, number, fontsize=T_PATH, color=GREY, family=MONO,
                weight="bold", va="center", zorder=5)
        tx += 2.5
    ax.text(tx, y, title, fontsize=T_TITLE, color=INK, weight="bold",
            va="center", zorder=5)
    y -= 1.5
    if path:
        ax.text(x0 + 1.6, y, path, fontsize=T_PATH, color=owner, family=MONO,
                va="center", zorder=5)
        y -= 1.5
    for line in subs:
        ax.text(x0 + 1.6, y, line, fontsize=T_SUB, color=MUTE, va="center",
                zorder=5)
        y -= LINE
    if emits:
        ax.text(x1 - 1.4, y0 + 1.15, emits, fontsize=T_TAG, color=INK,
                family=MONO, ha="right", va="center", zorder=5)


def arrow(ax, p0, p1, *, color="#b9b2a2", lw=1.6, rad=0.0, dashed=False):
    ax.add_patch(FancyArrowPatch(
        p0, p1, arrowstyle="-|>", mutation_scale=10, lw=lw, color=color,
        connectionstyle=f"arc3,rad={rad}", zorder=6, shrinkA=0, shrinkB=0,
        linestyle=(0, (3, 2)) if dashed else "solid"))


def elbow(ax, pts, *, color="#b9b2a2", lw=1.6, label=None, label_at=None):
    xs, ys = zip(*pts)
    ax.plot(xs, ys, color=color, lw=lw, solid_capstyle="round", zorder=6)
    arrow(ax, pts[-2], pts[-1], color=color, lw=lw)
    if label and label_at:
        ax.text(*label_at, label, fontsize=T_TAG, color=INK, family=MONO,
                ha="center", va="bottom", zorder=7,
                bbox=dict(facecolor=BG, edgecolor="none", pad=1.8))


def band_label(ax, x, y, text):
    ax.text(x, y, text, fontsize=T_TAG, color=GREY, family=MONO,
            weight="bold", va="center")


def build(ax):
    # ============================================ band 1: build one map
    band_label(ax, 1, 51.4, "1  BUILD ONE CONTACT MAP")

    ax.text(1, 49.4, "DNA sequence", fontsize=T_SUB, color=INK,
            weight="bold", va="center")
    for k, b in enumerate(RNG.integers(0, 4, 84)):
        ax.add_patch(Rectangle((1 + k * 14 / 84.0, 47.3), 14 / 84.0, 1.2,
                               color=[TEAL, AMBER, INK, "#c2bcae"][b],
                               zorder=3))

    ax.text(1, 45.9, "ATAC-seq", fontsize=T_SUB, color=INK, weight="bold",
            va="center")
    xg = np.linspace(0, N, 600)
    sig = RNG.lognormal(-1.1, 0.42, 600) * 0.5
    for c, h, sd in [(22, 2.1, 1.5), (38, 4.4, 1.9), (61, 1.5, 1.2),
                     (84, 5.6, 2.1), (103, 2.6, 1.6), (122, 1.8, 1.3),
                     (139, 3.2, 1.7), (151, 6.2, 2.0), (166, 2.0, 1.4)]:
        sig += h * np.exp(-((xg - c) ** 2) / (2 * sd ** 2))
    px, py = 1 + 14 * xg / N, 41.6 + 3.4 * sig / 7.2
    ax.fill_between(px, 41.6, py, color=TEAL, alpha=0.28, lw=0, zorder=3)
    ax.plot(px, py, color=TEAL, lw=1.0, zorder=4)

    ax.text(1, 40.2, "Hi-C", fontsize=T_SUB, color=INK, weight="bold",
            va="center")
    ax.text(5.4, 40.2, "a few cell types", fontsize=T_TINY, color=GREY,
            va="center")
    hic(ax, BASE, 1, 15, 34.2, 39.2)

    node(ax, 20, 43.0, 40, 50.4,
         number="01", title="UniversalEPI", path="pipeline/contact.py",
         owner=TEAL, subs=["ATAC + sequence in", "ensemble of 10"])
    node(ax, 20, 36.4, 40, 41.4,
         title="fallbacks", owner=GREY, fill=PANEL,
         subs=["megamap, power law"])

    arrow(ax, (15.3, 47.9), (19.7, 48.4))
    arrow(ax, (15.3, 43.4), (19.7, 45.4))
    arrow(ax, (15.3, 36.8), (19.7, 38.6))

    ax.plot([15.3, 51.5], [34.8, 34.8], color=TEAL, lw=2.3, zorder=5,
            solid_capstyle="round")
    arrow(ax, (50, 34.8), (52.7, 34.8), color=TEAL, lw=2.3)
    ax.text(33, 35.4, "measured · preferred", fontsize=T_TAG, color=TEAL,
            family=MONO, ha="center", va="bottom")

    arrow(ax, (40.3, 46.4), (52.7, 44.0))
    arrow(ax, (40.3, 38.9), (52.7, 37.5))

    ax.text(53, 51.4, "ONE CONTACT MAP", fontsize=T_SUB, color=INK,
            weight="bold", va="center")
    hic(ax, BASE, 53, 72, 34.2, 50.0, boundaries=True)

    node(ax, 76, 36.0, 99, 50.0,
         title="provenance", owner=AMBER, fill=PANEL,
         subs=["every record names the", "rung behind it:", "",
               "Hi-C         53.3%", "megamap      48.0%",
               "power law    44.3%", "UniversalEPI  n/a"])
    arrow(ax, (72.3, 43.0), (75.7, 43.0))

    elbow(ax, [(62, 34.0), (62, 32.2), (3, 32.2), (3, 30.2)],
          label="locus.contacts", label_at=(34, 32.6))

    # ====================================== band 2: baseline, then intervene
    band_label(ax, 7, 30.9, "2  BASELINE, THEN INTERVENE")

    node(ax, 1, 17.0, 30, 30.0,
         number="02", title="Baseline ABC", path="pipeline/abc.py",
         owner=TEAL,
         subs=["A_e = geomean(DNase, H3K27ac)",
               "ABC(E,G) = A_E×C / Σ A_e×C",
               "T_G = Σ A_e × C(e,G)",
               "",
               "T is the dial, ABC says",
               "which element moved it"],
         emits="baseline")

    node(ax, 33, 17.0, 62, 30.0,
         number="03", title="Anchor search", path="pipeline/anchors.py",
         owner=TEAL,
         subs=["ATAC peaks, minus the element,",
               "TSS ±1 kb, gene body, CTCF",
               "MA0139.1, and ±1 kb around",
               "every other element",
               "",
               "Zero survivors is valid."],
         emits="staples[]")

    node(ax, 65, 17.0, 99, 30.0,
         number="04", title="Intervene, then sweep",
         path="pipeline/sweep.py", owner=TEAL,
         subs=["C' = max( C, q · K(i−A, j−B) )",
               "K Gaussian, σ ≈ 2–3 bins. A point edit",
               "would make hijacking impossible.",
               "Row sums preserved: a locus has a",
               "finite contact budget.",
               "q at p50…p99 of real contacts"],
         emits="Δ sweep")

    arrow(ax, (30.3, 23.5), (32.7, 23.5))
    arrow(ax, (62.3, 23.5), (64.7, 23.5))

    elbow(ax, [(82, 16.8), (82, 15.6), (3, 15.6), (3, 14.2)],
          label="staples[].delta", label_at=(70, 16.0))

    # ============================================ band 3: score and rank
    band_label(ax, 7, 14.9, "3  SCORE, GRADE, RANK")

    node(ax, 1, 1.0, 26, 14.0,
         number="05", title="Re-score window", path="scoring/predict.py",
         owner=AMBER,
         subs=["Recompute ABC and T on C'",
               "for every gene in ±1 Mb",
               "ΔT/T   relative dose",
               "ΔABC   which element",
               "up collateral, down theft"],
         emits="Δ per gene")

    node(ax, 29, 1.0, 54, 14.0,
         number="06", title="Annotate", path="scoring/risk.py", owner=AMBER,
         subs=["COSMIC CGC oncogenes",
               "ClinGen TS3: must not rise",
               "ClinGen HI3: must not fall",
               "DepMap common-essential",
               "Can raise, never lower"],
         emits="tiers")

    gates = [
        ("TRACTABLE", "an anchor survived", TEAL),
        ("CONFIDENT", "rung · ensemble FC > 0", AMBER),
        ("SAFE", "nothing annotated rose", TEAL),
        ("DELIVERABLE", "cell type · reversible", TEAL),
    ]
    centres = []
    for k, (name, test, col) in enumerate(gates):
        top = 14.0 - k * 3.35
        bot = top - 2.9
        ax.add_patch(FancyBboxPatch(
            (57, bot), 22, 2.9, boxstyle="round,pad=0,rounding_size=0.6",
            facecolor=BG, edgecolor=EDGE, lw=1.2, zorder=3))
        ax.add_patch(Rectangle((57, bot), 0.5, 2.9, facecolor=col,
                               edgecolor="none", zorder=4))
        ax.text(58.3, top - 0.95, name, fontsize=T_SUB, color=col,
                family=MONO, weight="bold", va="center", zorder=5)
        ax.text(58.3, top - 2.15, test, fontsize=T_TINY, color=MUTE,
                va="center", zorder=5)
        centres.append((top + bot) / 2.0)
        arrow(ax, (54.3, 7.5), (56.7, (top + bot) / 2.0))

    node(ax, 82, 4.5, 99, 14.0,
         number="07", title="Rank", path="scoring/rank.py", owner=AMBER,
         subs=["1 Safe", "2 Confident", "3 ΔT/T", "4 Deliverable"])
    for c in centres:
        arrow(ax, (79.3, c), (81.7, 9.0), rad=0.05)

    ax.add_patch(FancyBboxPatch(
        (82, 0.5), 17, 3.3, boxstyle="round,pad=0,rounding_size=0.6",
        facecolor=PANEL, edgecolor=ROSE, lw=1.3, zorder=3))
    ax.text(90.5, 2.8, "ABSTAIN", fontsize=T_SUB, color=ROSE, family=MONO,
            weight="bold", ha="center", va="center", zorder=5)
    ax.text(90.5, 1.4, "max-conf FC = 0", fontsize=T_TINY, color=MUTE,
            ha="center", va="center", zorder=5)
    ab = [(79.3, centres[1]), (80.5, centres[1]), (80.5, 2.1), (81.3, 2.1)]
    ax.plot(*zip(*ab), color=ROSE, lw=1.4, ls=(0, (3, 2)), zorder=6,
            solid_capstyle="round")
    arrow(ax, (80.8, 2.1), (81.9, 2.1), color=ROSE, lw=1.4)


def main(out="docs/figures/pipeline_flow.png"):
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


if __name__ == "__main__":
    main()

"""Procedural 2.5D render of storyboard.md (no Blender): numpy + scipy + Pillow, encoded with ffmpeg.

Look (after the Peter Mac "What Goes Wrong in Cancer?" film): layered slate-blue volumetric backdrop
with shallow depth of field, glossy pastel molecular surfaces shaded from height fields (warm key,
cool fill and rim, ambient occlusion, soft blue-tinted cast shadows), bloom, grade and grain, bold
captions with leader lines. Everything is illustrative and not to scale; no structure is taken from PDB.

    pip install -r animation/requirements.txt
    python animation/procedural.py --width 1920 --out animation/renders/procedural_1080p
    python animation/procedural.py --width 1280 --out animation/renders/preview --frames 0,150,600
"""
import argparse
import math
import os
import subprocess
from functools import lru_cache
from multiprocessing import Pool

import numpy as np
from PIL import Image, ImageChops, ImageDraw, ImageFilter, ImageFont
from scipy.ndimage import gaussian_filter

FPS, NFRAMES = 24, 720
BW, BH = 1280, 720            # world units for the TAD scenes == pixels at zoom 1 and 720p
W = H = S = None              # output size and pixels per world unit, set by init()
CACHE = None
FONT_BOLD = "C:/Windows/Fonts/segoeuib.ttf"
FONT_SEMI = "C:/Windows/Fonts/seguisb.ttf"
CENTER = np.array([650.0, 380.0])

PAL = {
    "bg_top": (78, 100, 128), "bg_bot": (34, 46, 66), "haze": (255, 236, 214),
    "dnaA": (244, 248, 255), "dnaB": (212, 196, 246), "nuc": (198, 178, 238),
    "pol2": (255, 158, 200), "gtf": (210, 178, 252), "ctcf": (124, 238, 218), "coh": (104, 214, 238),
    "enh": (255, 188, 136), "rep": (150, 176, 210), "stap": (52, 214, 186), "mrna": (255, 174, 112),
    "glow": (255, 208, 176),
}
BASES = [((255, 226, 132), (140, 206, 255)), ((255, 160, 192), (170, 236, 150))]


def init(width, cache):
    global W, H, S, CACHE
    W, H, S, CACHE = width, round(width * 9 / 16), width / BW, cache


# ---------------------------------------------------------------- helpers
def clamp01(x):
    return max(0.0, min(1.0, x))


def ss(a, b, x):
    t = clamp01((x - a) / (b - a))
    return t * t * (3 - 2 * t)


def lerp(a, b, t):
    return a + (b - a) * t


def unit(v):
    v = np.asarray(v, float)
    n = math.hypot(v[0], v[1])
    return v / n if n > 1e-6 else np.array([1.0, 0.0])


def lerp_angle(a, b, t):
    return a + ((b - a + math.pi) % (2 * math.pi) - math.pi) * t


def catmull(P, per=24):
    P = np.asarray(P, float)
    ext = np.vstack([2 * P[0] - P[1], P, 2 * P[-1] - P[-2]])
    p0, p1, p2, p3 = (ext[i:len(ext) - 3 + i][:, None, :] for i in range(4))
    t = np.linspace(0, 1, per, endpoint=False)[None, :, None]
    out = 0.5 * (2 * p1 + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t ** 2
                 + (-p0 + 3 * p1 - 3 * p2 + p3) * t ** 3)
    return np.vstack([out.reshape(-1, 2), P[-1:]])


def at_u(curve, u, per=24):
    f = u * per
    i = int(min(max(math.floor(f), 0), len(curve) - 2))
    return lerp(curve[i], curve[i + 1], f - i)


class View:
    """World -> pixel mapping for a canvas (full res: sf=S, half res: sf=S/2)."""

    def __init__(self, cx, cy, z, sf, w, h):
        self.cx, self.cy, self.z, self.sf, self.w, self.h = cx, cy, z, sf, w, h
        self.k = z * sf

    def pts(self, p):
        p = np.asarray(p, float)
        return np.column_stack([(p[:, 0] - self.cx) * self.k + self.w / 2, (p[:, 1] - self.cy) * self.k + self.h / 2])

    def pt1(self, p):
        return (p[0] - self.cx) * self.k + self.w / 2, (p[1] - self.cy) * self.k + self.h / 2


# ---------------------------------------------------------------- lighting and sprites
def _u3(*v):
    v = np.array(v, np.float32)
    return v / np.linalg.norm(v)


KEY, FILL = _u3(-0.5, -0.65, 0.75), _u3(0.65, 0.45, 0.45)
HALF = _u3(*(KEY + np.array([0, 0, 1], np.float32)))
KEYC = np.array([1.0, 0.95, 0.88], np.float32)
FILLC = np.array([0.70, 0.82, 1.0], np.float32)
RIMC = np.array([0.80, 0.93, 1.0], np.float32)


def shade(nx, ny, nz, col, ao, rim_k=0.42, spec_k=0.5):
    """Warm wrapped key, cool fill, Blinn spec and a cool fresnel rim."""
    kd = np.clip((nx * KEY[0] + ny * KEY[1] + nz * KEY[2] + 0.35) / 1.35, 0, 1)[..., None]
    fd = np.clip(nx * FILL[0] + ny * FILL[1] + nz * FILL[2], 0, 1)[..., None]
    sp = (np.clip(nx * HALF[0] + ny * HALF[1] + nz * HALF[2], 0, 1) ** 48)[..., None]
    rim = (np.clip(1 - nz, 0, 1) ** 2.2)[..., None]
    ao = np.asarray(ao, np.float32)[..., None]
    light = 0.30 * ao + kd * 0.80 * KEYC * (0.55 + 0.45 * ao) + fd * 0.28 * FILLC
    return col * light + sp * spec_k * ao + rim * rim_k * RIMC * (0.6 + 0.4 * col)


def to_rgba(rgb, a):
    return Image.fromarray(np.dstack([np.clip(rgb, 0, 1) * 255, np.clip(a, 0, 1) * 255]).astype(np.uint8), "RGBA")


@lru_cache(maxsize=None)
def sphere(q, rgb, k=1.0):
    """Shaded sphere sprite of diameter q px."""
    r = max(q / 2, 0.6)
    n = int(math.ceil(2 * r)) + 2
    c = n / 2
    y, x = np.mgrid[0:n, 0:n].astype(np.float32) + 0.5
    dx, dy = (x - c) / r, (y - c) / r
    d = np.sqrt(dx * dx + dy * dy)
    dd = np.maximum(d, 1.0)
    nx, ny = dx / dd, dy / dd
    nz = np.sqrt(np.clip(1 - nx * nx - ny * ny, 0, 1))
    col = np.array(rgb, np.float32) / 255 * k
    return to_rgba(shade(nx, ny, nz, col, 0.7 + 0.3 * nz), (1 - d) * r + 0.5)


def build_surface(spheres, ps, smooth=0.3):
    """Molecular-surface sprite: union of spheres as a height field, smoothed normals, crevice AO."""
    ext = max(max(abs(x), abs(y)) + r for x, y, z, r, c in spheres) + 2
    n = int(2 * ext * ps) + 8
    c0 = n / 2
    Hf = np.full((n, n), -1e9, np.float32)
    C = np.zeros((n, n, 3), np.float32)
    A = np.zeros((n, n), np.float32)
    rs = []
    for x, y, z, r, rgb in spheres:
        R = r * ps
        px, py = c0 + x * ps, c0 + y * ps
        x0, x1 = max(int(px - R - 2), 0), min(int(px + R + 3), n)
        y0, y1 = max(int(py - R - 2), 0), min(int(py + R + 3), n)
        yy, xx = np.mgrid[y0:y1, x0:x1].astype(np.float32) + 0.5
        d2 = (xx - px) ** 2 + (yy - py) ** 2
        h = z * ps + np.sqrt(np.clip(R * R - d2, 0, None))
        sub = Hf[y0:y1, x0:x1]
        upd = (d2 < R * R) & (h > sub)
        sub[upd] = h[upd]
        C[y0:y1, x0:x1][upd] = np.array(rgb, np.float32) / 255
        A[y0:y1, x0:x1] = np.maximum(A[y0:y1, x0:x1], np.clip(R - np.sqrt(d2) + 0.5, 0, 1))
        rs.append(R)
    fg = Hf > -1e8
    Hf[~fg] = Hf[fg].min() - 0.25 * ext * ps
    C[~fg] = C[fg].mean(axis=0)
    mr = float(np.mean(rs))
    gy, gx = np.gradient(gaussian_filter(Hf, smooth * mr))
    nrm = np.sqrt(gx * gx + gy * gy + 1)
    ao = np.clip(1 - np.clip(gaussian_filter(Hf, 0.22 * ext * ps) - Hf, 0, None) / (0.30 * ext * ps), 0.3, 1)
    C = gaussian_filter(C, (0.2 * mr, 0.2 * mr, 0))
    return to_rgba(shade(-gx / nrm, -gy / nrm, 1 / nrm, C, ao), A)


def blob_spheres(rgb, size, n, seed, shape=(1.0, 1.0), cleft=False, accent=None, jitter=0.05):
    rng = np.random.default_rng(seed)
    rx, ry = size / 2 * shape[0], size / 2 * shape[1]
    out = []
    while len(out) < n:
        x, y = rng.uniform(-1, 1, 2)
        rr = x * x + y * y
        if rr > 1 or (cleft and y < -0.05 and abs(x) < 0.22 - y * 0.25):
            continue
        rad = size * rng.uniform(0.07, 0.13)
        z = math.sqrt(1 - rr) * min(rx, ry) * 0.75 + rng.uniform(-0.08, 0.08) * size
        col = tuple(np.clip(np.array(rgb, float) * (1 + rng.normal(0, jitter)), 0, 255))
        if accent is not None and rng.random() < 0.12:
            col = accent
        out.append((x * (rx - rad * 0.6), y * (ry - rad * 0.6), z, rad, col))
    return out


def ring_spheres(rgb, R, ry, tube, front, seed=7):
    rng = np.random.default_rng(seed)
    out = []
    for a in np.linspace(0, 2 * math.pi, 150, endpoint=False):
        sa = math.sin(a)
        if (sa > -0.05) != front:
            continue
        for _ in range(2):
            out.append((math.cos(a) * R + rng.normal(0, tube * 0.25), sa * ry + rng.normal(0, tube * 0.25),
                        sa * ry * 1.2 + rng.normal(0, 1), tube * rng.uniform(0.55, 0.95),
                        tuple(np.clip(np.array(rgb, float) * (1 + rng.normal(0, 0.05)), 0, 255))))
    return out


def grna_spheres(size, z0):
    return [(math.cos(a) * size * 0.42, math.sin(a) * size * 0.42 - 2, z0, 3.0, PAL["mrna"])
            for a in np.linspace(-1.1, 0.7, 16)]


def nucleosome_spheres():
    sp = blob_spheres(PAL["nuc"], 24, 45, 8, (1.0, 0.8))
    for a in np.linspace(0, 2 * math.pi, 34, endpoint=False):
        sp.append((math.cos(a) * 13.5, math.sin(a) * 10, 4 + 6 * math.sin(a), 2.8, PAL["dnaA"]))
    return sp


@lru_cache(maxsize=None)
def prot(name):
    P = PAL
    sp = {
        "pol2": lambda: blob_spheres(P["pol2"], 80, 230, 1, (1.0, 0.85), cleft=True, accent=(255, 200, 224)),
        "gtf": lambda: blob_spheres(P["gtf"], 52, 120, 2, accent=(232, 214, 255)),
        "ctcf": lambda: blob_spheres(P["ctcf"], 36, 70, 3, (0.8, 1.2)),
        "enh": lambda: blob_spheres(P["enh"], 58, 140, 4, accent=(255, 220, 180)),
        "rep": lambda: blob_spheres(P["rep"], 68, 170, 5, (1.0, 0.9), accent=(190, 206, 230)),
        "stapA": lambda: blob_spheres(P["stap"], 46, 110, 6, (1.1, 0.9)) + grna_spheres(46, 14),
        "stapB": lambda: blob_spheres(P["stap"], 46, 110, 9, (1.1, 0.9)) + grna_spheres(46, 14),
        "bridge": lambda: blob_spheres(P["stap"], 24, 30, 10),
        "nuc": nucleosome_spheres,
        "coh_b": lambda: ring_spheres(P["coh"], 48, 17, 7.5, False),
        "coh_f": lambda: ring_spheres(P["coh"], 48, 17, 7.5, True),
    }[name]()
    return build_surface(sp, S * 1.25)


SHADOW_R = {"pol2": 34, "gtf": 22, "ctcf": 15, "enh": 25, "rep": 29, "stapA": 21, "stapB": 21,
            "bridge": 10, "nuc": 12, "coh_f": 30}


@lru_cache(maxsize=None)
def prot_scaled(name, zq, aq=0):
    base = prot(name)
    f = math.exp(zq / 40) / (S * 1.25)
    im = base.convert("RGBa").resize((max(2, round(base.width * f)), max(2, round(base.height * f))), Image.LANCZOS)
    if aq:
        im = im.rotate(aq * 10, Image.BICUBIC, expand=True)
    return im.convert("RGBA")


_faded = {}


def blit(cv, spr, x, y, alpha=1.0):
    if alpha <= 0.02:
        return
    if alpha < 0.98:
        k = max(1, int(alpha * 32))
        if spr.width * spr.height <= 90000:
            key = (id(spr), k)
            hit = _faded.get(key)
            if hit is None or hit[0] is not spr:
                f = spr.copy()
                f.putalpha(spr.getchannel("A").point(lambda v: v * k // 32))
                hit = _faded[key] = (spr, f)
            spr = hit[1]
        else:
            f = spr.copy()
            f.putalpha(spr.getchannel("A").point(lambda v: v * k // 32))
            spr = f
    cv.paste(spr, (int(round(x - spr.width / 2)), int(round(y - spr.height / 2))), spr)


def draw_prot(cv, view, name, wp, alpha=1.0, angle=0.0):
    x, y = view.pt1(wp)
    m = 300 * view.sf
    if -m < x < view.w + m and -m < y < view.h + m:
        blit(cv, prot_scaled(name, int(round(40 * math.log(view.k))), int(round(angle / 10)) % 36), x, y, alpha)


@lru_cache(maxsize=32)
def soft_disc(rq, rgb, blur):
    """Blurred disc sprite; big ones are built small and upscaled (they are blurry anyway)."""
    n = int(2 * rq + 4 * blur + 4)
    ds = max(1, n // 160)
    m = n // ds
    r, b = rq / ds, blur / ds
    a = Image.new("L", (m, m), 0)
    ImageDraw.Draw(a).ellipse([m / 2 - r, m / 2 - r, m / 2 + r, m / 2 + r], fill=255)
    a = a.filter(ImageFilter.GaussianBlur(b))
    if ds > 1:
        a = a.resize((n, n), Image.BILINEAR)
    im = Image.new("RGBA", a.size, rgb + (0,))
    im.putalpha(a)
    return im


@lru_cache(maxsize=8)
def glow_mask(soft):
    m = 192
    r = m / 2 / (1 + 2 * soft)
    a = Image.new("L", (m, m), 0)
    ImageDraw.Draw(a).ellipse([m / 2 - r, m / 2 - r, m / 2 + r, m / 2 + r], fill=255)
    return a.filter(ImageFilter.GaussianBlur(r * soft))


def draw_soft(cv, x, y, r, col, alpha, soft=0.5):
    """Large soft disc: only the visible part of a small mask is resized, so huge glows stay cheap."""
    if alpha <= 0.01 or r < 1:
        return
    n = 2 * r * (1 + 2 * soft)
    X0, Y0 = x - n / 2, y - n / 2
    vx0, vy0 = max(0, math.ceil(X0)), max(0, math.ceil(Y0))
    vx1, vy1 = min(cv.width, math.floor(X0 + n)), min(cv.height, math.floor(Y0 + n))
    if vx1 - vx0 < 1 or vy1 - vy0 < 1:
        return
    gm = glow_mask(soft)
    s = gm.width / n
    m = gm.resize((vx1 - vx0, vy1 - vy0), Image.BILINEAR,
                  box=((vx0 - X0) * s, (vy0 - Y0) * s, (vx1 - X0) * s, (vy1 - Y0) * s))
    a = clamp01(alpha)
    cv.paste(col, (vx0, vy0, vx1, vy1), m.point(lambda v: int(v * a)))


def draw_glow(cv, view, wp, radius, intensity, col):
    x, y = view.pt1(wp)
    draw_soft(cv, x, y, radius * view.k, col, intensity)


# ---------------------------------------------------------------- DNA
def draw_dna(cv, view, curve, t, colA, colB, bright=1.0, alpha=1.0, phase=0.0):
    scr = view.pts(curve)
    m = 60 * view.sf
    if scr[:, 0].max() < -m or scr[:, 0].min() > view.w + m or scr[:, 1].max() < -m or scr[:, 1].min() > view.h + m:
        return
    k = view.k
    r = 3.3 * k
    helix = r >= 2.2
    cumw = np.r_[0, np.cumsum(np.hypot(*np.diff(curve, axis=0).T))]
    sw = 3.8 if helix else max(1.4 / k, 3.3 * 1.5)       # bead spacing in world units, stable across zoom
    s = np.arange(0, cumw[-1], sw)
    x, y = np.interp(s, cumw, scr[:, 0]), np.interp(s, cumw, scr[:, 1])
    vis = (x > -24) & (x < view.w + 24) & (y > -24) & (y < view.h + 24)
    if not vis.any():
        return
    cA = tuple(min(255, int(v * bright)) for v in colA)
    cB = tuple(min(255, int(v * bright)) for v in colB)
    if not helix:
        spr = sphere(max(3, int(round(r * 4.4))), tuple((a + b) // 2 for a, b in zip(cA, cB)), 1.0)
        for xi, yi in zip(x[vis], y[vis]):
            blit(cv, spr, xi, yi, alpha)
        return
    tx, ty = np.gradient(x), np.gradient(y)
    n = np.hypot(tx, ty) + 1e-6
    nx, ny = -ty / n, tx / n
    th = s * 0.33 + t * 1.5 + phase
    amp = 3.8 * k
    sn, cs = np.sin(th), np.cos(th)
    q = max(3, int(round(r * 2)))
    idx = np.nonzero(vis)[0]
    ax, ay, bx, by = x + nx * amp * sn, y + ny * amp * sn, x - nx * amp * sn, y - ny * amp * sn
    back, front = [], []
    for i in idx:
        (front if cs[i] >= 0 else back).append((ax[i], ay[i], cA))
        (front if cs[i] < 0 else back).append((bx[i], by[i], cB))
    for px, py, col in back:
        blit(cv, sphere(q, col, 0.66), px, py, alpha)
    qr = max(2, int(round(r * 1.3)))
    for i in idx[::2]:
        if abs(sn[i]) > 0.4:
            pair = BASES[((i * 2654435761) >> 9) % 2]
            flip = ((i * 40503) >> 6) % 2 == 0
            for f, col in ((0.33, pair[flip]), (0.67, pair[not flip])):
                blit(cv, sphere(qr, tuple(min(255, int(v * bright)) for v in col), 0.95),
                     lerp(ax[i], bx[i], f), lerp(ay[i], by[i], f), alpha)
    for px, py, col in front:
        blit(cv, sphere(q, col, 1.0), px, py, alpha)


# ---------------------------------------------------------------- TAD and its timeline
CTRL = np.array([
    (-700, 1000), (-200, 830), (250, 700), (600, 615), (470, 530), (340, 400), (360, 265), (470, 165),
    (640, 125), (800, 165), (905, 265), (930, 400), (845, 535), (680, 615), (1030, 700), (1480, 830), (1980, 1000),
], float)
E, R, P = 5, 9, 11
GENE_U = (11.05, 12.3)
JIG_TAD = np.ones(len(CTRL))
JIG_TAD[[0, 1, 15, 16]] = 0
NUC_U = (1.6, 2.2, 3.9, 4.35, 6.0, 6.5, 7.0, 7.6, 8.25, 10.1, 13.4, 13.9, 14.4, 14.9, 15.5)
E_WTS = {4: .5, 5: 1, 6: .55, 7: .2}
R_WTS = {8: .4, 9: 1, 10: .3}


def enh_nat(t):
    return ss(12.0, 13.3, t) - ss(16.0, 17.0, t)


def enh_staple(t):
    return ss(23.0, 25.3, t)


def rep_contact(t):
    return ss(16.6, 17.8, t) - ss(19.8, 21.0, t)


def tad_ctrl(t):
    p = CTRL.copy()
    for amt, tgt, idx, wts in ((enh_nat(t), (835, 400), E, E_WTS), (enh_staple(t), (856, 394), E, E_WTS),
                               (rep_contact(t), (850, 318), R, R_WTS)):
        if amt > 0:
            d = (np.array(tgt) - CTRL[idx]) * amt
            for k, w in wts.items():
                p[k] += d * w
    i = np.arange(len(p))
    p[:, 0] += 3.5 * JIG_TAD * np.sin(t * 1.7 + i * 1.3)
    p[:, 1] += 3.5 * JIG_TAD * np.sin(t * 1.3 + i * 2.1)
    return p


def rate(t):
    """Pol II initiations per second (illustrative): baseline, enhancer up, repressor down, staple up."""
    if t < 7.2:
        return 0.0
    r = 0.4 + 0.4 * enh_nat(t) + 0.5 * ss(24.8, 26.2, t)
    return r * (1 - 0.68 * rep_contact(t))


def schedule():
    spawns, acc, rs, sm = [], 0.95, [], 0.0
    dt = 1 / 240
    for k in range(int(31 / dt)):
        t = k * dt
        r = rate(t)
        acc += r * dt
        if acc >= 1:
            acc -= 1
            spawns.append(t)
        sm += (r - sm) * dt / 1.2
        rs.append(sm)
    return spawns, np.array(rs)


SPAWNS, RATE_S = schedule()


def rate_s(t):
    return float(RATE_S[min(int(t * 240), len(RATE_S) - 1)])


def anchor_cam(wx, wy, sx, sy, z):
    """Camera that puts world point (wx, wy) at 720p-screen point (sx, sy)."""
    return wx - (sx - 640) / z, wy - (sy - 360) / z, z


def cam_s35(t):
    z = lerp(0.94, 0.98, ss(6, 20, t))
    cx, cy, _ = anchor_cam(640, 370, 740, 420, z)
    k = ss(6, 8, t)
    return cx + 8 * math.sin(t * 0.3) * k, cy + 5 * math.sin(t * 0.23) * k, z


def cam(t):
    if t < 6:
        k = ss(2.2, 6.0, t)
        end = cam_s35(6.0)
        return lerp(640, end[0], k), lerp(380, end[1], k), math.exp(lerp(math.log(0.2), math.log(end[2]), k))
    if t < 20:
        return cam_s35(t)
    pull = anchor_cam(640, 370, 790, 440, 0.86)
    push = anchor_cam(720, 380, 820, 450, 0.98)
    if t < 23:
        a, b, k = cam_s35(20.0), pull, ss(20, 22.6, t)
    elif t < 29:
        a, b, k = pull, push, ss(23, 28.6, t)
    else:
        a, b, k = push, anchor_cam(700, 380, 800, 440, 0.92), ss(29, 30, t)
    return tuple(lerp(u, v, k) for u, v in zip(a, b))


def stapler(t, enh, prom):
    """Fused dCas9-dCas9: half A docks beside the enhancer, both drift to the promoter, half B locks."""
    if t < 19.8:
        return None
    d = 36.0
    enh_site, prom_site = enh + np.array([-6, -34]), prom + np.array([-44, -40])
    hover = np.array([lerp(1150, 700, ss(19.8, 22.0, t)), lerp(60, 260, ss(19.8, 22.0, t))])
    tumble = 2.0 + t * 0.8
    uh = np.array([math.cos(tumble), math.sin(tumble)])
    k1 = ss(22.0, 23.0, t)
    A = lerp(hover - uh * d / 2, enh_site, k1)
    ut = unit(prom_site - enh_site)
    ang = lerp_angle(tumble, math.atan2(ut[1], ut[0]), k1)
    B = A + np.array([math.cos(ang), math.sin(ang)]) * d
    B = lerp(B, prom_site, ss(25.3, 25.9, t))
    return A, B, ss(19.8, 20.8, t), enh_site, prom_site


# ---------------------------------------------------------------- chromatin around the TAD: one continuous fibre
TAD_LOOP = np.array([(-40, 0), (-170, -85), (-300, -215), (-280, -350), (-170, -450), (0, -490), (160, -450),
                     (265, -350), (290, -215), (205, -80), (40, 0)], float)


def build_chains():
    rng = np.random.default_rng(42)
    rows = [-1450, -550, 350, 1250, 2150]
    seq, bases, jig = [np.array([-5200.0, -1450.0])], [], [1.0]
    chains = []
    prev = seq[0]
    for ri, ry in enumerate(rows):
        row = []
        for cx in range(-3900, 5101, 1000):
            p = np.array([cx + rng.uniform(-220, 220), ry + rng.uniform(-200, 200)])
            if not (-500 < p[0] < 1800 and -350 < p[1] < 1150):
                row.append(p)
        if ri % 2:
            row = row[::-1]
        if ri == 2:
            k = sum(1 for p in row if p[0] < 640)
            row = row[:k] + ["TAD"] + row[k:]
        for p in row:
            if isinstance(p, str):
                seq += [CTRL[0], CTRL[1]]
                jig += [0.0, 0.0]
                chains.append((np.array(seq), np.array(jig), bases))
                seq, jig, bases = [CTRL[15], CTRL[16]], [0.0, 0.0], []
                prev = CTRL[16]
                continue
            ang = rng.uniform(0, 2 * math.pi)
            sc = rng.uniform(0.55, 0.9)
            u = np.array([math.cos(ang), math.sin(ang)])
            e = np.array([-math.sin(ang), math.cos(ang)])
            B = p - u * 245 * sc
            to_w = lambda lx, ly: B + lx * e - ly * u
            entry, exit_ = to_w(-160 * sc, 120 * sc), to_w(160 * sc, 120 * sc)
            seq.append((prev + entry) / 2 + rng.normal(0, 120, 2))
            jig.append(1.0)
            seq.append(entry)
            jig.append(1.0)
            iL = len(seq)
            for j, (lx, ly) in enumerate(TAD_LOOP):
                f = 1.0 if j in (0, len(TAD_LOOP) - 1) else sc * (1 + rng.normal(0, 0.08))
                seq.append(to_w(lx * f, ly * f))
                jig.append(1.0)
            bases.append((iL, len(seq) - 1, math.degrees(math.atan2(-e[1], e[0]))))
            seq.append(exit_)
            jig.append(1.0)
            prev = exit_
    seq.append(np.array([6200.0, 2150.0]))
    jig.append(1.0)
    chains.append((np.array(seq), np.array(jig), bases))
    return chains


CHAINS = build_chains()


def draw_neighbors(cv, view, t, dim):
    cA = tuple(int(v * dim) for v in PAL["dnaA"])
    cB = tuple(int(v * dim) for v in PAL["dnaB"])
    for ctrl, jig, bases in CHAINS:
        i = np.arange(len(ctrl))
        p = ctrl + 4 * jig[:, None] * np.column_stack([np.sin(t * 1.5 + i * 1.3), np.cos(t * 1.2 + i * 2.3)])
        on = [b for b in bases if -400 * view.sf < view.pt1((p[b[0]] + p[b[1]]) / 2)[0] < view.w + 400 * view.sf]
        mids = [((p[a] + p[b]) / 2, ang) for a, b, ang in on]
        for mid, ang in mids:
            draw_prot(cv, view, "coh_b", mid, dim, ang)
        draw_dna(cv, view, catmull(p), t, cA, cB, phase=1.7)
        for (a, b, ang), (mid, _) in zip(on, mids):
            draw_prot(cv, view, "ctcf", p[a + 1] * 0.75 + p[a] * 0.25, dim)
            draw_prot(cv, view, "ctcf", p[b - 1] * 0.75 + p[b] * 0.25, dim)
            draw_prot(cv, view, "coh_f", mid, dim, ang)


# ---------------------------------------------------------------- background plates (half resolution)
def build_far(path):
    rng = np.random.default_rng(3)
    lw, lh = W // 4, H // 4                   # built small, upscaled to 2x the half-res frame
    y, x = np.mgrid[0:lh, 0:lw].astype(np.float32)
    g = (y / lh)[..., None]
    img = np.array(PAL["bg_top"], np.float32) * (1 - g) + np.array(PAL["bg_bot"], np.float32) * g
    warm = np.exp(-(((x - 0.3 * lw) / (0.45 * lw)) ** 2 + ((y - 0.2 * lh) / (0.45 * lh)) ** 2))[..., None]
    cool = np.exp(-(((x - 0.8 * lw) / (0.5 * lw)) ** 2 + ((y - 0.85 * lh) / (0.4 * lh)) ** 2))[..., None]
    img = img + warm * np.array([62, 56, 48]) + cool * np.array([10, 30, 44])
    im = Image.fromarray(np.clip(img, 0, 255).astype(np.uint8), "RGB")
    layer = Image.new("RGB", (lw, lh), (0, 0, 0))
    d = ImageDraw.Draw(layer)
    for _ in range(26):
        pts = catmull(np.cumsum(rng.normal(0, 0.09 * lw, (6, 2)), axis=0) + rng.uniform([0, 0], [lw, lh]))
        c = int(rng.uniform(14, 30))
        d.line([tuple(p) for p in pts], fill=(c, c + 4, c + 10), width=int(rng.uniform(0.012, 0.03) * lw))
    for k in range(5):                         # light shafts from the upper left
        x0 = lw * (0.05 + 0.17 * k + rng.uniform(-0.03, 0.03))
        wdt = lw * rng.uniform(0.03, 0.07)
        d.polygon([(x0, -10), (x0 + wdt, -10), (x0 + wdt + 0.55 * lw, lh + 10), (x0 + 0.55 * lw, lh + 10)],
                  fill=(26, 24, 20))
    layer = layer.filter(ImageFilter.GaussianBlur(0.025 * lw))
    im = ImageChops.add(im, layer)
    im.resize((W, H), Image.BICUBIC).filter(ImageFilter.GaussianBlur(4)).save(path)


def build_mid(path):
    rng = np.random.default_rng(4)
    sf = S / 2
    im = Image.new("RGBA", (W, H), (170, 176, 210, 0))
    cols = [PAL["dnaB"], PAL["dnaA"], (176, 214, 226), (226, 190, 214), (190, 200, 236)]
    for _ in range(40):
        pts = catmull(np.cumsum(rng.normal(0, 150 * sf, (7, 2)), axis=0) + rng.uniform([0, 0], [W, H]))
        col = tuple(int(v * 0.85) for v in cols[rng.integers(len(cols))])
        cum = np.r_[0, np.cumsum(np.hypot(*np.diff(pts, axis=0).T))]
        rad = rng.uniform(5, 10) * sf
        spr = sphere(int(rad * 2), col, 0.9)
        for s in np.arange(0, cum[-1], rad * 1.1):
            x, y = np.interp(s, cum, pts[:, 0]) - rad, np.interp(s, cum, pts[:, 1]) - rad
            if 0 <= x < W - 2 * rad - 2 and 0 <= y < H - 2 * rad - 2:
                im.alpha_composite(spr, (int(x), int(y)))
    for _ in range(12):
        spr = prot_scaled(["pol2", "rep", "gtf", "enh"][rng.integers(4)], int(round(40 * math.log(sf * rng.uniform(0.6, 1.3)))))
        x, y = rng.uniform(0, W - spr.width), rng.uniform(0, H - spr.height)
        im.alpha_composite(spr, (int(x), int(y)))
    im = im.filter(ImageFilter.GaussianBlur(7 * sf))
    a = im.getchannel("A").point(lambda v: int(v * 0.45))
    im.putalpha(a)
    im.save(path)


@lru_cache(maxsize=None)
def plate(name):
    return Image.open(os.path.join(CACHE, name)).copy()


def crop_plate(im, view, parallax, zexp, w, h):
    sc = view.z ** zexp
    bw, bh = w / sc, h / sc
    cx = im.width / 2 - (view.cx - 640) * view.k * parallax
    cy = im.height / 2 - (view.cy - 370) * view.k * parallax
    cx = min(max(cx, bw / 2), im.width - bw / 2)
    cy = min(max(cy, bh / 2), im.height - bh / 2)
    return im.resize((w, h), Image.BILINEAR, box=(cx - bw / 2, cy - bh / 2, cx + bw / 2, cy + bh / 2))


DUST = np.random.default_rng(12).uniform(0, 1, (60, 4))


def draw_dust(cv, view, t):
    for fx, fy, dep, ph in DUST:
        par = 0.3 + 0.7 * dep
        x = (fx * cv.width - (view.cx - 640) * view.k * par + 12 * view.sf * math.sin(t * 0.3 + ph * 9)) % cv.width
        y = (fy * cv.height - (view.cy - 370) * view.k * par + 9 * view.sf * math.cos(t * 0.25 + ph * 7)) % cv.height
        r = (1.2 + 2.5 * dep) * view.sf * 2
        blit(cv, soft_disc(max(2, int(r)), (240, 246, 255), max(1, int(r * (1.2 - dep)))), x, y, 0.12 + 0.22 * dep)


BOKEH = [((0.10, 0.90), 150, PAL["dnaB"], 0.6), ((0.95, 0.10), 120, PAL["ctcf"], 1.7),
         ((0.88, 0.96), 190, PAL["pol2"], 2.9), ((0.03, 0.30), 95, PAL["enh"], 4.1)]
FG_MOLS = [("rep", (1.0, 0.42), 3.4, 1.1), ("pol2", (0.02, 0.72), 3.0, 2.3)]


@lru_cache(maxsize=16)
def fg_mol(name, sq):
    base = prot(name)
    size = int(math.exp(sq / 25))
    small = base.convert("RGBa").resize((max(2, size // 6), max(2, size // 6)), Image.BILINEAR)
    small = small.filter(ImageFilter.GaussianBlur(max(1, size // 70)))
    return small.resize((size, size), Image.BILINEAR).convert("RGBA")


def draw_foreground(cv, view, t):
    for (fx, fy), r, col, ph in BOKEH:
        x = fx * W - (view.cx - 640) * view.k * 1.5 + 25 * S * math.sin(t * 0.4 + ph)
        y = fy * H - (view.cy - 370) * view.k * 1.5 + 18 * S * math.cos(t * 0.33 + ph)
        draw_soft(cv, x, y, r * S, col, 0.2, 0.25)
    for name, (fx, fy), sc, ph in FG_MOLS:
        size = 80 * S * sc * view.z
        x = fx * W - (view.cx - 640) * view.k * 1.8 + 30 * S * math.sin(t * 0.21 + ph)
        y = fy * H - (view.cy - 370) * view.k * 1.8 + 20 * S * math.cos(t * 0.17 + ph)
        blit(cv, fg_mol(name, int(round(25 * math.log(size)))), x, y, 0.4)


def cast_shadow(cv, view, curve, items, off, blur, strength):
    """Soft, blue-tinted shadow of the given silhouettes, offset away from the key light."""
    q = 4
    m = Image.new("L", (W // q, H // q), 0)
    d = ImageDraw.Draw(m)
    k = view.k / q
    ox, oy = off[0] * k, off[1] * k
    if curve is not None:
        p = view.pts(curve[::3]) / q
        d.line([(x + ox, y + oy) for x, y in p], fill=200, width=max(1, int(9 * k)))
    for name, pos, a, *_ in items:
        x, y = view.pt1(pos)
        r = SHADOW_R.get(name, 15) * k
        x, y = x / q + ox, y / q + oy
        d.ellipse([x - r, y - r, x + r, y + r], fill=int(255 * a))
    m = m.filter(ImageFilter.GaussianBlur(max(1, blur * k)))
    mul = Image.merge("RGB", [m.point(lambda v, w=w: 255 - int(v * strength * w)) for w in (1.0, 0.93, 0.78)])
    return ImageChops.multiply(cv, mul.resize((W, H), Image.BILINEAR))


# ---------------------------------------------------------------- the TAD scenes
def tad_scene(t):
    cx, cy, z = cam(t)
    vf = View(cx, cy, z, S, W, H)
    vh = View(cx, cy, z, S / 2, W // 2, H // 2)
    hl = ss(4.0, 5.3, t)

    bg = crop_plate(plate("far.png"), vh, 0.08, 0.08, vh.w, vh.h)
    mid = crop_plate(plate("mid.png"), vh, 0.3, 0.3, vh.w, vh.h)
    bg.paste(mid, (0, 0), mid)
    draw_neighbors(bg, vh, t, 1 - 0.42 * hl)
    focus = ss(3.6, 6.0, t)
    bg = bg.filter(ImageFilter.GaussianBlur((0.2 + 2.6 * focus) * S))
    draw_dust(bg, vh, t)
    cv = bg.resize((W, H), Image.BILINEAR)

    g = 0.06 + 0.10 * hl + 0.16 * min(1.0, rate_s(t) / 0.9) + 0.10 * ss(28.6, 30, t)
    draw_glow(cv, vf, CENTER, 470, g, PAL["glow"])

    ctrl = tad_ctrl(t)
    curve = catmull(ctrl)
    pt = lambda u: at_u(curve, u)
    prom, enh, sil = pt(P), pt(E), pt(R)

    items = []
    for u in NUC_U:
        p0, p1 = pt(u), pt(u + 0.05)
        items.append(("nuc", p0 + unit([-(p1 - p0)[1], (p1 - p0)[0]]) * 5, 1.0))
    items.append(("ctcf", pt(3.25), 1.0))
    items.append(("ctcf", pt(12.75), 1.0))
    coh_pos = (pt(3) + pt(13)) / 2 + np.array([0, 4])
    items.append(("coh_f", coh_pos, 1.0))
    items.append(("enh", enh + unit(prom - enh) * 18, 1.0))
    items.append(("rep", sil + unit(prom - sil) * 22, 1.0))
    items.append(("gtf", prom + np.array([-18, 4]) + np.array([60, -90]) * (1 - ss(6.6, 7.2, t)), ss(6.5, 7.0, t)))

    glows, chains, pols = [], [], []
    st = stapler(t, enh, prom)
    if st:
        A, B, sa, enh_site, prom_site = st
        items += [("stapA", A, sa), ("bridge", (A + B) / 2 + np.array([0, -2]), sa), ("stapB", B, sa)]
        for (a, b), p in (((22.75, 23.45), enh_site), ((25.55, 26.35), prom_site)):
            if a < t < b:
                glows.append((p, 60, 0.55 * math.sin(math.pi * (t - a) / (b - a)), (170, 255, 236)))

    for t0 in SPAWNS:
        tt = t - t0
        if tt < 0 or tt > 8.6:
            continue
        seed = t0 * 7.3
        wait = prom + np.array([10, -8])
        end = pt(GENE_U[1])
        if tt < 0.7:
            pos, pa = lerp(prom + np.array([90, -110]), wait, ss(0, 0.5, tt)), ss(0, 0.4, tt)
        elif tt < 3.0:
            pos = pt(lerp(GENE_U[0], GENE_U[1], (tt - 0.7) / 2.3)) + np.array([10, -8]) * (1 - ss(0.7, 1.1, tt))
            pa = 1.0
        else:
            pos, pa = end + np.array([75, 10]) * ss(3.0, 3.9, tt), 1 - ss(3.0, 3.9, tt)
        if 0.55 < tt < 1.4:
            glows.append((prom, 70, 0.45 * math.sin(math.pi * (tt - 0.55) / 0.85), (255, 226, 190)))
        if 0.7 < tt < 3.3:
            glows.append((pos, 62, 0.22 * pa, (255, 170, 210)))
        nb = int(44 * clamp01((tt - 0.7) / 2.3))
        if nb > 1:
            anchor = pos if tt < 3.0 else end + np.array([37, 5]) * ss(3.0, 3.9, tt) + np.array([55, -28]) * (tt - 3.0)
            js = np.arange(nb)
            curl = 0.45 + 0.25 * clamp01(tt - 3.0)
            angs = math.atan2(-0.45, 1.0) + 0.25 * math.sin(seed) + curl * np.sin(js * 0.32 + t * 1.6 + seed) \
                + 0.04 * js * math.sin(seed * 1.7)
            steps = np.column_stack([np.cos(angs), np.sin(angs)]) * 4.6
            chains.append((anchor + np.cumsum(steps, axis=0) - steps[0], 1 - ss(6.6, 8.6, tt)))
        if pa > 0:
            pols.append(("pol2", pos, pa))
    items += pols

    cv = cast_shadow(cv, vf, curve, items, (14, 20), 12, 0.36)
    draw_prot(cv, vf, "coh_b", coh_pos)
    draw_dna(cv, vf, curve, t, PAL["dnaA"], PAL["dnaB"], bright=1 + 0.08 * hl)
    cv = cast_shadow(cv, vf, None, items, (6, 9), 5, 0.42)
    for p, r, a, col in glows:
        draw_glow(cv, vf, p, r, a, col)
    mspr = sphere(max(3, int(round(4.2 * vf.k * 2))), PAL["mrna"], 1.0)
    for pts, a in chains:
        for x, y in vf.pts(pts):
            blit(cv, mspr, x, y, a)
    for name, pos, a in items:
        draw_prot(cv, vf, name, pos, a)
    draw_foreground(cv, vf, t)
    return cv, vf, curve


# ---------------------------------------------------------------- scene 1: the cell
CW, CH = 3200, 1800
NUC_C, NUC_R = np.array([1720.0, 860.0]), 330.0
CEN = np.array([2110.0, 1090.0])
GOL_C = np.array([1990.0, 870.0])
PLATES = [((0, 0, 3200, 1800), 1.25), ((850, 388, 2550, 1344), 2.3), ((1340, 646, 2100, 1074), 5.2)]
LOWQ = 0.25


def cell_shape(x, y):
    dx, dy = (x - 1600) / 1250, (y - 900) / 780
    ang = np.arctan2(dy, dx)
    wob = 1 + 0.045 * np.sin(3 * ang + 0.7) + 0.03 * np.sin(5 * ang + 2.1) + 0.015 * np.sin(9 * ang + 0.3)
    return np.sqrt(dx * dx + dy * dy) / wob


def cell_geometry():
    """All cell detail in cell units, generated once so every plate resolution agrees."""
    rng = np.random.default_rng(21)
    G = {}
    p = rng.uniform((350, 120), (2850, 1680), (30000, 2))
    keep = (cell_shape(p[:, 0], p[:, 1]) < 0.96) & (np.hypot(*(p - NUC_C).T) > NUC_R + 14)
    G["specks"] = [(x, y, rng.uniform(1.0, 1.8), rng.integers(3)) for x, y in p[keep]]

    mts = []
    for a in np.linspace(0, 2 * math.pi, 56, endpoint=False):
        hd, q, line = a + rng.normal(0, 0.05), CEN.copy(), [CEN.copy()]
        kappa = rng.normal(0, 0.004)
        for _ in range(400):
            hd += kappa
            q = q + 8 * np.array([math.cos(hd), math.sin(hd)])
            if cell_shape(q[0], q[1]) > 0.965 or math.hypot(*(q - NUC_C)) < NUC_R + 8:
                break
            line.append(q)
        if len(line) > 4:
            mts.append(np.array(line))
    G["mts"] = mts

    sheets, ribos = [], []
    for k in range(6):
        rk = NUC_R + 52 + k * 34
        a = math.radians(75) + rng.uniform(0, 0.12)
        while a < math.radians(306):
            a1 = min(a + rng.uniform(0.35, 0.8), math.radians(306))
            ang = np.linspace(a, a1, int((a1 - a) * rk / 6) + 2)
            r = rk + 10 * np.sin(3 * ang + k) + 5 * np.sin(7 * ang + 2 * k)
            pts = NUC_C + np.column_stack([np.cos(ang), np.sin(ang)]) * r[:, None]
            sheets.append(pts)
            tg = np.gradient(pts, axis=0)
            nrm = np.column_stack([-tg[:, 1], tg[:, 0]]) / (np.hypot(tg[:, 0], tg[:, 1])[:, None] + 1e-6)
            cum = np.r_[0, np.cumsum(np.hypot(*np.diff(pts, axis=0).T))]
            for s in np.arange(3, cum[-1], 7.0):
                i = min(int(np.searchsorted(cum, s)), len(pts) - 1)
                for side in (-1, 1):
                    ribos.append(pts[i] + nrm[i] * side * 7.6 + rng.normal(0, 0.6, 2))
            a = a1 + rng.uniform(0.04, 0.12)
    bridges = []
    for k in range(-1, 5):
        for _ in range(6):
            a = rng.uniform(math.radians(80), math.radians(300))
            r0 = NUC_R + 2 if k < 0 else NUC_R + 52 + k * 34
            bridges.append(NUC_C + np.outer([r0, NUC_R + 52 + (k + 1) * 34], [math.cos(a), math.sin(a)]))
    nodes = [NUC_C + rng.uniform(NUC_R + 290, NUC_R + 430) * np.array([math.cos(a), math.sin(a)])
             for a in rng.uniform(2.0, 4.4, 36)]
    smooth = []
    for i, nd in enumerate(nodes):
        dists = sorted((math.hypot(*(nd - m)), j) for j, m in enumerate(nodes) if j != i)
        for _, j in dists[:2]:
            if j > i:
                mid = (nd + nodes[j]) / 2 + rng.normal(0, 10, 2)
                smooth.append(catmull([nd, mid, nodes[j]], 8))
    G.update(sheets=sheets, ribos=ribos, bridges=bridges, smooth=smooth)

    cis = []
    for k in range(6):
        rk, hk = 205 + k * 19, math.radians(32 - k * 2.6)
        ang = np.linspace(-hk, hk, 40)
        r = rk + 3 * np.sin(ang * 9 + k)
        cis.append(GOL_C + np.column_stack([np.cos(ang), np.sin(ang)]) * r[:, None])
    gv = []
    for _ in range(18):
        k = rng.integers(6)
        side = rng.choice([-1, 1])
        a = side * (math.radians(32 - k * 2.6) + rng.uniform(0.04, 0.14)) if rng.random() < 0.6 else rng.uniform(-0.5, 0.5)
        rr = 205 + k * 19 + (0 if abs(a) > 0.5 else rng.uniform(120, 150))
        gv.append((GOL_C + rr * np.array([math.cos(a), math.sin(a)]), rng.uniform(5, 8.5)))
    G.update(cis=cis, gves=gv)

    fib = []
    for _ in range(11):
        q = NUC_C + rng.uniform(-0.6, 0.6, 2) * NUC_R
        hd, pts = rng.uniform(0, 2 * math.pi), [q]
        for _ in range(650):
            hd += rng.normal(0, 0.32)
            nq = q + 5.5 * np.array([math.cos(hd), math.sin(hd)])
            if math.hypot(*(nq - NUC_C)) > NUC_R * 0.9:
                hd = math.atan2(*(NUC_C - q)[::-1]) + rng.normal(0, 0.5)
                nq = q + 5.5 * np.array([math.cos(hd), math.sin(hd)])
            q = nq
            pts.append(q)
        fib.append(catmull(np.array(pts)[::3], 6))
    G["fibres"] = fib
    nuc_pts = []
    for _ in range(240):
        v = rng.normal(0, 1, 3)
        v /= np.linalg.norm(v)
        if v[2] > 0.12:
            nuc_pts.append(v)
    G["pores"] = nuc_pts
    ang = np.linspace(0, 2 * math.pi, 60, endpoint=False)
    G["nucleolus"] = np.array([1690, 905]) + np.column_stack([np.cos(ang), np.sin(ang)]) * (88 + 8 * np.sin(3 * ang) + 5 * np.sin(5 * ang + 1))[:, None]
    G["granules"] = [(np.array([1690, 905]) + rng.uniform(0, 80) * np.array([math.cos(a), math.sin(a)]), rng.uniform(2, 4), rng.integers(2))
                     for a in rng.uniform(0, 2 * math.pi, 520)]
    return G


def cell_lowres():
    """Whole-domain low-frequency layers (background with out-of-focus neighbour cells, mottling)."""
    rng = np.random.default_rng(9)
    w, h = int(CW * LOWQ), int(CH * LOWQ)
    y, x = np.mgrid[0:h, 0:w].astype(np.float32) / LOWQ
    g = (y / CH)[..., None]
    img = np.array(PAL["bg_top"], np.float32) * (1 - g) + np.array(PAL["bg_bot"], np.float32) * g
    warm = np.exp(-(((x - 700) / 1500) ** 2 + ((y - 200) / 900) ** 2))[..., None]
    img = img + warm * np.array([60, 52, 42])
    low = Image.fromarray(np.clip(img, 0, 255).astype(np.uint8), "RGB")
    lay = Image.new("RGBA", (w, h), (236, 196, 214, 0))
    d = ImageDraw.Draw(lay)
    for cx, cy, rx, ry, col in ((-250, -150, 900, 620, (226, 186, 210)), (3500, 250, 700, 600, (206, 190, 236)),
                                (3350, 2050, 950, 700, (232, 196, 214)), (-200, 2000, 820, 560, (204, 196, 232))):
        d.ellipse([(cx - rx) * LOWQ, (cy - ry) * LOWQ, (cx + rx) * LOWQ, (cy + ry) * LOWQ], fill=col + (120,))
    lay = lay.filter(ImageFilter.GaussianBlur(28))
    low.paste(lay, (0, 0), lay)
    noise = Image.fromarray((rng.random((h // 6, w // 6)) * 255).astype(np.uint8)).resize((w, h), Image.BICUBIC)
    noise = ImageChops.add(noise.filter(ImageFilter.GaussianBlur(2)),
                           Image.fromarray((rng.random((h // 2, w // 2)) * 120).astype(np.uint8)).resize((w, h), Image.BICUBIC), 1.6)
    return low, noise


def stroke(d, pts, tf, ppu, layers):
    p = [tf(q) for q in pts]
    for col, wu, (ox, oy) in layers:
        d.line([(x + ox * ppu, y + oy * ppu) for x, y in p], fill=col, width=max(1, int(round(wu * ppu))), joint="curve")


def build_cell_plate(reg, ppu, path, G, low, tex):
    x0, y0, x1, y1 = reg
    w, h = int(round((x1 - x0) * ppu)), int(round((y1 - y0) * ppu))
    box = (x0 * LOWQ, y0 * LOWQ, x1 * LOWQ, y1 * LOWQ)
    bg = np.asarray(low.resize((w, h), Image.BICUBIC, box=box), np.float32) / 255
    tx = np.asarray(tex.resize((w, h), Image.BICUBIC, box=box), np.float32) / 255
    X = (x0 + (np.arange(w, dtype=np.float32) + 0.5) / ppu)[None, :]
    Y = (y0 + (np.arange(h, dtype=np.float32) + 0.5) / ppu)[:, None]
    d = cell_shape(X, Y)
    rx, ry = X - 1600, Y - 900
    rr = np.sqrt(rx * rx + ry * ry) + 1e-3
    din = rr * (1 / np.maximum(d, 1e-3) - 1)
    inside = np.clip(din * ppu + 0.5, 0, 1)[..., None]
    dc = np.clip(d, 0, 0.995)
    slope = np.minimum(3 * dc ** 5 / np.sqrt(1 - dc ** 6) * 0.35, 6)
    ux, uy = rx / rr, ry / rr
    nrm = np.sqrt(1 + slope ** 2)
    lit = shade(ux * slope / nrm, uy * slope / nrm, 1 / nrm, np.array([232, 186, 204], np.float32) / 255,
                np.full(d.shape, 0.9, np.float32), rim_k=0.38, spec_k=0.3)
    lit *= 0.9 + 0.2 * tx[..., None]
    lit = lit * 0.88 + bg * 0.12
    lit *= 1 - 0.07 * np.exp(-((din - 26) / 16) ** 2)[..., None]
    lit += (np.exp(-((din - 2.5) / 1.1) ** 2) + np.exp(-((din - 6.5) / 1.1) ** 2))[..., None] * 0.3 * np.array([1.0, 0.96, 1.0])
    lit += (np.exp(-((din - 13) / 7) ** 2) * np.clip(-0.6 * ux - 0.8 * uy, 0, 1) ** 3)[..., None] * 0.55
    bg = bg + (np.exp(-(np.clip(-din, 0, None) / 22) ** 2) * (din < 0))[..., None] * 0.12 * np.array([1.0, 0.92, 0.97])
    img = bg * (1 - inside) + lit * inside
    del bg, tx, d, din, lit, slope, ux, uy, nrm, rx, ry, rr
    im = Image.fromarray((np.clip(img, 0, 1) * 255).astype(np.uint8), "RGB")
    del img
    tf = lambda p: ((p[0] - x0) * ppu, (p[1] - y0) * ppu)
    inreg = lambda p, m=30: x0 - m < p[0] < x1 + m and y0 - m < p[1] < y1 + m
    dr = ImageDraw.Draw(im)

    spk = [(214, 160, 184), (232, 190, 206), (252, 236, 242)]
    for x, y, r, c in G["specks"]:
        if inreg((x, y), 3):
            X0, Y0 = tf((x, y))
            rp = max(0.6, r * ppu)
            dr.ellipse([X0 - rp, Y0 - rp, X0 + rp, Y0 + rp], fill=spk[c])

    ov = Image.new("RGBA", im.size, (240, 248, 255, 0))
    od = ImageDraw.Draw(ov)
    for line in G["mts"]:
        p = [tf(q) for q in line]
        od.line(p, fill=(236, 246, 255, 28), width=max(1, int(5 * ppu)), joint="curve")
        od.line(p, fill=(248, 252, 255, 85), width=max(1, int(1.6 * ppu)), joint="curve")
    im.paste(ov, (0, 0), ov)
    del ov

    dr = ImageDraw.Draw(im)
    er_layers = [((198, 150, 182), 17, (1.6, 2.4)), ((140, 126, 200), 15, (0, 0)), ((206, 196, 248), 11, (0, 0)),
                 ((246, 242, 255), 3, (-1.4, -2.0))]
    for pts in G["smooth"]:
        if inreg(pts[0], 200):
            stroke(dr, pts, tf, ppu, [((198, 150, 182), 9, (1.2, 1.8)), ((160, 146, 212), 7.5, (0, 0)),
                                      ((220, 212, 252), 5, (0, 0)), ((248, 246, 255), 1.6, (-0.8, -1.0))])
    for pts in G["bridges"]:
        stroke(dr, pts, tf, ppu, er_layers[1:3])
    for pts in G["sheets"]:
        if inreg(pts[len(pts) // 2], 400):
            stroke(dr, pts, tf, ppu, er_layers)
    for q in G["ribos"]:
        if inreg(q, 4):
            X0, Y0 = tf(q)
            rp = max(0.7, 2.1 * ppu)
            dr.ellipse([X0 - rp, Y0 - rp, X0 + rp, Y0 + rp], fill=(112, 92, 168))
            dr.ellipse([X0 - rp * 0.45 - rp * 0.25, Y0 - rp * 0.45 - rp * 0.3, X0 + rp * 0.45 - rp * 0.25, Y0 + rp * 0.45 - rp * 0.3],
                       fill=(176, 160, 222))
    gol_layers = [((200, 150, 160), 15, (1.6, 2.4)), ((214, 146, 72), 13.5, (0, 0)), ((255, 208, 136), 9.5, (0, 0)),
                  ((255, 242, 206), 2.6, (-1.2, -1.8))]
    for pts in G["cis"]:
        stroke(dr, pts, tf, ppu, gol_layers)
        for end in (pts[0], pts[-1]):
            X0, Y0 = tf(end)
            for col, rad in (((214, 146, 72), 9.5), ((255, 208, 136), 7.5), ((255, 242, 206), 2.5)):
                rp = rad * ppu
                dr.ellipse([X0 - rp, Y0 - rp, X0 + rp, Y0 + rp], fill=col)
    for q, rad in G["gves"]:
        X0, Y0 = tf(q)
        spr = sphere(max(3, int(2 * rad * ppu)), (255, 214, 150) if rad > 6.5 else (226, 210, 255), 1.0)
        im.paste(spr, (int(X0 - spr.width / 2), int(Y0 - spr.height / 2)), spr)

    halo = soft_disc(max(2, int(46 * ppu)), (250, 244, 255), max(1, int(18 * ppu)))
    X0, Y0 = tf(CEN)
    im.paste(halo, (int(X0 - halo.width / 2), int(Y0 - halo.height / 2)), halo.getchannel("A").point(lambda v: v // 2))
    for (ox, oy, wu, hu) in ((-12, -3, 18, 7), (8, -9, 7, 18)):
        a, b = tf(CEN + np.array([ox, oy])), tf(CEN + np.array([ox + wu, oy + hu]))
        dr.rounded_rectangle([a, b], radius=3 * ppu, fill=(214, 228, 244), outline=(150, 166, 204), width=max(1, int(1.2 * ppu)))
        dr.line([(a[0] + 2 * ppu, a[1] + 2 * ppu), (b[0] - 2 * ppu, a[1] + 2 * ppu)], fill=(244, 248, 255), width=max(1, int(1.0 * ppu)))

    # nucleus
    nb = [max(0, int((NUC_C[0] - NUC_R - 20 - x0) * ppu)), max(0, int((NUC_C[1] - NUC_R - 20 - y0) * ppu)),
          min(w, int((NUC_C[0] + NUC_R + 20 - x0) * ppu)), min(h, int((NUC_C[1] + NUC_R + 20 - y0) * ppu))]
    if nb[2] > nb[0] and nb[3] > nb[1]:
        sub = im.crop(nb)
        arr = np.asarray(sub, np.float32) / 255
        Xn = (x0 + (np.arange(nb[0], nb[2], dtype=np.float32) + 0.5) / ppu)[None, :]
        Yn = (y0 + (np.arange(nb[1], nb[3], dtype=np.float32) + 0.5) / ppu)[:, None]
        ndx, ndy = (Xn - NUC_C[0]) / NUC_R, (Yn - NUC_C[1]) / NUC_R
        nd = np.sqrt(ndx * ndx + ndy * ndy)
        nin = np.clip((1 - nd) * NUC_R * ppu + 0.5, 0, 1)[..., None]
        dd = np.maximum(nd, 1.0)
        nx, ny = ndx / dd, ndy / dd
        nz = np.sqrt(np.clip(1 - nx * nx - ny * ny, 0, 1))
        lit = shade(nx, ny, nz, np.array([176, 154, 224], np.float32) / 255, 0.8 + 0.2 * nz, rim_k=0.5, spec_k=0.05)
        dn = (1 - nd) * NUC_R
        lit *= 1 - 0.10 * np.exp(-((dn - 22) / 14) ** 2)[..., None]
        lit *= 1 - 0.06 * ((dn > 3.2) & (dn < 5.4))[..., None]
        lit += (np.exp(-((dn - 2.2) / 1.2) ** 2) + np.exp(-((dn - 6.4) / 1.2) ** 2))[..., None] * 0.28
        arr = arr * (1 - nin) + lit * nin
        sub = Image.fromarray((np.clip(arr, 0, 1) * 255).astype(np.uint8), "RGB")
        ox_, oy_ = nb[0], nb[1]
        stf = lambda p: (tf(p)[0] - ox_, tf(p)[1] - oy_)
        ov = Image.new("RGBA", sub.size, (230, 222, 255, 0))
        od = ImageDraw.Draw(ov)
        for f in G["fibres"]:
            p = [stf(q) for q in f]
            od.line(p, fill=(150, 126, 206, 120), width=max(1, int(4.6 * ppu)), joint="curve")
            od.line(p, fill=(226, 214, 252, 150), width=max(1, int(2.2 * ppu)), joint="curve")
            for q in p[::3]:
                rp = 2.3 * ppu
                od.ellipse([q[0] - rp, q[1] - rp, q[0] + rp, q[1] + rp], fill=(240, 234, 255, 160))
        od.polygon([stf(q) for q in G["nucleolus"]], fill=(150, 118, 206, 225))
        for q, rad, c in G["granules"]:
            X0, Y0 = stf(q)
            rp = rad * ppu
            od.ellipse([X0 - rp, Y0 - rp, X0 + rp, Y0 + rp], fill=((118, 88, 180, 210), (190, 168, 236, 210))[c])
        for v in G["pores"]:
            X0, Y0 = stf(NUC_C + NUC_R * v[:2])
            rp = 7.5 * ppu * (0.45 + 0.55 * v[2])
            od.ellipse([X0 - rp, Y0 - rp, X0 + rp, Y0 + rp], outline=(126, 104, 184, 230), width=max(1, int(2.2 * ppu)))
            rc = rp * 0.32
            od.ellipse([X0 - rc, Y0 - rc, X0 + rc, Y0 + rc], fill=(238, 230, 255, 230))
        ov = ov.filter(ImageFilter.GaussianBlur(0.35 * ppu))
        sub.paste(ov, (0, 0), ov)
        arr = np.asarray(sub, np.float32) / 255
        sheen = (np.clip(nx * HALF[0] + ny * HALF[1] + nz * HALF[2], 0, 1) ** 60 * 0.12
                 + np.exp(-(((nx + 0.42) ** 2 + (ny + 0.5) ** 2) / 0.07)) * 0.035)[..., None]
        arr = arr + sheen * nin
        im.paste(Image.fromarray((np.clip(arr, 0, 1) * 255).astype(np.uint8), "RGB"), (nb[0], nb[1]))
    im.save(path)


MPPU_4K = 6.8


def organelles():
    rng = np.random.default_rng(5)
    out, placed = [], []

    def ok(p, rad):
        if cell_shape(p[0], p[1]) > 0.86 or math.hypot(*(p - NUC_C)) < NUC_R + 250 + rad:
            return False
        if math.hypot(*(p - np.array([2230, 870]))) < 170 + rad or math.hypot(*(p - CEN)) < 80 + rad:
            return False
        return all(math.hypot(*(p - q)) > rad + r2 + 14 for q, r2 in placed)

    for kind, count in (("mito", 11), ("lyso", 6), ("perox", 5), ("ves", 22)):
        n = 0
        while n < count:
            p = rng.uniform((350, 150), (2850, 1650))
            if kind == "mito":
                prm = (round(rng.uniform(130, 190)), round(rng.uniform(50, 60)), round(rng.uniform(0, math.pi), 2))
                rad = prm[0] / 2
            else:
                rad = {"lyso": rng.uniform(16, 24), "perox": rng.uniform(11, 15), "ves": rng.uniform(6, 12)}[kind]
                prm = (round(rad, 1), int(rng.integers(3)))
            if ok(p, rad):
                placed.append((p, rad))
                out.append((kind, p, prm, rng.uniform(0, 6)))
                n += 1
    return out


ORG = organelles()


@lru_cache(maxsize=None)
def mito_sprite(L, wd, ang, seed):
    ppu = MPPU_4K * S / 3
    r, seg = wd / 2, L / 2 - wd / 2
    n = int(2 * (L / 2 + 4) * ppu) + 4
    g = (np.arange(n, dtype=np.float32) + 0.5) / ppu - n / (2 * ppu)
    X, Y = g[None, :], g[:, None]
    ca, sa = math.cos(ang), math.sin(ang)
    lx, ly = X * ca + Y * sa, -X * sa + Y * ca
    dxl = lx - np.clip(lx, -seg, seg)
    dist = np.sqrt(dxl * dxl + ly * ly)
    din = r - dist
    tt = np.clip(dist / r, 0, 0.999)
    inv = 1 / np.maximum(dist, 1e-4)
    nlx, nly = dxl * inv * tt, ly * inv * tt
    nz = np.sqrt(1 - tt * tt)
    nx, ny = nlx * ca - nly * sa, nlx * sa + nly * ca
    col = np.empty(dist.shape + (3,), np.float32)
    col[:] = np.array([236, 124, 112]) / 255
    rng = np.random.default_rng(seed)
    for k, xk in enumerate(np.arange(-seg - r * 0.45 + 7, seg + r * 0.45 - 6, 15.0)):
        side = 1 if k % 2 == 0 else -1
        dxk = np.abs(lx - xk - 2.2 * np.sin(ly * 0.22 + k * 1.7))
        m = (dxk < 3.4) & (side * ly > -rng.uniform(0.25, 0.6) * r) & (din > 4.5)
        prof = (1 - (dxk / 3.4) ** 2)[..., None]
        fold = np.array([255, 214, 192], np.float32) / 255 * (0.72 + 0.28 * prof)
        lumen = (dxk < 0.9)[..., None]
        fold = np.where(lumen, np.array([228, 146, 132], np.float32) / 255, fold)
        col = np.where(m[..., None], fold, col)
    for lo, hi, c in ((-1, 2.6, (255, 204, 182)), (2.6, 4.2, (246, 158, 140)), (4.2, 5.6, (255, 214, 192))):
        col = np.where(((din >= lo) & (din < hi))[..., None], np.array(c, np.float32) / 255, col)
    lit = shade(nx, ny, nz, col, 0.75 + 0.25 * nz, rim_k=0.45, spec_k=0.45)
    return to_rgba(lit, din * ppu + 0.5)


@lru_cache(maxsize=None)
def ball_sprite(kind, rad, var):
    ppu = MPPU_4K * S / 3
    n = int(2 * rad * ppu) + 4
    g = (np.arange(n, dtype=np.float32) + 0.5 - n / 2) / (rad * ppu)
    X, Y = g[None, :], g[:, None]
    d = np.sqrt(X * X + Y * Y)
    dd = np.maximum(d, 1.0)
    nx, ny = X / dd, Y / dd
    nz = np.sqrt(np.clip(1 - nx * nx - ny * ny, 0, 1))
    col = {"lyso": (198, 150, 226), "perox": (214, 238, 190),
           "ves": [(236, 228, 255), (255, 224, 236), (222, 242, 255)][var]}[kind]
    base = np.empty(d.shape + (3,), np.float32)
    base[:] = np.array(col, np.float32) / 255
    rng = np.random.default_rng(int(rad * 10) + var)
    if kind == "lyso":
        for _ in range(26):
            cx, cy = rng.uniform(-0.6, 0.6, 2)
            base = np.where((((X - cx) ** 2 + (Y - cy) ** 2) < rng.uniform(0.004, 0.012))[..., None],
                            np.array([140, 92, 176], np.float32) / 255, base)
    if kind == "perox":
        base = np.where(((np.abs(X) < 0.32) & (np.abs(Y) < 0.28))[..., None], np.array([150, 194, 138], np.float32) / 255, base)
    lit = shade(nx, ny, nz, base, 0.7 + 0.3 * nz, rim_k=0.5)
    return to_rgba(lit, (1 - d) * rad * ppu + 0.5)


@lru_cache(maxsize=None)
def plate_cell(i):
    return Image.open(os.path.join(CACHE, f"cell_{i}.png")).copy()


def sample_cell(box):
    x0, y0, x1, y1 = box
    best, blend = 0, 1.0
    for i, ((rx0, ry0, rx1, ry1), _) in enumerate(PLATES):
        m = min(x0 - rx0, y0 - ry0, rx1 - x1, ry1 - y1)
        if m >= 0:
            best, blend = i, ss(0, 0.04 * (rx1 - rx0), m)

    def crop(i):
        (rx0, ry0, _, _), ppu = PLATES[i]
        ppu *= S / 3
        im = plate_cell(i)
        bx = (max(0, (x0 - rx0) * ppu), max(0, (y0 - ry0) * ppu), min(im.width, (x1 - rx0) * ppu), min(im.height, (y1 - ry0) * ppu))
        return im.resize((W, H), Image.BICUBIC, box=bx)

    img = crop(best)
    if best > 0 and blend < 1:
        img = Image.blend(crop(best - 1), img, blend)
    return img


FLY = np.random.default_rng(17).uniform(0, 1, (9, 5))


def cell_scene(t):
    z = math.exp(math.log(5.6) * clamp01(t / 3.35) ** 1.6)
    k = ss(0, 3.0, t)
    cx, cy = lerp(1600, NUC_C[0], k), lerp(900, NUC_C[1], k)
    cw, ch = CW / z, CH / z
    x0, y0 = cx - cw / 2, cy - ch / 2
    f = W / cw
    cv = sample_cell((x0, y0, x0 + cw, y0 + ch))
    mppu = MPPU_4K * S / 3
    for kind, p, prm, ph in ORG:
        q = p + 10 * np.array([math.sin(t * 0.7 + ph), math.cos(t * 0.5 + ph)])
        sx, sy = (q[0] - x0) * f, (q[1] - y0) * f
        ext = (prm[0] / 2 if kind == "mito" else prm[0]) * f + 4
        if not (-ext < sx < W + ext and -ext < sy < H + ext):
            continue
        base = mito_sprite(*prm, int(ph * 100)) if kind == "mito" else ball_sprite(kind, *prm)
        sc = f / mppu
        spr = base.convert("RGBa").resize((max(2, round(base.width * sc)), max(2, round(base.height * sc))), Image.LANCZOS).convert("RGBA")
        blit(cv, spr, sx, sy, 0.97)
    b = ss(2.4, 3.3, t)
    if b > 0.02:
        half = cv.resize((W // 2, H // 2), Image.BILINEAR).filter(ImageFilter.GaussianBlur(b * 3 * S))
        cv = half.resize((W, H), Image.BILINEAR)
    dolly = 1.25 * ss(0, 3.35, t)
    cols = [(255, 226, 236), (226, 216, 255), (255, 214, 190)]
    for ox, oy, dep, rr, ci in FLY:
        s = 1 / (1.5 + 2.2 * dep - dolly)
        x, y = W / 2 + (ox - 0.5) * 1.4 * W * s, H / 2 + (oy - 0.5) * 1.4 * H * s
        r = (40 + 70 * rr) * S * s
        draw_soft(cv, x, y, r, cols[int(ci * 3)], 0.28 * (1 - ss(2.8, 3.35, t)), 0.35)
    return cv, ((NUC_C[0] - x0) * f, (NUC_C[1] - y0) * f, NUC_R * f)


# ---------------------------------------------------------------- text and post
@lru_cache(maxsize=16)
def font(path, px):
    return ImageFont.truetype(path, px)


CAPTIONS = [
    (12.4, 16.0, "Enhancer elements can increase\ntranscriptional output"),
    (16.6, 20.2, "Repressor elements can decrease\ntranscriptional output"),
    (21.0, 28.8, "Our stapler molecule enforces\nco-localisation of promoters\nand regulatory elements"),
]


def white_elements(cv, box, draw_fn, alpha):
    """Draw white shapes (via draw_fn on an L mask) with a soft drop shadow onto cv."""
    x0, y0, x1, y1 = (int(v) for v in box)
    x0, y0, x1, y1 = max(0, x0), max(0, y0), min(W, x1), min(H, y1)
    if x1 <= x0 or y1 <= y0:
        return
    m = Image.new("L", (x1 - x0, y1 - y0), 0)
    draw_fn(ImageDraw.Draw(m), -x0, -y0)
    sh = m.filter(ImageFilter.GaussianBlur(3.5 * S)).point(lambda v: int(v * 0.6 * alpha))
    dx, dy = int(2 * S), int(3 * S)
    shb = (x0 + dx, y0 + dy, min(W, x1 + dx), min(H, y1 + dy))
    cv.paste((14, 20, 32), shb, sh.crop((0, 0, shb[2] - shb[0], shb[3] - shb[1])))
    cv.paste((255, 255, 255), (x0, y0, x1, y1), m.point(lambda v: int(v * alpha)))


def overlay(cv, t, view, curve):
    fb = font(FONT_BOLD, int(44 * S))
    for a, b, txt in CAPTIONS:
        al = ss(a, a + 0.5, t) * (1 - ss(b - 0.5, b, t))
        if al > 0:
            pos, sp = (70 * S, 52 * S), int(10 * S)
            bb = ImageDraw.Draw(Image.new("L", (1, 1))).multiline_textbbox(pos, txt, font=fb, spacing=sp)
            white_elements(cv, (bb[0] - 20 * S, bb[1] - 20 * S, bb[2] + 20 * S, bb[3] + 20 * S),
                           lambda d, ox, oy: d.multiline_text((pos[0] + ox, pos[1] + oy), txt, font=fb, fill=255, spacing=sp), al)
    al = ss(7.6, 8.1, t) * (1 - ss(11.4, 11.9, t))
    if al > 0 and curve is not None:
        gp = np.array(view.pt1(at_u(curve, 11.9)))
        lp = gp + np.array([100, 74]) * S
        r = 9 * S
        e = gp + unit(lp - gp) * r
        fl = font(FONT_BOLD, int(42 * S))

        def draw(d, ox, oy):
            o = np.array([ox, oy])
            g, ee, l = gp + o, e + o, lp + o
            d.ellipse([g[0] - r, g[1] - r, g[0] + r, g[1] + r], outline=255, width=max(1, int(2.6 * S)))
            d.ellipse([g[0] - r / 3, g[1] - r / 3, g[0] + r / 3, g[1] + r / 3], fill=255)
            d.line([tuple(ee), tuple(l), (l[0] + 24 * S, l[1])], fill=255, width=max(1, int(2.4 * S)))
            d.text((l[0] + 34 * S, l[1] - 30 * S), "Gene", font=fl, fill=255)

        white_elements(cv, (gp[0] - 20 * S, gp[1] - 20 * S, lp[0] + 200 * S, lp[1] + 50 * S), draw, al)
    al = ss(6.5, 7.5, t)
    if al > 0:
        fs = font(FONT_SEMI, int(17 * S))
        txt = "Illustrative, not to scale. In silico concept."
        tw = fs.getbbox(txt)[2]
        x, y = W - 36 * S - tw, H - 46 * S
        white_elements(cv, (x - 10 * S, y - 10 * S, x + tw + 10 * S, y + 34 * S),
                       lambda d, ox, oy: d.text((x + ox, y + oy), txt, font=fs, fill=255), 0.72 * al)
    return cv


@lru_cache(maxsize=1)
def vignette():
    y, x = np.mgrid[0:H, 0:W].astype(np.float32)
    d = np.sqrt(((x - W / 2) / (W / 2)) ** 2 + ((y - H / 2) / (H / 2)) ** 2)
    v = np.clip(1 - 0.30 * np.clip(d - 0.5, 0, None) ** 1.5, 0, 1)
    g = Image.fromarray((v * 255).astype(np.uint8), "L")
    return Image.merge("RGB", (g, g, g))


def _grade_lut():
    lut = []
    for gain, lift in ((1.03, 0.0), (1.0, 0.0), (0.98, 0.018)):
        for v in range(256):
            x = v / 255
            y = (x * 0.7 + x * x * (3 - 2 * x) * 0.3) * gain + lift * (1 - x)
            lut.append(int(min(1, max(0, y)) * 255 + 0.5))
    return lut


GRADE = _grade_lut()
BRIGHT = [0 if v < 205 else min(255, int((v - 205) * 4.0)) for v in range(256)] * 3


def post(cv, i):
    sm = cv.resize((W // 6, H // 6), Image.BILINEAR).point(BRIGHT)
    bloom = ImageChops.add(sm.filter(ImageFilter.GaussianBlur(3 * S / 3 * 2)), sm.filter(ImageFilter.GaussianBlur(12 * S / 3 * 2)))
    bloom = bloom.point(lambda v: int(v * 0.55)).resize((W, H), Image.BILINEAR)
    cv = ImageChops.screen(cv, bloom).point(GRADE)
    cv = ImageChops.multiply(cv, vignette())
    g = Image.fromarray(np.random.default_rng(i).integers(0, 11, (H, W), dtype=np.uint8), "L")
    return ImageChops.add(cv, Image.merge("RGB", (g, g, g)), 1.0, -5)


def render(i):
    t = i / FPS
    path = os.path.join(CACHE, "frames", f"f_{i:04d}.jpg")
    if os.path.exists(path):
        return
    view = curve = None
    if t < 3.35:
        cv, (nx, ny, nr) = cell_scene(t)
        a2 = ss(2.2, 3.0, t)
        if a2 > 0:
            tad, view, curve = tad_scene(t)
            nr *= 1 + 0.8 * ss(2.8, 3.35, t)
            q = 4
            mask = Image.new("L", (W // q, H // q), 0)
            ImageDraw.Draw(mask).ellipse([(nx - nr) / q, (ny - nr) / q, (nx + nr) / q, (ny + nr) / q], fill=int(255 * a2))
            mask = mask.filter(ImageFilter.GaussianBlur(6 * S / q)).resize((W, H), Image.BILINEAR)
            cv = Image.composite(tad, cv, mask)
            ra = 0.65 * ss(2.2, 2.7, t) * (1 - ss(3.05, 3.35, t))
            if ra > 0:
                ring = Image.new("L", (W // q, H // q), 0)
                ImageDraw.Draw(ring).ellipse([(nx - nr) / q, (ny - nr) / q, (nx + nr) / q, (ny + nr) / q],
                                             outline=int(255 * ra), width=max(1, int(10 * S / q)))
                ring = ring.filter(ImageFilter.GaussianBlur(5 * S / q)).resize((W, H), Image.BILINEAR)
                tint = Image.merge("RGB", [ring.point(lambda v, c=c: v * c // 255) for c in (240, 228, 255)])
                cv = ImageChops.screen(cv, tint)
    else:
        cv, view, curve = tad_scene(t)
    if t > 6:
        cv = overlay(cv, t, view, curve)
    post(cv, i).save(path, quality=96, subsampling=0)


def build_plates(cache):
    jobs = [("far.png", lambda p: build_far(p)), ("mid.png", lambda p: build_mid(p))]
    need_cell = any(not os.path.exists(os.path.join(cache, f"cell_{i}.png")) for i in range(len(PLATES)))
    for name, fn in jobs:
        if not os.path.exists(os.path.join(cache, name)):
            print("building", name, flush=True)
            fn(os.path.join(cache, name))
    if need_cell:
        G = cell_geometry()
        low, tex = cell_lowres()
        for i, (reg, ppu) in enumerate(PLATES):
            path = os.path.join(cache, f"cell_{i}.png")
            if not os.path.exists(path):
                print("building", path, flush=True)
                build_cell_plate(reg, ppu * S / 3, path, G, low, tex)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--width", type=int, default=1920)
    ap.add_argument("--out", default="animation/renders/procedural_1080p")
    ap.add_argument("--frames", default="", help="comma list of frame indices for stills, e.g. 0,100,300")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    args = ap.parse_args()
    cache = os.path.abspath(args.out)
    os.makedirs(os.path.join(cache, "frames"), exist_ok=True)
    init(args.width, cache)
    build_plates(cache)
    idx = [int(v) for v in args.frames.split(",")] if args.frames else list(range(NFRAMES))
    with Pool(args.workers, initializer=init, initargs=(args.width, cache)) as pool:
        for k, _ in enumerate(pool.imap_unordered(render, idx, chunksize=2)):
            if k % 24 == 0:
                print(f"{k}/{len(idx)} frames", flush=True)
    if not args.frames:
        import imageio_ffmpeg
        ff = imageio_ffmpeg.get_ffmpeg_exe()
        mp4 = cache.rstrip("/\\") + ".mp4"
        src = os.path.join(cache, "frames", "f_%04d.jpg")
        subprocess.run([ff, "-y", "-framerate", str(FPS), "-i", src, "-c:v", "libx264", "-pix_fmt", "yuv420p",
                        "-crf", "16", "-preset", "slow", "-tune", "film", "-colorspace", "bt709", "-color_primaries", "bt709",
                        "-color_trc", "bt709", "-movflags", "+faststart", mp4], check=True)
        print("wrote", mp4)


if __name__ == "__main__":
    main()

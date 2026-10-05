"""Animation A1: MRF fingerprints (IR-FISP, EPG) for WM / GM / CSF at 3T.

Renders frames with matplotlib at native pixel size, then encodes
GIF (two-pass palette) and MP4 with ffmpeg.

Output: SCRATCH/assets/anim/anim_fingerprint.{gif,mp4}
"""
import os
import shutil
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a1_fisp_epg import simulate, N_TR  # noqa: E402

SCRATCH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(SCRATCH, "assets", "anim")
FRAME_DIR = os.path.join(SCRATCH, "anim_src", "_a1_frames")
NAME = "anim_fingerprint"

# ---------------------------------------------------------------- style
INK = "#1B1B1F"
MUTED = "#5B5F66"
HAIR = "#D9DADC"
RED = "#A8322D"
BLUE = "#1F4E79"
LBLUE = "#8DB3D9"
LRED = "#E6B3AE"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Liberation Serif"],
    "mathtext.fontset": "stix",
    "font.size": 16,
    "axes.edgecolor": MUTED,
    "axes.linewidth": 1.0,
    "axes.labelcolor": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK,
    "ytick.labelcolor": INK,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "axes.labelsize": 18,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})

W, H, DPI = 1100, 620, 100
FPS = 15
T_INTRO = 0.6      # s, empty axes with inversion marker
T_DRAW = 10.0      # s, progressive reveal of 500 TRs
T_HOLD = 2.0       # s, final frame hold

# ---------------------------------------------------------------- layout (pixels)
L_PX, R_PX = 104, 936          # axes left / right edge
FA_TOP, FA_H = 606, 98
TR_GAP, TR_H = 22, 52
MAIN_GAP, MAIN_BOT = 26, 86


def px_axes(fig, left, bottom, width, height):
    return fig.add_axes([left / W, bottom / H, width / W, height / H])


def place_labels(y, gap, lo, hi):
    """1-D label de-overlap: keep order, enforce min spacing, clamp to [lo, hi]."""
    order = np.argsort(y)
    p = np.array(y, float)[order]
    for _ in range(200):
        moved = False
        for i in range(1, len(p)):
            d = p[i] - p[i - 1]
            if d < gap - 1e-12:
                s = (gap - d) / 2
                p[i - 1] -= s
                p[i] += s
                moved = True
        if p[0] < lo:
            p += lo - p[0]
        if p[-1] > hi:
            p -= p[-1] - hi
        if not moved:
            break
    out = np.empty_like(p)
    out[order] = p
    return out


def build():
    d = simulate()
    fa, tr = d["fa"], d["tr"]
    curves = d["curves"]
    n = N_TR
    x = np.arange(n)

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    width = R_PX - L_PX
    fa_bot = FA_TOP - FA_H
    tr_top = fa_bot - TR_GAP
    tr_bot = tr_top - TR_H
    main_top = tr_bot - MAIN_GAP
    ax_fa = px_axes(fig, L_PX, fa_bot, width, FA_H)
    ax_tr = px_axes(fig, L_PX, tr_bot, width, TR_H)
    ax = px_axes(fig, L_PX, MAIN_BOT, width, main_top - MAIN_BOT)

    xlim = (-14, n + 2)
    for a in (ax_fa, ax_tr, ax):
        a.set_xlim(*xlim)
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)

    # --- flip-angle strip
    ax_fa.set_ylim(0, 76)
    ax_fa.set_yticks([0, 35, 70])
    ax_fa.set_xticks([])
    ax_fa.spines["bottom"].set_color(HAIR)
    segs = np.stack([np.stack([x, np.zeros(n)], 1), np.stack([x, fa], 1)], 1)
    fa_bars = LineCollection([], colors=LBLUE, linewidths=1.05)
    ax_fa.add_collection(fa_bars)
    fa_head = LineCollection([], colors=BLUE, linewidths=2.2)
    ax_fa.add_collection(fa_head)
    # inversion pulse (180 deg, off scale): red bar before the train + label
    inv_x = -7
    ax_fa.plot([inv_x, inv_x], [0, 76], color=RED, lw=2.4, solid_capstyle="butt", clip_on=False)
    inv_txt = ax_fa.text(inv_x + 6, 73, "180° inversion", color=RED, fontsize=15,
                         ha="left", va="top")
    fig.text((R_PX + 14) / W, (fa_bot + FA_H * 0.5) / H, "Flip angle\n(deg)",
             ha="left", va="center", fontsize=17, color=INK, linespacing=1.05)

    # --- TR strip
    ax_tr.set_ylim(10.4, 15.6)
    ax_tr.set_yticks([11, 15])
    ax_tr.set_xticks([])
    ax_tr.spines["bottom"].set_visible(False)
    tr_dots, = ax_tr.plot([], [], ls="none", marker="o", ms=2.3, mfc=MUTED, mec="none")
    tr_head, = ax_tr.plot([], [], ls="none", marker="o", ms=6, mfc=INK, mec="none")
    fig.text((R_PX + 14) / W, (tr_bot + TR_H * 0.5) / H, "TR (ms)",
             ha="left", va="center", fontsize=17, color=INK)

    # --- main panel
    ymax = 0.19
    ax.set_ylim(0, ymax)
    ax.set_yticks([0, 0.05, 0.10, 0.15])
    ax.set_yticklabels(["0", "0.05", "0.10", "0.15"])
    ax.set_xticks([0, 100, 200, 300, 400, 500])
    ax.set_xlabel("Time point (TR index)", labelpad=6)
    ax.set_ylabel("Signal (a.u.)", labelpad=10)
    ax.grid(axis="y", color="#EFEFF0", lw=0.8)
    ax.set_axisbelow(True)

    lines, heads, labels, leaders = [], [], [], []
    for name, T1, T2, col, s in curves:
        ln, = ax.plot([], [], color=col, lw=2.0, solid_joinstyle="round", zorder=3)
        hd, = ax.plot([], [], ls="none", marker="o", ms=7, mfc=col, mec="white", mew=1.0, zorder=4)
        lb = ax.text(0, 0, name, color=col, fontsize=18, ha="left", va="center",
                     clip_on=False, zorder=5)
        lb.set_visible(False)
        ld, = ax.plot([], [], color=col, lw=1.0, zorder=2, clip_on=False)
        lines.append(ln)
        heads.append(hd)
        labels.append(lb)
        leaders.append(ld)

    # cursor through all three panels
    cursors = [a.axvline(0, color="#B9BBBF", lw=0.9, zorder=1) for a in (ax_fa, ax_tr, ax)]

    # label spacing in data units (main axes height in px)
    main_h_px = main_top - MAIN_BOT
    gap = 25.0 * ymax / main_h_px
    sig = np.array([c[4] for c in curves])

    n_intro = int(round(T_INTRO * FPS))
    n_draw = int(round(T_DRAW * FPS))
    n_hold = int(round(T_HOLD * FPS))
    total = n_intro + n_draw + n_hold

    os.makedirs(FRAME_DIR, exist_ok=True)
    for f in os.listdir(FRAME_DIR):
        os.remove(os.path.join(FRAME_DIR, f))

    lab_y = None
    lab_anchor = None
    for fi in range(total):
        if fi < n_intro:
            k = 0
        elif fi < n_intro + n_draw:
            u = (fi - n_intro + 1) / n_draw
            k = int(np.ceil(u * n))
        else:
            k = n
        k = max(0, min(n, k))
        final = fi >= n_intro + n_draw

        fa_bars.set_segments(segs[:k])
        fa_head.set_segments(segs[k - 1:k] if (k > 0 and not final) else [])
        tr_dots.set_data(x[:k], tr[:k])
        if k > 0 and not final:
            tr_head.set_data([x[k - 1]], [tr[k - 1]])
        else:
            tr_head.set_data([], [])

        if k > 0:
            # Smooth the label *anchors* (not the placed positions), then
            # de-overlap from scratch every frame: labels never overlap; when
            # two curves cross, their labels swap places in one step.
            head_y = sig[:, k - 1]
            if lab_anchor is None:
                lab_anchor = head_y.copy()
            else:
                lab_anchor += 0.3 * (head_y - lab_anchor)
            if final and fi >= n_intro + n_draw + 8:
                lab_anchor = head_y.copy()
            lab_y = place_labels(lab_anchor, gap, 0.012, ymax - 0.008)
        for j in range(len(curves)):
            lines[j].set_data(x[:k], sig[j, :k])
            if k > 0:
                heads[j].set_data([x[k - 1]], [sig[j, k - 1]])
                labels[j].set_visible(True)
                labels[j].set_position((x[k - 1] + 8, lab_y[j]))
                yh = sig[j, k - 1]
                if abs(lab_y[j] - yh) > 0.35 * gap:
                    leaders[j].set_data([x[k - 1] + 2.5, x[k - 1] + 7], [yh, lab_y[j]])
                else:
                    leaders[j].set_data([], [])
            else:
                heads[j].set_data([], [])
                labels[j].set_visible(False)
                leaders[j].set_data([], [])
        for c in cursors:
            c.set_xdata([k - 1, k - 1])
            c.set_visible(0 < k and not final)

        fig.canvas.draw()
        buf = np.asarray(fig.canvas.buffer_rgba())[..., :3]
        Image.fromarray(buf).save(os.path.join(FRAME_DIR, f"f_{fi:04d}.png"))
    plt.close(fig)
    return total, d


def encode():
    os.makedirs(OUT_DIR, exist_ok=True)
    pat = os.path.join(FRAME_DIR, "f_%04d.png")
    pal = os.path.join(FRAME_DIR, "palette.png")
    gif = os.path.join(OUT_DIR, NAME + ".gif")
    mp4 = os.path.join(OUT_DIR, NAME + ".mp4")
    ff = ["/usr/bin/ffmpeg", "-y", "-loglevel", "error"]
    subprocess.run(ff + ["-framerate", str(FPS), "-i", pat,
                         "-vf", "palettegen=stats_mode=diff:max_colors=96", pal], check=True)
    subprocess.run(ff + ["-framerate", str(FPS), "-i", pat, "-i", pal,
                         "-lavfi", "paletteuse=dither=bayer:bayer_scale=5:diff_mode=rectangle",
                         "-loop", "0", gif], check=True)
    subprocess.run(ff + ["-framerate", str(FPS), "-i", pat,
                         "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20",
                         "-preset", "slow", "-movflags", "+faststart", mp4], check=True)
    return gif, mp4


def export_stills(indices, gif):
    """Extract frames from the encoded GIF (so palette/dither artifacts are visible)."""
    fdir = os.path.join(OUT_DIR, "frames")
    os.makedirs(fdir, exist_ok=True)
    for f in os.listdir(fdir):
        if f.startswith(NAME + "_f"):
            os.remove(os.path.join(fdir, f))
    im = Image.open(gif)
    out = []
    for i in indices:
        im.seek(i)
        dst = os.path.join(fdir, f"{NAME}_f{i:03d}.png")
        im.convert("RGB").save(dst)
        out.append(dst)
    return out


if __name__ == "__main__":
    total, d = build()
    gif, mp4 = encode()
    n_intro = int(round(T_INTRO * FPS))
    n_draw = int(round(T_DRAW * FPS))
    picks = [2, n_intro + 12, n_intro + 40, n_intro + 80, n_intro + 120, total - 1]
    stills = export_stills(picks, gif)
    print("frames", total, "duration_s", total / FPS)
    for p in (gif, mp4):
        print(p, os.path.getsize(p))
    for s in stills:
        print(s)

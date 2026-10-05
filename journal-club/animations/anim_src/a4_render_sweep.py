"""Animation A4: GESFIDE curve (38 echoes) while SO2, R, CBV sweep one at a time.

Reads anim_src/_a4_sweep_cache.npz (run a4_sweep_data.py first).
Writes SCRATCH/assets/anim/anim_sweep.{gif,mp4} and representative frames.
"""
import os
import shutil
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = os.path.dirname(HERE)
OUT_DIR = os.path.join(SCRATCH, "assets", "anim")
FRAMES_OUT = os.path.join(OUT_DIR, "frames")
TMP = os.path.join(HERE, "_a4_frames")
NAME = "anim_sweep"

INK, MUTED, HAIR = "#1B1B1F", "#5B5F66", "#D9DADC"
RED, BLUE, LBLUE, LRED = "#A8322D", "#1F4E79", "#8DB3D9", "#E6B3AE"
REF_GRAY = "#B4B7BC"
DIM = "#A3A6AC"

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Liberation Serif"], "mathtext.fontset": "stix",
    "font.size": 16, "axes.edgecolor": MUTED, "axes.linewidth": 1.0, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
    "xtick.labelsize": 16, "ytick.labelsize": 16, "axes.labelsize": 18,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
})

W, H, DPI = 1100, 620, 100
FPS = 12
N_HOLD0, N_SWEEP, N_HOLD1, N_FINAL = 5, 36, 8, 24     # frames

d = np.load(os.path.join(HERE, "_a4_sweep_cache.npz"))
TE = d["TE"]


def _norm(S):
    return S / S[..., :1]


S_so2, so2_grid = _norm(d["S_so2"]), d["SO2_fine"]          # SO2 descending 0.90 -> 0.40
S_R, R_grid = _norm(d["S_R"]), d["R_sweep"]
S_cbv, cbv_grid = _norm(d["S_cbv"]), d["CBV_sweep"]


def _interp_rows(x, xs, rows):
    i = int(np.clip(np.searchsorted(xs, x) - 1, 0, len(xs) - 2))
    w = (x - xs[i]) / (xs[i + 1] - xs[i])
    return (1 - w) * rows[i] + w * rows[i + 1]


def curve_so2(v):
    return _interp_rows(v, so2_grid[::-1], S_so2[::-1])


def curve_R(v):
    return _interp_rows(np.log(v), np.log(R_grid), S_R)


def curve_cbv(v):
    return _interp_rows(v, cbv_grid, S_cbv)


def _fmt_R(v):
    return f"{v:.1f} µm" if v < 9.95 else f"{v:.0f} µm"


SEGMENTS = [
    dict(key="SO₂", value=lambda s: 0.90 - 0.50 * s, curve=curve_so2,
         label=lambda v: f"SO₂ = {100*v:.0f}%", ends=("90%", "40%"),
         fixed=["R = 5 µm", "CBV = 3%"], start="SO₂ = 90%"),
    dict(key="R", value=lambda s: 2.0 * 15.0 ** s, curve=curve_R,
         label=lambda v: f"R = {_fmt_R(v)}", ends=("2 µm", "30 µm"),
         fixed=["SO₂ = 60%", "CBV = 3%"], start="R = 2 µm"),
    dict(key="CBV", value=lambda s: 0.01 + 0.05 * s, curve=curve_cbv,
         label=lambda v: f"CBV = {100*v:.1f}%", ends=("1%", "6%"),
         fixed=["SO₂ = 60%", "R = 5 µm"], start="CBV = 1%"),
]


def ease(u):
    return u * u * (3 - 2 * u)


def frame_plan():
    plan = []
    for k, seg in enumerate(SEGMENTS):
        last = k == len(SEGMENTS) - 1
        plan += [(k, 0.0)] * N_HOLD0
        plan += [(k, ease(i / (N_SWEEP - 1))) for i in range(N_SWEEP)]
        plan += [(k, 1.0)] * (N_FINAL if last else N_HOLD1)
    return plan


def y_limits():
    curves = [seg["curve"](seg["value"](s)) for seg in SEGMENTS for s in np.linspace(0, 1, 41)]
    lo, hi = min(c.min() for c in curves), max(c.max() for c in curves)
    return np.floor((lo - 0.03) / 0.05) * 0.05, max(1.05, np.ceil((hi + 0.02) / 0.05) * 0.05)


def build_figure():
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    ax = fig.add_axes([0.088, 0.135, 0.585, 0.765])
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(0, 134)
    ax.set_xticks([0, 25, 50, 75, 100, 125])
    y0, y1 = y_limits()
    ax.set_ylim(y0, y1)
    ax.set_yticks(np.round(np.arange(np.ceil(y0 * 10) / 10, y1 + 1e-9, 0.1), 2))
    ax.set_xlabel("Echo time (ms)", labelpad=6)
    ax.set_ylabel("Signal (normalized)", labelpad=8)
    ax.tick_params(length=5, width=1.0)
    ax.grid(axis="y", color="#EEEFF1", lw=0.8)
    ax.set_axisbelow(True)
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    for x, lab in [(50, "180°"), (100, "spin echo")]:
        ax.axvline(x, color=MUTED, lw=1.0, ls=(0, (5, 4)), zorder=1)
        ax.text(x, 1.012, lab, transform=tr, ha="center", va="bottom", fontsize=16, color=MUTED)

    art = {}
    art["ref"], = ax.plot([], [], "-o", color=REF_GRAY, lw=1.6, ms=4.2, mec=REF_GRAY, zorder=2)
    art["cur"], = ax.plot([], [], "-o", color=BLUE, lw=1.9, ms=6.2, mec=BLUE, zorder=4)
    art["fill"] = None
    art["ax"] = ax

    # ---- right panel
    x0, x1 = 0.725, 0.955
    art["steps"] = []
    fig.text(x0, 0.905, "Sweep", fontsize=16, color=MUTED, va="baseline")
    for xs, seg in zip([0.806, 0.866, 0.899], SEGMENTS):
        art["steps"].append(fig.text(xs, 0.905, seg["key"], fontsize=16, color=DIM, va="baseline"))
    art["value"] = fig.text(x0, 0.775, "", fontsize=32, color=BLUE, va="baseline")
    ys = 0.705
    fig.add_artist(Line2D([x0 + 0.004, x1 - 0.004], [ys, ys], color=HAIR, lw=4, solid_capstyle="round"))
    art["prog"] = Line2D([x0 + 0.004, x0 + 0.004], [ys, ys], color=LBLUE, lw=4, solid_capstyle="round")
    fig.add_artist(art["prog"])
    art["knob"] = Line2D([x0], [ys], marker="o", ms=13, color=BLUE, mec="white", mew=1.5)
    fig.add_artist(art["knob"])
    art["end0"] = fig.text(x0, 0.645, "", fontsize=16, color=MUTED, ha="left", va="baseline")
    art["end1"] = fig.text(x1, 0.645, "", fontsize=16, color=MUTED, ha="right", va="baseline")
    art["slider"] = (x0 + 0.004, x1 - 0.004, ys)

    fig.text(x0, 0.535, "Held fixed", fontsize=16, color=MUTED, style="italic", va="baseline")
    art["fixed"] = [fig.text(x0, 0.470, "", fontsize=19, color=INK, va="baseline"),
                    fig.text(x0, 0.410, "", fontsize=19, color=INK, va="baseline")]

    fig.add_artist(Line2D([x0, x0 + 0.04], [0.316, 0.316], color=REF_GRAY, lw=1.6, marker="o",
                          ms=4.2, markevery=[1], mec=REF_GRAY))
    art["leg_ref"] = fig.text(x0 + 0.052, 0.305, "", fontsize=16, color=MUTED, va="baseline")
    fig.add_artist(Line2D([x0, x0 + 0.04], [0.256, 0.256], color=BLUE, lw=1.9, marker="o",
                          ms=6.2, markevery=[1], mec=BLUE))
    fig.text(x0 + 0.052, 0.245, "current", fontsize=16, color=INK, va="baseline")

    fig.text(x0, 0.175, "2D Monte Carlo, 38 echoes", fontsize=16, color=MUTED, va="baseline")
    fig.text(x0, 0.128, "EV + IV signal, no T2 decay", fontsize=16, color=MUTED, va="baseline")
    return fig, art


def draw_frame(fig, art, k, s):
    seg = SEGMENTS[k]
    ax = art["ax"]
    v = seg["value"](s)
    ref = seg["curve"](seg["value"](0.0))
    cur = seg["curve"](v)
    # no echoes around the 180 pulse: break the line between part A and part B
    gap = int(np.searchsorted(TE, 50.0))
    x = np.insert(TE, gap, np.nan)
    ref, cur = np.insert(ref, gap, np.nan), np.insert(cur, gap, np.nan)
    art["ref"].set_data(x, ref)
    art["cur"].set_data(x, cur)
    if art["fill"] is not None:
        art["fill"].remove()
    ok = ~np.isnan(x)
    art["fill"] = ax.fill_between(x, ref, cur, where=ok, interpolate=False, color=LBLUE, alpha=0.28,
                                  lw=0, zorder=1.5)
    art["value"].set_text(seg["label"](v))
    for j, t in enumerate(art["steps"]):
        t.set_color(INK if j == k else DIM)
        t.set_fontweight("bold" if j == k else "normal")
    a, b, ys = art["slider"]
    xk = a + (b - a) * s
    art["knob"].set_data([xk], [ys])
    art["prog"].set_data([a, xk], [ys, ys])
    art["end0"].set_text(seg["ends"][0])
    art["end1"].set_text(seg["ends"][1])
    for t, txt in zip(art["fixed"], seg["fixed"]):
        t.set_text(txt)
    art["leg_ref"].set_text("start: " + seg["start"])


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(FRAMES_OUT, exist_ok=True)
    shutil.rmtree(TMP, ignore_errors=True)
    os.makedirs(TMP)
    fig, art = build_figure()
    plan = frame_plan()
    for i, (k, s) in enumerate(plan):
        draw_frame(fig, art, k, s)
        fig.canvas.draw()
        img = np.asarray(fig.canvas.buffer_rgba())[..., :3]
        assert img.shape[:2] == (H, W), img.shape
        Image.fromarray(img).save(os.path.join(TMP, f"f_{i:04d}.png"))
    plt.close(fig)
    n = len(plan)
    print(f"{n} frames at {FPS} fps = {n/FPS:.1f} s")

    gif = os.path.join(OUT_DIR, NAME + ".gif")
    mp4 = os.path.join(OUT_DIR, NAME + ".mp4")
    pal = os.path.join(TMP, "palette.png")
    src = ["-framerate", str(FPS), "-i", os.path.join(TMP, "f_%04d.png")]
    subprocess.run(["ffmpeg", "-y", "-v", "error", *src, "-vf", "palettegen=stats_mode=diff:max_colors=128",
                    pal], check=True)
    subprocess.run(["ffmpeg", "-y", "-v", "error", *src, "-i", pal, "-lavfi",
                    "paletteuse=dither=bayer:bayer_scale=5:diff_mode=rectangle", "-loop", "0", gif], check=True)
    subprocess.run(["ffmpeg", "-y", "-v", "error", *src, "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20",
                    "-preset", "slow", "-movflags", "+faststart", mp4], check=True)
    for p in (gif, mp4):
        print(f"{p}: {os.path.getsize(p)/1e6:.2f} MB")

    # representative frames
    for f in os.listdir(FRAMES_OUT):
        if f.startswith(NAME + "_f"):
            os.remove(os.path.join(FRAMES_OUT, f))
    seg_len = N_HOLD0 + N_SWEEP + N_HOLD1
    picks = [N_HOLD0 + N_SWEEP // 2, seg_len - 1, seg_len + N_HOLD0 + N_SWEEP // 2, 2 * seg_len - 1,
             2 * seg_len + N_HOLD0 + N_SWEEP // 2, n - 1]
    for p in picks:
        shutil.copy(os.path.join(TMP, f"f_{p:04d}.png"), os.path.join(FRAMES_OUT, f"{NAME}_f{p:03d}.png"))
    print("frames:", picks)
    shutil.rmtree(TMP, ignore_errors=True)


if __name__ == "__main__":
    main()

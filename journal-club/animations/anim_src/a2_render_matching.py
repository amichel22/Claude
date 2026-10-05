"""Animation A2: dictionary matching of one 'measured' GESFIDE voxel (MRvF).

Inputs : SCRATCH/assets/anim/mrvf_dict.npz       (dictionary, T2-off total signal)
         anim_src/_a2_measured_cache.npz          (fresh-geometry voxel, run a2_measured_data.py)
Outputs: SCRATCH/assets/anim/anim_matching.{gif,mp4} + representative frames.

Matching is computed, never hard-coded: unit-norm measured vector x, unit-norm dictionary
atoms a_j, similarity_j = <x, a_j>, best = argmax_j.  The animation shows a random order of
candidates (first N_SHOW one by one, the rest fast-forwarded) and the running maximum.
Usage: python a2_render_matching.py [--snr 1000] [--name anim_matching]
"""
import argparse
import os
import shutil
import subprocess

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import blended_transform_factory
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = os.path.dirname(HERE)
OUT_DIR = os.path.join(SCRATCH, "assets", "anim")
FRAMES_OUT = os.path.join(OUT_DIR, "frames")

INK, MUTED, HAIR = "#1B1B1F", "#5B5F66", "#D9DADC"
RED, BLUE, LBLUE, LRED = "#A8322D", "#1F4E79", "#8DB3D9", "#E6B3AE"
CAND = "#A3A6AC"
UNVISITED = "#EEEFF1"
CMAP = LinearSegmentedColormap.from_list("sim", ["#D6E4F2", LBLUE, BLUE])

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Liberation Serif"], "mathtext.fontset": "stix",
    "font.size": 17, "axes.edgecolor": MUTED, "axes.linewidth": 1.0, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
    "xtick.labelsize": 16, "ytick.labelsize": 16, "axes.labelsize": 18,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
})

W, H, DPI = 1200, 620, 100
FPS = 12
N_INTRO, N_SLOW, F_SLOW, N_MED, F_MED, N_FAST, N_FINAL = 12, 6, 6, 26, 2, 14, 30
N_SHOW = N_SLOW + N_MED               # candidates shown one at a time (32)
Q_LO, Q_HI = 2.6, 5.4                 # colour range of -log10(1 - similarity)
N_SHADES = 8


# ------------------------------------------------------------------ data
def load(snr, noise_seed=0, order_seed=0):
    D = np.load(os.path.join(OUT_DIR, "mrvf_dict.npz"))
    M = np.load(os.path.join(HERE, "_a2_measured_cache.npz"))
    TE, SO2, R, CBV = D["TE"], D["SO2"], D["R"], D["CBV"]
    S = D["signals"]                                      # (nSO2, nR, nCBV, nTE), T2 off
    shape = S.shape[:3]
    A = S.reshape(-1, S.shape[-1])
    A = A / np.linalg.norm(A, axis=1, keepdims=True)      # unit-norm atoms
    idx = np.array(np.unravel_index(np.arange(len(A)), shape)).T
    Mt = M["M_total"]                                     # complex, noise-free
    sigma = abs(Mt[0]) / snr                              # SNR relative to the first echo
    rng = np.random.default_rng([noise_seed, 0])
    meas = np.abs(Mt + sigma * (rng.standard_normal(Mt.shape) + 1j * rng.standard_normal(Mt.shape)))
    x = meas / np.linalg.norm(meas)
    sims = A @ x
    best = int(np.argmax(sims))
    order = np.random.default_rng([order_seed, 1]).permutation(len(A))
    truth = (float(M["SO2_true"]), float(M["CBV_true"]), float(M["R_true"]))
    return dict(TE=TE, SO2=SO2, R=R, CBV=CBV, A=A, idx=idx, shape=shape, x=x, sims=sims,
                best=best, order=order, truth=truth, snr=snr, noise_free=np.abs(Mt))


# ------------------------------------------------------------------ plan
def frame_plan(nd):
    """List of (n_checked, shown_candidate or None, phase)."""
    plan = [(0, None, "intro")] * N_INTRO
    for i in range(N_SHOW):
        reps = F_SLOW if i < N_SLOW else F_MED
        plan += [(i + 1, i, "slow")] * reps
    rest = nd - N_SHOW
    for f in range(N_FAST):
        n = N_SHOW + int(round(rest * (f + 1) / N_FAST))
        plan.append((n, n - 1, "fast"))
    plan += [(nd, None, "final")] * N_FINAL
    return plan


# ------------------------------------------------------------------ figure
def gap_xy(TE, y):
    g = int(np.searchsorted(TE, 50.0))
    return np.insert(TE, g, np.nan), np.insert(y, g, np.nan)


def fmt_R(r):
    return f"{r:g} µm"


def build_figure(d):
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    TE = d["TE"]
    # ---------------- left: signals
    ax = fig.add_axes([0.075, 0.125, 0.535, 0.645])
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(0, 134)
    ax.set_xticks([0, 25, 50, 75, 100, 125])
    ax.set_ylim(0.085, 0.27)
    ax.set_yticks([0.10, 0.15, 0.20, 0.25])
    ax.set_xlabel("Echo time (ms)", labelpad=6)
    ax.set_ylabel("Signal / ‖signal‖", labelpad=8)
    ax.tick_params(length=5, width=1.0)
    ax.grid(axis="y", color="#EEEFF1", lw=0.8)
    ax.set_axisbelow(True)
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    for xv, lab in [(50, "180°"), (100, "spin echo")]:
        ax.axvline(xv, color=MUTED, lw=1.0, ls=(0, (5, 4)), zorder=1)
        ax.text(xv, 1.012, lab, transform=tr, ha="center", va="bottom", fontsize=16, color=MUTED)
    art = dict(ax=ax)
    art["fast"] = [ax.plot([], [], "-", color=CAND, lw=1.2, alpha=0.45, zorder=2)[0] for _ in range(8)]
    art["cand"], = ax.plot([], [], "-", color=CAND, lw=2.2, zorder=3)
    art["best"], = ax.plot([], [], "-", color=RED, lw=3.0, zorder=4)
    mx, my = gap_xy(TE, d["x"])
    ax.plot(mx, my, "o", color=INK, ms=5.6, mec="white", mew=0.5, zorder=5)
    # legend for the dots (lower left is free of curves)
    ax.plot([5.5], [0.0975], "o", color=INK, ms=5.6, mec="white", mew=0.5, zorder=5, clip_on=False)
    ax.text(9.0, 0.0975, f"measured voxel (SNR {d['snr']:g})", fontsize=16, color=INK, va="center")

    # ---------------- header table above the left axes
    cols = dict(SO2=0.300, R=0.370, CBV=0.455, sim=0.520)
    yh, y1, y2 = 0.945, 0.890, 0.835
    for key, lab in [("SO2", "SO₂"), ("R", "R"), ("CBV", "CBV"), ("sim", "similarity")]:
        fig.text(cols[key], yh, lab, fontsize=16, color=MUTED, va="baseline")
    art["sw_cand"] = Line2D([0.075, 0.108], [y1 + 0.012] * 2, color=CAND, lw=2.6)
    art["sw_best"] = Line2D([0.075, 0.108], [y2 + 0.012] * 2, color=RED, lw=2.6)
    fig.add_artist(art["sw_cand"])
    fig.add_artist(art["sw_best"])
    art["lab_cand"] = fig.text(0.118, y1, "candidate", fontsize=18, color=INK, va="baseline")
    art["lab_best"] = fig.text(0.118, y2, "best so far", fontsize=18, color=RED, va="baseline")
    art["row_cand"] = {k: fig.text(cols[k], y1, "", fontsize=18, color=INK, va="baseline") for k in cols}
    art["row_best"] = {k: fig.text(cols[k], y2, "", fontsize=18, color=RED, va="baseline") for k in cols}

    # ---------------- right: SO2 x R map of the dictionary nodes
    SO2p, R = 100 * d["SO2"], d["R"]
    hx = fig.add_axes([0.778, 0.435, 0.192, 0.42])
    hx.set_xscale("log")
    lr = np.log(R)
    mid = (lr[1:] + lr[:-1]) / 2
    redges = np.exp(np.concatenate([[lr[0] - (mid[0] - lr[0])], mid, [lr[-1] + (lr[-1] - mid[-1])]]))
    sedges = np.concatenate([SO2p - 2.5, [SO2p[-1] + 2.5]])
    hx.set_xlim(redges[0], redges[-1])
    hx.set_ylim(sedges[0], sedges[-1])
    hx.set_xticks([2, 3, 5, 10, 20, 30])
    hx.set_xticklabels(["2", "3", "5", "10", "20", "30"])
    hx.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    hx.set_yticks([30, 40, 50, 60, 70, 80, 90])
    hx.set_xlabel("R (µm, log scale)", labelpad=4, fontsize=17)
    hx.set_ylabel("SO₂ (%)", labelpad=4, fontsize=17)
    hx.tick_params(length=4, width=1.0)
    for s in ["top", "right", "left", "bottom"]:
        hx.spines[s].set_visible(False)
    cells = np.empty((len(SO2p), len(R)), dtype=object)
    for i in range(len(SO2p)):
        for j in range(len(R)):
            rect = Rectangle((redges[j], sedges[i]), redges[j + 1] - redges[j], sedges[i + 1] - sedges[i],
                             facecolor=UNVISITED, edgecolor="white", lw=1.4, zorder=1)
            hx.add_patch(rect)
            cells[i, j] = rect
    art["cells"] = cells
    art["cur_box"] = Rectangle((1, 1), 0, 0, fill=False, edgecolor=INK, lw=2.0, zorder=3, visible=False)
    hx.add_patch(art["cur_box"])
    art["best_ring"], = hx.plot([], [], "o", ms=20, mfc="none", mec=RED, mew=2.6, zorder=4, clip_on=False)
    art["hx"], art["redges"], art["sedges"] = hx, redges, sedges
    fig.text(0.735, 0.945, f"Dictionary: {len(d['A'])} entries", fontsize=18, color=INK, va="baseline")
    fig.text(0.735, 0.893, "darker = higher similarity", fontsize=16, color=MUTED, va="baseline")

    # ---------------- right bottom: progress during the search, readout at the end
    art["prog_txt"] = fig.text(0.735, 0.215, "", fontsize=16, color=MUTED, va="baseline")
    art["prog_bg"] = Line2D([0.737, 0.965], [0.175, 0.175], color=HAIR, lw=4, solid_capstyle="round")
    art["prog"] = Line2D([0.737, 0.737], [0.175, 0.175], color=LBLUE, lw=4, solid_capstyle="round")
    fig.add_artist(art["prog_bg"])
    fig.add_artist(art["prog"])
    rc = dict(name=0.735, best=0.805, true=0.905)
    ry = [0.245, 0.180, 0.120, 0.060]
    ro = []
    ro.append(fig.text(rc["best"], ry[0], "best match", fontsize=17, color=RED, va="baseline"))
    ro.append(fig.text(rc["true"], ry[0], "true", fontsize=17, color=MUTED, va="baseline"))
    art["ro_vals"] = []
    for k, lab in enumerate(["SO₂", "CBV", "R"]):
        ro.append(fig.text(rc["name"], ry[k + 1], lab, fontsize=20, color=INK, va="baseline"))
        vb = fig.text(rc["best"], ry[k + 1], "", fontsize=20, color=INK, va="baseline")
        vt = fig.text(rc["true"], ry[k + 1], "", fontsize=20, color=MUTED, va="baseline")
        ro += [vb, vt]
        art["ro_vals"].append((vb, vt))
    art["readout"] = ro
    return fig, art


def params_of(d, j):
    i_s, i_r, i_c = d["idx"][j]
    return d["SO2"][i_s], d["R"][i_r], d["CBV"][i_c]


def set_row(row, d, j, sim_val=None):
    if j is None:
        for t in row.values():
            t.set_text("–" if t is not row["sim"] else "–")
        return
    so2, r, cbv = params_of(d, j)
    row["SO2"].set_text(f"{100*so2:.0f}%")
    row["R"].set_text(fmt_R(r))
    row["CBV"].set_text(f"{100*cbv:.0f}%")
    row["sim"].set_text(f"{d['sims'][j]:.6f}")


def cell_color(s):
    # 8 discrete shades (clean GIF palette); q = -log10(1 - similarity)
    q = -np.log10(max(1.0 - s, 1e-12))
    k = int(np.clip(np.floor((q - Q_LO) / (Q_HI - Q_LO) * N_SHADES), 0, N_SHADES - 1))
    return CMAP(k / (N_SHADES - 1))


def draw_frame(fig, art, d, n_checked, shown, phase, state):
    TE, A, order, sims = d["TE"], d["A"], d["order"], d["sims"]
    checked = order[:n_checked]
    # ---- running best (computed from the similarities of the entries checked so far)
    best_j = int(checked[np.argmax(sims[checked])]) if n_checked else None
    final = phase == "final"
    # ---- heat map: visited cells coloured by the best similarity over CBV seen so far
    new = order[state["n_prev"]:n_checked]
    for j in new:
        i_s, i_r, _ = d["idx"][j]
        state["cellmax"][i_s, i_r] = max(state["cellmax"][i_s, i_r], sims[j])
        art["cells"][i_s, i_r].set_facecolor(cell_color(state["cellmax"][i_s, i_r]))
    state["n_prev"] = n_checked
    # current-candidate box
    cur = order[shown] if shown is not None else None
    if cur is not None and phase == "slow":
        i_s, i_r, _ = d["idx"][cur]
        re, se = art["redges"], art["sedges"]
        art["cur_box"].set_bounds(re[i_r], se[i_s], re[i_r + 1] - re[i_r], se[i_s + 1] - se[i_s])
        art["cur_box"].set_visible(True)
    else:
        art["cur_box"].set_visible(False)
    if best_j is not None:
        so2, r, _ = params_of(d, best_j)
        art["best_ring"].set_data([r], [100 * so2])
        art["best_ring"].set_markeredgewidth(3.2 if final else 2.4)
    else:
        art["best_ring"].set_data([], [])
    # ---- curves
    for ln in art["fast"]:
        ln.set_data([], [])
    if phase == "slow":
        art["cand"].set_data(*gap_xy(TE, A[cur]))
    elif phase == "fast":
        art["cand"].set_data(*gap_xy(TE, A[cur]))
        batch = order[max(N_SHOW, n_checked - 52):n_checked]
        pick = batch[np.linspace(0, len(batch) - 1, len(art["fast"])).astype(int)]
        for ln, j in zip(art["fast"], pick):
            ln.set_data(*gap_xy(TE, A[j]))
    else:
        art["cand"].set_data([], [])
    if best_j is not None:
        art["best"].set_data(*gap_xy(TE, A[best_j]))
    else:
        art["best"].set_data([], [])
    # ---- header table
    if final:
        art["sw_cand"].set_visible(False)
        art["lab_cand"].set_text(f"all {len(A)} entries checked")
        art["lab_cand"].set_color(MUTED)
        for t in art["row_cand"].values():
            t.set_text("")
        art["lab_best"].set_text("best match")
    else:
        art["sw_cand"].set_visible(True)
        art["lab_cand"].set_text("candidate")
        art["lab_cand"].set_color(INK)
        set_row(art["row_cand"], d, cur)
        art["lab_best"].set_text("best so far")
    set_row(art["row_best"], d, best_j)
    # ---- progress / readout
    if final:
        art["prog_txt"].set_text("")
        art["prog_bg"].set_visible(False)
        art["prog"].set_visible(False)
        so2, r, cbv = params_of(d, best_j)
        ts, tc, tr_ = d["truth"]
        vals = [(f"{100*so2:.0f}%", f"{100*ts:.0f}%"), (f"{100*cbv:.0f}%", f"{100*tc:.0f}%"),
                (fmt_R(r), fmt_R(tr_))]
        for (vb, vt), (a, b) in zip(art["ro_vals"], vals):
            vb.set_text(a)
            vt.set_text(b)
        for t in art["readout"]:
            t.set_visible(True)
    else:
        art["prog_txt"].set_text(f"checked {n_checked} of {len(A)} entries")
        art["prog_bg"].set_visible(True)
        art["prog"].set_visible(True)
        a0, a1 = 0.737, 0.965
        art["prog"].set_data([a0, a0 + (a1 - a0) * n_checked / len(A)], [0.175, 0.175])
        for t in art["readout"]:
            t.set_visible(False)
    return best_j


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snr", type=float, default=1000.0)
    ap.add_argument("--name", default="anim_matching")
    args = ap.parse_args()
    name = args.name
    tmp = os.path.join(HERE, f"_{name}_frames")
    d = load(args.snr)
    bj = d["best"]
    so2, r, cbv = params_of(d, bj)
    print(f"SNR {args.snr:g}: best match SO2 {so2:.2f}, R {r:g} um, CBV {cbv:.2f}, similarity {d['sims'][bj]:.7f}")
    top = np.argsort(-d["sims"])[:5]
    for j in top:
        print("   ", tuple(float(v) for v in params_of(d, j)), f"{d['sims'][j]:.7f}")
    pos = int(np.flatnonzero(d["order"] == bj)[0])
    print(f"global best is entry #{pos+1} in the display order")

    os.makedirs(FRAMES_OUT, exist_ok=True)
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    fig, art = build_figure(d)
    plan = frame_plan(len(d["A"]))
    state = dict(n_prev=0, cellmax=np.full(d["shape"][:2], -np.inf))
    for i, (n, shown, phase) in enumerate(plan):
        draw_frame(fig, art, d, n, shown, phase, state)
        fig.canvas.draw()
        img = np.asarray(fig.canvas.buffer_rgba())[..., :3]
        assert img.shape[:2] == (H, W), img.shape
        Image.fromarray(img).save(os.path.join(tmp, f"f_{i:04d}.png"))
    plt.close(fig)
    n = len(plan)
    print(f"{n} frames at {FPS} fps = {n/FPS:.1f} s")

    gif = os.path.join(OUT_DIR, name + ".gif")
    mp4 = os.path.join(OUT_DIR, name + ".mp4")
    pal = os.path.join(tmp, "palette.png")
    src = ["-framerate", str(FPS), "-i", os.path.join(tmp, "f_%04d.png")]
    subprocess.run(["ffmpeg", "-y", "-v", "error", *src, "-vf", "palettegen=stats_mode=diff:max_colors=256",
                    pal], check=True)
    subprocess.run(["ffmpeg", "-y", "-v", "error", *src, "-i", pal, "-lavfi",
                    "paletteuse=dither=bayer:bayer_scale=5:diff_mode=rectangle", "-loop", "0", gif], check=True)
    subprocess.run(["ffmpeg", "-y", "-v", "error", *src, "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20",
                    "-preset", "slow", "-movflags", "+faststart", mp4], check=True)
    for p in (gif, mp4):
        print(f"{p}: {os.path.getsize(p)/1e6:.2f} MB")

    for f in os.listdir(FRAMES_OUT):
        if f.startswith(name + "_f"):
            os.remove(os.path.join(FRAMES_OUT, f))
    s0 = N_INTRO
    picks = [N_INTRO - 1, s0 + 2 * F_SLOW + 3, s0 + N_SLOW * F_SLOW + 20 * F_MED,
             s0 + N_SLOW * F_SLOW + N_MED * F_MED + N_FAST // 2, n - 1]
    for p in picks:
        shutil.copy(os.path.join(tmp, f"f_{p:04d}.png"), os.path.join(FRAMES_OUT, f"{name}_f{p:03d}.png"))
    print("frames:", picks)
    shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()

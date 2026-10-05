"""Animation A3 (anim_spins): extravascular spin dephasing around a small vs a large vessel.

Physics in a3_spins_sim.py.  This script: run / cache the two simulations, render frames
with matplotlib at 1200x680, encode GIF (two-pass palette) + MP4 with ffmpeg.
"""
import os
import shutil
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from a3_spins_sim import simulate, dw_wall, dchi_si  # noqa: E402

SCRATCH = os.path.dirname(HERE)
OUT = os.path.join(SCRATCH, "assets", "anim")
FR_DIR = os.path.join(HERE, "_a3_frames")
CACHE = os.path.join(HERE, "_a3_sim_cache.npz")
NAME = "anim_spins"

# ---------------- parameters ----------------
SO2 = 0.60          # spec default (Hct 0.4 -> dchi_SI = 5.43e-7)
D = 1.0             # um^2/ms
DT = 0.05           # ms
N_SPINS = 40000     # per vessel, for the signal trace
N_SHOW = 300        # drawn
FOV = 5.0           # half-width of the periodic square, in units of R
RADII = (3.0, 30.0)
FPS = 12
DOT_AA = os.environ.get("A3_DOT_AA", "0") == "1"
DOT_S = float(os.environ.get("A3_DOT_S", "32"))
W_PX, H_PX, DPI = 1200, 680, 100

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
    "text.color": INK,
    "axes.edgecolor": "#3A3D42",
    "axes.labelcolor": INK,
    "xtick.color": "#3A3D42",
    "ytick.color": "#3A3D42",
    "xtick.labelcolor": INK,
    "ytick.labelcolor": INK,
    "axes.linewidth": 1.0,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})


# ---------------- simulation (cached) ----------------
def get_sims(force=False):
    if os.path.exists(CACHE) and not force:
        z = np.load(CACHE, allow_pickle=True)
        return [z["small"].item(), z["large"].item()]
    sims = []
    for i, R in enumerate(RADII):
        s = simulate(R, so2=SO2, D=D, dt=DT, n=N_SPINS, fov_factor=FOV,
                     n_keep=N_SHOW, seed=11 + i)
        s["m_ge"] = np.abs(s["m_ge"])
        s["m_se"] = np.abs(s["m_se"])
        sims.append(s)
    np.savez_compressed(CACHE, small=np.array(sims[0], dtype=object),
                        large=np.array(sims[1], dtype=object))
    return sims


# ---------------- frame schedule ----------------
def frame_schedule():
    """List of (step_index, flipped_flag, tag).  dt grid = 0.05 ms."""
    def half(t0, t1, slow_at):
        # speed (ms/frame): 0.8 normally, dipping to ~0.22 near the 180 deg pulse
        ts = [t0]
        t = t0
        while t < t1 - 1e-9:
            v = 1.00 - 0.72 * np.exp(-((t - slow_at) / 5.0) ** 2)
            t = t + v
            ts.append(min(t, t1))
        ts = np.array(ts)
        ts = t0 + (ts - t0) * (t1 - t0) / (ts[-1] - t0)
        return ts
    a = half(0.0, 50.0, 50.0)
    b = half(50.0, 100.0, 50.0)
    sched = []
    for _ in range(7):
        sched.append((0, False, "start"))
    for t in a[1:]:
        sched.append((int(round(t / DT)), False, "run"))
    for _ in range(6):
        sched.append((int(round(50 / DT)), False, "preflip"))
    for _ in range(9):
        sched.append((int(round(50 / DT)), True, "flip"))
    for t in b[1:]:
        sched.append((int(round(t / DT)), True, "run2"))
    return sched


# ---------------- drawing helpers ----------------
def ax_px(fig, x, y_top, w, h, **kw):
    """Axes from pixel box (origin top-left)."""
    return fig.add_axes([x / W_PX, 1 - (y_top + h) / H_PX, w / W_PX, h / H_PX], **kw)


def field_image(n=362):
    u = np.linspace(-FOV, FOV, n)          # in units of R
    X, Y = np.meshgrid(u, u)
    r2 = X ** 2 + Y ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        f = (Y ** 2 - X ** 2) / r2 ** 2     # (R/r)^2 cos(2 phi), B0 along +y
    f[r2 < 1] = 0
    g = np.sign(f) * np.sqrt(np.clip(np.abs(f), 0, 1))   # display compression only
    g = np.round(g * 5) / 5                               # 11 flat bands (GIF-friendly)
    cmap = plt.get_cmap("RdBu_r")
    rgb = cmap(0.5 + 0.5 * g)[..., :3]
    a = 0.42
    rgb = (1 - a) * 1.0 + a * rgb
    return rgb


def phase_cmap():
    """Cyclic phase map: in phase (0) = near-black; + phase -> red, - phase -> blue;
    +-pi = pale sand.  Lightness rises quickly with |phase| so small phase spreads show."""
    anchors = [(0.00, "#F3E6B5"), (0.17, "#A9C9E8"), (0.31, "#3D7DC8"), (0.42, "#1F4E79"),
               (0.50, "#1B1B1F"), (0.58, "#8A2A25"), (0.69, "#D0473D"), (0.83, "#EFAA9B"),
               (1.00, "#F3E6B5")]
    return LinearSegmentedColormap.from_list("phase_cyc", anchors, N=36)


def wrap(p):
    return (p + np.pi) % (2 * np.pi) - np.pi


def build_figure(sims):
    fig = plt.figure(figsize=(W_PX / DPI, H_PX / DPI), dpi=DPI)
    cm = phase_cmap()
    fimg = field_image()
    PANEL = 362
    px = [64, 774]
    py = 52
    tx = [104, 814]
    tw = PANEL - 40 + 0
    ty, th = 482, 128

    H = {"panels": [], "traces": []}
    titles = ["Small vessel, R = 3 µm", "Large vessel, R = 30 µm"]
    scale = [(10.0, "10 µm"), (100.0, "100 µm")]
    for j, s in enumerate(sims):
        R, L = s["R"], s["L"]
        ax = ax_px(fig, px[j], py, PANEL, PANEL)
        ax.imshow(fimg, extent=[-L, L, -L, L], origin="lower", interpolation="nearest",
                  zorder=0)
        ax.add_patch(Circle((0, 0), R, facecolor=RED, edgecolor="none", zorder=2))
        ax.set_xlim(-L, L)
        ax.set_ylim(-L, L)
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color(HAIR)
            sp.set_linewidth(1.2)
        p0 = s["keep_pos"][0]
        sc = ax.scatter(p0[:, 0], p0[:, 1], s=DOT_S, c=np.zeros(len(p0)), cmap=cm,
                        vmin=-np.pi, vmax=np.pi, edgecolors=INK, linewidths=0.4, zorder=3,
                        antialiased=DOT_AA)
        fig.text((px[j] + PANEL / 2) / W_PX, 1 - 30 / H_PX, titles[j], ha="center",
                 va="center", fontsize=20, color=INK)
        # scale bar under panel (outside)
        sl, slab = scale[j]
        bar_px = PANEL * sl / (2 * L)
        yb = 1 - (py + PANEL + 14) / H_PX
        x0 = px[j] / W_PX
        fig.add_artist(plt.Line2D([x0, x0 + bar_px / W_PX], [yb, yb], color=INK, lw=2.4,
                                  transform=fig.transFigure, solid_capstyle="butt"))
        fig.text(x0 + bar_px / W_PX + 8 / W_PX, yb, slab, ha="left", va="center",
                 fontsize=15, color=INK)
        H["panels"].append(dict(ax=ax, sc=sc))

        # trace
        at = ax_px(fig, tx[j], ty, tw, th)
        at.set_xlim(0, 100)
        at.set_ylim(0.4, 1.02)
        at.set_xticks([0, 25, 50, 75, 100])
        at.set_yticks([0.4, 0.6, 0.8, 1.0])
        at.set_yticklabels(["0.4", "0.6", "0.8", "1.0"])
        at.tick_params(labelsize=15, length=4)
        at.spines["top"].set_visible(False)
        at.spines["right"].set_visible(False)
        at.set_xlabel("Time (ms)", fontsize=16, labelpad=2)
        at.grid(axis="y", color="#EEEEF0", lw=0.8)
        at.set_axisbelow(True)
        v50 = at.axvline(50, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)
        at.axvline(100, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)
        lab50 = at.text(50, 1.035, "180°", ha="center", va="bottom", fontsize=15,
                        color=MUTED, transform=at.transData, clip_on=False)
        at.text(100, 1.035, "spin echo", ha="right", va="bottom", fontsize=15,
                color=MUTED, clip_on=False)
        at.text(0.0, 1.035, "Signal", ha="left", va="bottom", fontsize=15,
                color=INK, clip_on=False)
        (ge_line,) = at.plot([], [], color="#9A9DA3", lw=1.6, ls=(0, (1.5, 2.0)), zorder=2)
        (se_line,) = at.plot([], [], color=BLUE, lw=2.4, zorder=3, solid_capstyle="round",
                              clip_on=False)
        (dot,) = at.plot([], [], "o", ms=7, color=BLUE, zorder=4, clip_on=False)
        ge_lab = at.text(0, 0, "", fontsize=14, color=MUTED, ha="right", va="top")
        end_lab = at.text(0, 0, "", fontsize=16, color=BLUE, ha="left", va="center",
                          clip_on=False)
        H["traces"].append(dict(ax=at, ge=ge_line, se=se_line, dot=dot, v50=v50,
                                lab50=lab50, ge_lab=ge_lab, end_lab=end_lab))

    # ---------- centre column ----------
    cx = 600 / W_PX
    H["time"] = fig.text(cx, 1 - 78 / H_PX, "", ha="center", va="center", fontsize=28,
                         color=INK)
    H["event"] = fig.text(cx, 1 - 122 / H_PX, "", ha="center", va="center", fontsize=17,
                          color=MUTED)
    # B0 arrow
    arr = FancyArrowPatch((cx, 1 - 232 / H_PX), (cx, 1 - 168 / H_PX),
                          transform=fig.transFigure, arrowstyle="-|>", mutation_scale=18,
                          color=INK, lw=1.8)
    fig.add_artist(arr)
    fig.text(cx + 14 / W_PX, 1 - 196 / H_PX, r"$B_0$", ha="left", va="center",
             fontsize=20, color=INK)
    fig.text(cx, 1 - 250 / H_PX, r"vessel axis $\perp$ $B_0$ (into screen)",
             ha="center", va="center", fontsize=14, color=MUTED)

    # phase colour strip
    cb_w, cb_h = 200, 16
    cbx = 600 - cb_w / 2
    fig.text(cx, 1 - 286 / H_PX, "Spin phase", ha="center", va="center", fontsize=16,
             color=INK)
    acb = ax_px(fig, cbx, 300, cb_w, cb_h)
    grad = np.linspace(-np.pi, np.pi, 256)[None, :]
    acb.imshow(grad, aspect="auto", cmap=cm, vmin=-np.pi, vmax=np.pi,
               extent=[-np.pi, np.pi, 0, 1])
    acb.set_yticks([])
    acb.set_xticks([-np.pi, 0, np.pi])
    acb.set_xticklabels([r"$-\pi$", "0", r"$+\pi$"], fontsize=16)
    acb.tick_params(length=3)
    for sp in acb.spines.values():
        sp.set_color(MUTED)
        sp.set_linewidth(0.8)

    # field strip
    fig.text(cx, 1 - 366 / H_PX, "Field offset (background)", ha="center", va="center",
             fontsize=16, color=INK)
    afb = ax_px(fig, cbx, 380, cb_w, cb_h)
    g = np.clip(np.round(np.linspace(-1.1, 1.1, 200) * 5) / 5, -1, 1)
    rgb = plt.get_cmap("RdBu_r")(0.5 + 0.5 * g)[None, :, :3]
    rgb = 0.58 + 0.42 * rgb
    afb.imshow(rgb, aspect="auto", extent=[-1, 1, 0, 1])
    afb.set_yticks([])
    afb.set_xticks([-1, 0, 1])
    afb.set_xticklabels(["−", "0", "+"], fontsize=16)
    afb.tick_params(length=3)
    for sp in afb.spines.values():
        sp.set_color(MUTED)
        sp.set_linewidth(0.8)

    # trace legend
    lx0 = 600 - 92
    ly = [512, 542]
    for k, (ls, col, lw, lab) in enumerate([("-", BLUE, 2.4, "with 180° (spin echo)"),
                                           ((0, (1.5, 2.0)), "#9A9DA3", 1.6, "no 180°")]):
        y = 1 - ly[k] / H_PX
        fig.add_artist(plt.Line2D([lx0 / W_PX, (lx0 + 34) / W_PX], [y, y], color=col, lw=lw,
                                  ls=ls, transform=fig.transFigure))
        fig.text((lx0 + 42) / W_PX, y, lab, ha="left", va="center", fontsize=15,
                 color=INK if k == 0 else MUTED)
    fig.text(cx, 1 - 470 / H_PX, r"Signal = $|\langle e^{i\varphi}\rangle|$", ha="center",
             va="center", fontsize=16, color=INK)

    # parameter note
    dchi = dchi_si(SO2)
    fig.text(cx, 1 - 604 / H_PX,
             f"SO₂ = {SO2*100:.0f}%, Hct = 0.4, B₀ = 3 T",
             ha="center", va="center", fontsize=14, color=MUTED)
    fig.text(cx, 1 - 628 / H_PX,
             "D = 1 µm²/ms, vessel = 3.1% of area",
             ha="center", va="center", fontsize=14, color=MUTED)
    fig.text(cx, 1 - 652 / H_PX,
             f"{N_SPINS//1000}k spins simulated, 300 shown",
             ha="center", va="center", fontsize=14, color=MUTED)
    return fig, H


def update(fig, H, sims, k, flipped, tag):
    t = k * DT
    for j, s in enumerate(sims):
        pos = s["keep_pos"][k]
        ph = s["keep_phi"][k]
        if flipped:
            ph = ph - 2 * s["phi180_keep"]
        P = H["panels"][j]
        P["sc"].set_offsets(pos)
        P["sc"].set_array(wrap(ph))
        flash = tag in ("preflip", "flip")
        for sp in P["ax"].spines.values():
            sp.set_color(RED if flash else HAIR)
            sp.set_linewidth(2.4 if flash else 1.2)

        T = H["traces"][j]
        tt = s["t"][: k + 1]
        if flipped:
            se = s["m_se"][: k + 1]
        else:
            se = s["m_ge"][: k + 1]
        T["se"].set_data(tt, se)
        T["dot"].set_data([t], [se[-1]])
        if flipped and k > s["i180"]:
            i0 = s["i180"]
            T["ge"].set_data(s["t"][i0: k + 1], s["m_ge"][i0: k + 1])
        else:
            T["ge"].set_data([], [])
        T["v50"].set_color(RED if flash else MUTED)
        T["v50"].set_linewidth(1.8 if flash else 1.0)
        T["lab50"].set_color(RED if flash else MUTED)
        if k == len(s["t"]) - 1:
            T["end_lab"].set_position((103.5, se[-1]))
            T["end_lab"].set_text(f"{se[-1]:.2f}")
        else:
            T["end_lab"].set_text("")
            T["ge_lab"].set_text("")

    H["time"].set_text(f"t = {t:.0f} ms")
    if tag == "start":
        H["event"].set_text("90° pulse: all spins in phase")
        H["event"].set_color(INK)
    elif tag == "preflip":
        H["event"].set_text("180° pulse ...")
        H["event"].set_color(RED)
    elif tag == "flip":
        H["event"].set_text("180° pulse: every phase flips sign")
        H["event"].set_color(RED)
    elif tag == "run":
        H["event"].set_text("spins dephase in the field")
        H["event"].set_color(MUTED)
    elif tag == "run2":
        if k >= len(sims[0]["t"]) - 1:
            H["event"].set_text("spin echo at 100 ms")
            H["event"].set_color(INK)
        else:
            H["event"].set_text("spins that stay put rephase")
            H["event"].set_color(MUTED)


def render(sims):
    if os.path.isdir(FR_DIR):
        shutil.rmtree(FR_DIR)
    os.makedirs(FR_DIR)
    fig, H = build_figure(sims)
    sched = frame_schedule()
    for i, (k, fl, tag) in enumerate(sched):
        update(fig, H, sims, k, fl, tag)
        fig.savefig(os.path.join(FR_DIR, f"{i:04d}.png"), dpi=DPI)
    plt.close(fig)
    return len(sched)


def encode(nframes, hold_s=2.0):
    os.makedirs(OUT, exist_ok=True)
    gif = os.path.join(OUT, NAME + ".gif")
    mp4 = os.path.join(OUT, NAME + ".mp4")
    pal = os.path.join(HERE, "_a3_palette.png")
    inp = ["-framerate", str(FPS), "-i", os.path.join(FR_DIR, "%04d.png")]
    hold = f"tpad=stop_mode=clone:stop_duration={hold_s}"
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", *inp, "-vf",
                    f"{hold},palettegen=max_colors=128:stats_mode=diff", pal], check=True)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", *inp, "-i", pal, "-lavfi",
                    f"[0:v]{hold}[v];[v][1:v]paletteuse=dither=bayer:bayer_scale=5:diff_mode=rectangle",
                    "-loop", "0", gif], check=True)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", *inp, "-vf", hold,
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20",
                    "-preset", "slow", "-movflags", "+faststart", mp4], check=True)
    return gif, mp4


if __name__ == "__main__":
    force = "--resim" in sys.argv
    sims = get_sims(force=force)
    for s in sims:
        print(f"R={s['R']:g}: GE50={s['m_ge'][s['i180']]:.3f} GE100={s['m_ge'][-1]:.3f} "
              f"SE100={s['m_se'][-1]:.3f} dw0={s['dw0']:.4f} rad/ms stuck={s['n_stuck']}")
    if "--sim-only" in sys.argv:
        sys.exit()
    n = render(sims)
    print("frames", n)
    if "--no-encode" not in sys.argv:
        gif, mp4 = encode(n)
        for f in (gif, mp4):
            print(f, os.path.getsize(f) / 1e6, "MB")

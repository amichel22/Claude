"""IR-FISP MRF signal model (extended phase graph) for animation A1.

Sequence (after Jiang et al., MRM 2015, FISP-MRF):
  * non-selective 180 deg inversion, perfect, followed by TI before the first pulse
  * N_TR RF pulses (phase 0) with smoothly varying pseudo-random flip angles
  * readout at constant TE; unbalanced (FISP) gradient dephases by one
    EPG order per TR; repetition times follow a smooth pseudo-random pattern

The signal at each TR is F+_0 at TE.  Returned as complex; plot |S|.
Relaxation only (no B1 / slice profile / off-resonance), M0 = 1.
"""
import numpy as np

N_TR = 500
TI_MS = 40.0
TE_MS = 2.0
SEED = 7

TISSUES = [
    # name,          T1 (ms), T2 (ms), color
    ("White matter", 830.0, 80.0, "#1F4E79"),
    ("Gray matter", 1330.0, 110.0, "#A8322D"),
    ("CSF", 4000.0, 2000.0, "#5B5F66"),
]


def _smooth_noise(n, n_knots, rng):
    """Smooth pseudo-random curve in [-1, 1] (cubic interpolation of random knots)."""
    from scipy.interpolate import CubicSpline
    xk = np.linspace(0, n - 1, n_knots)
    yk = rng.uniform(-1, 1, n_knots)
    y = CubicSpline(xk, yk, bc_type="natural")(np.arange(n))
    return y / max(1e-9, np.abs(y).max())


# Flip-angle lobes: (length in TRs, peak flip angle in deg, fraction of the half-sine used).
# The last lobe is cut near its peak so the three curves end well separated.
LOBES = [(136, 24.0, 1.0), (137, 63.0, 1.0), (158, 70.0, 1.0), (69, 22.0, 0.6)]


def make_schedule(n=N_TR, seed=SEED):
    """Flip angles (deg) and TRs (ms).

    Flip angles: consecutive half-sine lobes of unequal length and unequal
    peak (Jiang 2015 style), plus a small smooth perturbation, clipped 5-70 deg.
    TR: smooth (Perlin-like) pseudo-random pattern around 13 ms (11-15 ms)
    with a small per-TR jitter.
    """
    rng = np.random.default_rng(seed)
    assert sum(L for L, _, _ in LOBES) == n
    fa = []
    for L, pk, fr in LOBES:
        k = np.arange(L)
        fa.append(5.0 + (pk - 5.0) * np.sin(fr * np.pi * (k + 0.5) / L) ** 1.3)
    fa = np.concatenate(fa)
    # smooth perturbation, faded in smoothly above 5 deg so the train stays continuous
    env = np.clip((fa - 5.0) / 12.0, 0.0, 1.0)
    fa = fa + 3.0 * _smooth_noise(n, 40, rng) * env
    fa = np.clip(fa, 5.0, 70.0)

    tr = 13.0 + 1.5 * _smooth_noise(n, 28, rng)
    tr = tr + 0.15 * rng.standard_normal(n)  # small jitter on top of the smooth pattern
    tr = np.clip(tr, 11.0, 15.0)
    return fa, tr


def epg_ir_fisp(fa_deg, tr_ms, T1, T2, te_ms=TE_MS, ti_ms=TI_MS):
    n = len(fa_deg)
    K = n + 2
    Fp = np.zeros(K, complex)
    Fm = np.zeros(K, complex)
    Z = np.zeros(K, complex)

    # inversion + spoiled delay TI
    Z[0] = -1.0
    E1 = np.exp(-ti_ms / T1)
    Z[0] = Z[0] * E1 + (1 - E1)

    sig = np.zeros(n, complex)
    for i in range(n):
        a = np.deg2rad(fa_deg[i])
        c2, s2, sa, ca = np.cos(a / 2) ** 2, np.sin(a / 2) ** 2, np.sin(a), np.cos(a)
        # RF about x (phase 0)
        Fp_n = c2 * Fp + s2 * Fm - 1j * sa * Z
        Fm_n = s2 * Fp + c2 * Fm + 1j * sa * Z
        Z_n = -0.5j * sa * Fp + 0.5j * sa * Fm + ca * Z
        Fp, Fm, Z = Fp_n, Fm_n, Z_n

        # relax to TE, sample
        E1, E2 = np.exp(-te_ms / T1), np.exp(-te_ms / T2)
        Fp *= E2; Fm *= E2; Z *= E1; Z[0] += 1 - E1
        sig[i] = Fp[0]

        # relax rest of TR
        t = tr_ms[i] - te_ms
        E1, E2 = np.exp(-t / T1), np.exp(-t / T2)
        Fp *= E2; Fm *= E2; Z *= E1; Z[0] += 1 - E1

        # unbalanced gradient: shift by one order
        Fp = np.concatenate(([0.0], Fp[:-1]))
        Fm = np.concatenate((Fm[1:], [0.0]))
        Fp[0] = np.conj(Fm[0])
    return sig


def simulate():
    fa, tr = make_schedule()
    out = {"fa": fa, "tr": tr, "curves": []}
    for name, T1, T2, col in TISSUES:
        s = epg_ir_fisp(fa, tr, T1, T2)
        out["curves"].append((name, T1, T2, col, np.abs(s)))
    return out


if __name__ == "__main__":
    d = simulate()
    print("FA range", d["fa"].min(), d["fa"].max(), "TR range", d["tr"].min(), d["tr"].max())
    for name, T1, T2, col, s in d["curves"]:
        print(f"{name:14s} max {s.max():.3f} mean {s.mean():.3f} end {s[-1]:.3f}")

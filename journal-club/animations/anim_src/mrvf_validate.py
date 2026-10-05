"""Validation of mrvf_sim.py (prints a report).

(0) split far-grid + near-exact field vs brute-force analytic sum
(a) diffusion off, T2 off: signal refocuses to 1 at the 100 ms spin echo
(b) GE decay rate grows with CBV and with (1 - SO2)            [from the dictionary]
(c) SE attenuation vs radius peaks at a few um                [dictionary + extra small R]
(d) long-time EV R2' vs static-dephasing formula (rough sanity check only)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mrvf_sim as ms  # noqa: E402

SCRATCH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DICT = os.path.join(SCRATCH, "assets", "anim", "mrvf_dict.npz")
np.set_printoptions(precision=4, suppress=True, linewidth=150)
TE = ms.TE_MS
A = ms.PART_A
SE = ms.SE_INDEX


def main():
    t0 = time.time()
    V = {s: ms.nested_voxels([0.01, 0.03, 0.06], seed=s) for s in (0, 1)}
    g = V[0][0.03]

    print(f"(0) field accuracy, CBV 3 pct voxel ({g.n} vessels, L = {g.L:.1f} R), units of gamma*B0*dchi(SO2=0):")
    rng = np.random.default_rng(11)
    v = rng.integers(0, g.n, 20000)
    r = 1 + 4 * rng.random(20000)
    a = 2 * np.pi * rng.random(20000)
    for name, P in [("uniform points", rng.random((20000, 2)) * g.L),
                    ("1R < r < 5R", (g.xy[v] + np.c_[r * np.cos(a), r * np.sin(a)]) % g.L)]:
        f, _ = g.evaluate(P)
        fd = g.field_direct(P)
        print(f"    {name:15s}: max |err| = {np.abs(f-fd).max():.1e}, rms = {np.sqrt(np.mean((f-fd)**2)):.1e}"
              f"  (field SD {fd.std():.3f}; EV wall amplitude up to 0.5)")

    print("\n(a) static spins (no diffusion), T2 off, R = 5 um, CBV = 3 %:")
    Hs = ms.combine_histories([ms.phase_history(V[s][0.03], 5.0, 20000, seed=s, static=True) for s in (0, 1)])
    for so2 in (0.3, 0.6, 0.9):
        S = ms.signals_from_history(Hs, so2)
        print(f"    SO2 = {so2:.1f}: S(TE = 100 ms)  total = {S['total'][0, SE]:.12f}  EV = {S['ev'][0, SE]:.12f}"
              f"  IV = {S['iv'][0, SE]:.12f}   (S at 41.6 ms = {S['total'][0, 12]:.3f})")

    d = np.load(DICT)
    SO2g, Rg, Cg, sig = d["SO2"], d["R"], d["CBV"], d["signals"]
    ir5, ic3, is6 = list(Rg).index(5.0), list(Cg).index(0.03), int(np.argmin(abs(SO2g - 0.6)))
    print("\n(b) GE decay rate R2* (1/s; log-linear fit to the 13 part-A echoes of the total signal, T2 off), R = 5 um")
    print("    vs CBV at SO2 = 60 %: " + "  ".join(f"{c*100:.0f}%: {ms.fit_rate(TE, sig[is6, ir5, ic], A):.2f}"
                                                for ic, c in enumerate(Cg)))
    print("    vs SO2 at CBV = 3 %:  " + "  ".join(f"{s*100:.0f}%: {ms.fit_rate(TE, sig[i, ir5, ic3], A):.2f}"
                                                 for i, s in enumerate(SO2g) if i % 2 == 0 or s == 0.95))

    print("\n(c) spin-echo attenuation at TE = 100 ms vs R (SO2 = 60 %, CBV = 3 %, T2 off)")
    extra = {}
    for R in (1.0, 1.5):
        hs = [ms.phase_history(V[s][0.03], R, 20000, seed=s, dt_max=0.05) for s in (0, 1)]
        extra[R] = ms.signals_from_history(ms.combine_histories(hs), 0.6)
    rows = [(R, extra[R]["total"][0, SE], extra[R]["ev"][0, SE], "extra run, dt <= 0.05 ms") for R in extra]
    rows += [(R, sig[is6, ir, ic3, SE], d["ev"][is6, ir, ic3, SE], "dictionary") for ir, R in enumerate(Rg)]
    best = max(rows, key=lambda x: 1 - x[1])
    for R, st, se, src in rows:
        print(f"    R = {R:5.1f} um: 1 - S_SE  total = {1-st:.4f}  EV = {1-se:.4f}   dR2(SE) = {-np.log(st)/0.1:.2f} 1/s   [{src}]")
    print(f"    -> largest SE attenuation at R = {best[0]} um")

    print("\n(d) long-time EV R2' vs static-dephasing formula CBV*2*pi*42.58e6*(4/3)*pi*dchi_cgs*B0")
    late = A & (TE >= 20.0)
    for cbv in (0.01, 0.03, 0.06):
        theory = ms.static_dephasing_R2p(cbv, 0.6)
        out = []
        hs = ms.combine_histories([ms.phase_history(V[s][cbv], 5.0, 20000, seed=s, static=True) for s in (0, 1)])
        ev = ms.signals_from_history(hs, 0.6)["ev"][0]
        out.append(("static", ms.fit_rate(TE, ev, late), ms.gesfide_rates(ev)[1]))
        ic = list(Cg).index(cbv)
        for R in (30.0, 10.0, 5.0):
            ev = d["ev"][is6, list(Rg).index(R), ic]
            out.append((f"R={R:.0f}um", ms.fit_rate(TE, ev, late), ms.gesfide_rates(ev)[1]))
        print(f"    CBV {cbv*100:.0f}%, SO2 60%: theory {theory:.2f} 1/s | " +
              " | ".join(f"{n}: slope(20-42 ms) {a:.2f}, GESFIDE (a-b)/2 {b:.2f}" for n, a, b in out))
    print(f"\nvalidation runtime {time.time()-t0:.0f} s")


if __name__ == "__main__":
    main()

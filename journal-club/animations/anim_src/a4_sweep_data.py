"""Data for animation A4 (anim_sweep): GESFIDE curves while one parameter sweeps.

Segment 1  SO2 90 -> 40 %   (R = 5 um, CBV = 3 %)  : exact, continuous (phase-history rescaling)
Segment 2  R   2 -> 30 um   (SO2 = 60 %, CBV = 3 %): 29 log-spaced radii, same 2 geometries and
                                                      common random numbers -> smooth in R
Segment 3  CBV 1 -> 6 %     (SO2 = 60 %, R = 5 um) : nested geometries, CBV step 0.25 %
                                                      (vessels added one batch at a time)
All runs use the dictionary's seeds (0, 1) and 20000 spins per seed, so the nodes at
(R, CBV) = (5 um, 3 %) etc. are identical to the dictionary entries.
Output: anim_src/_a4_sweep_cache.npz (complex compartment means, T2-free).
"""
import os
import sys
import time
import multiprocessing as mp

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mrvf_sim as ms  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "_a4_sweep_cache.npz")
SEEDS = (0, 1)
N_SPINS = 20000
SO2_FIX, R_FIX, CBV_FIX = 0.60, 5.0, 0.03
SO2_FINE = np.round(np.linspace(0.90, 0.40, 201), 6)
R_SWEEP = np.geomspace(2.0, 30.0, 29)
CBV_SWEEP = np.round(np.arange(0.01, 0.06001, 0.0025), 6)
N_WORKERS = int(os.environ.get("MRVF_WORKERS", "3"))
VOX = {}


def _run(args):
    kind, seed, val = args
    if kind == "R":
        H = ms.phase_history(VOX[seed][CBV_FIX], val, n_spins=N_SPINS, seed=seed)
        S = ms.signals_from_history(H, SO2_FIX)
    elif kind == "CBV":
        H = ms.phase_history(VOX[seed][val], R_FIX, n_spins=N_SPINS, seed=seed)
        S = ms.signals_from_history(H, SO2_FIX)
    else:  # SO2 segment: fine SO2 grid from one history
        H = ms.phase_history(VOX[seed][CBV_FIX], R_FIX, n_spins=N_SPINS, seed=seed)
        S = ms.signals_from_history(H, SO2_FINE)
    n_iv = int(H["iv"].sum())
    return kind, seed, val, S["M_ev"], S["M_iv"], len(H["iv"]) - n_iv, n_iv


def _pool(results, kind, values):
    """Pool seeds: spin-count weighted complex means -> total magnitude (T2 off)."""
    out = []
    for v in values:
        rs = [r for r in results if r[0] == kind and np.isclose(r[2], v)]
        nev = sum(r[5] for r in rs)
        niv = sum(r[6] for r in rs)
        Mev = sum(r[5] * r[3] for r in rs) / nev
        Miv = sum(r[6] * r[4] for r in rs) / max(niv, 1)
        f = niv / (nev + niv)
        out.append(np.abs((1 - f) * Mev + f * Miv))
    return np.array(out)


def main():
    t0 = time.time()
    for s in SEEDS:
        VOX[s] = ms.nested_voxels(list(CBV_SWEEP), seed=s, verbose=True)
    tasks = ([("SO2", s, 0.0) for s in SEEDS]
             + [("CBV", s, float(c)) for s in SEEDS for c in CBV_SWEEP[::-1]]
             + [("R", s, float(r)) for s in SEEDS for r in R_SWEEP])
    print(f"{len(tasks)} runs on {N_WORKERS} workers")
    with mp.get_context("fork").Pool(N_WORKERS) as pool:
        res = list(pool.imap_unordered(_run, tasks))
    so2_tot = _pool(res, "SO2", [0.0])[0]             # (nSO2_fine, 38)
    r_tot = _pool(res, "R", R_SWEEP)[:, 0]             # (nR, 38)
    cbv_tot = _pool(res, "CBV", CBV_SWEEP)[:, 0]       # (nCBV, 38)
    np.savez(OUT, TE=ms.TE_MS, SO2_fine=SO2_FINE, S_so2=so2_tot, R_sweep=R_SWEEP, S_R=r_tot,
             CBV_sweep=CBV_SWEEP, S_cbv=cbv_tot, SO2_fix=SO2_FIX, R_fix=R_FIX, CBV_fix=CBV_FIX)
    print(f"saved {OUT} in {time.time()-t0:.0f} s")


if __name__ == "__main__":
    main()

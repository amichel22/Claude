"""Data for animation A2 (anim_matching): a fresh 'measured' voxel for dictionary matching.

Simulates SO2 = 60 %, CBV = 3 %, R = 5 um with geometry/spin seeds that the dictionary
did NOT use (dictionary: seeds 0, 1), same spin count and time grid, T2 off (as the
dictionary's `signals`).  Stores the complex total signal so noise can be added later.
Output: anim_src/_a2_measured_cache.npz
"""
import os
import sys
import time
import multiprocessing as mp

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mrvf_sim as ms  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "_a2_measured_cache.npz")
SEEDS = (101, 102)
N_SPINS = 20000
SO2_TRUE, CBV_TRUE, R_TRUE = 0.60, 0.03, 5.0


def _run(seed):
    g = ms.make_voxel(CBV_TRUE, seed)
    return ms.phase_history(g, R_TRUE, n_spins=N_SPINS, seed=seed)


def main():
    t0 = time.time()
    with mp.get_context("fork").Pool(len(SEEDS)) as pool:
        hs = pool.map(_run, SEEDS)
    H = ms.combine_histories(hs)
    S = ms.signals_from_history(H, SO2_TRUE)
    np.savez(OUT, TE=ms.TE_MS, M_total=S["M_total"][0], total=S["total"][0], f_iv=S["f_iv"],
             SO2_true=SO2_TRUE, CBV_true=CBV_TRUE, R_true=R_TRUE, seeds=np.array(SEEDS),
             n_spins_per_seed=N_SPINS, cbv_actual=H["cbv"])
    print(f"saved {OUT} in {time.time()-t0:.0f} s; f_iv = {S['f_iv']:.4f}, CBV(actual) = {H['cbv']:.4f}")
    print("S:", np.round(S["total"][0], 4))


if __name__ == "__main__":
    main()

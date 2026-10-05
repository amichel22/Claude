"""Build the small GESFIDE MRvF dictionary with mrvf_sim.py.

Grid: SO2 0.30..0.95 (step 0.05, free via the phase-history trick),
      R in {2,3,4,5,7,10,15,20,30} um, CBV in {1..6} %.
Each (R, CBV) node = pooled spins of 2 geometry seeds (N_SPINS per seed).
Output: SCRATCH/assets/anim/mrvf_dict.npz   (axes order SO2, R, CBV, TE)
"""
import os
import sys
import time
import multiprocessing as mp

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mrvf_sim as ms  # noqa: E402

SCRATCH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(SCRATCH, "assets", "anim", "mrvf_dict.npz")

SO2 = np.round(np.arange(0.30, 0.951, 0.05), 2)
R_UM = np.array([2, 3, 4, 5, 7, 10, 15, 20, 30], float)
CBV = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.06])
SEEDS = (0, 1)
N_SPINS = 20000
N_WORKERS = int(os.environ.get("MRVF_WORKERS", "3"))

VOX = {}


def _task(args):
    seed, ic, ir = args
    g = VOX[seed][CBV[ic]]
    H = ms.phase_history(g, R_UM[ir], n_spins=N_SPINS, seed=seed)
    S = ms.signals_from_history(H, SO2)
    n_iv = int(H["iv"].sum())
    return seed, ic, ir, S["M_ev"], S["M_iv"], len(H["iv"]) - n_iv, n_iv, H["runtime"], H["n_stuck"]


def main():
    t_start = time.time()
    print(f"building geometries ({len(SEEDS)} seeds, nested CBV {CBV.tolist()}) ...")
    for s in SEEDS:
        VOX[s] = ms.nested_voxels(list(CBV), seed=s, verbose=True)
    t_geom = time.time() - t_start
    tasks = [(s, ic, ir) for s in SEEDS for ic in range(len(CBV)) for ir in range(len(R_UM))]
    # heaviest first (high CBV -> larger near table)
    tasks.sort(key=lambda a: -a[1])
    nS, nR, nC, nT = len(SO2), len(R_UM), len(CBV), len(ms.TE_MS)
    Mev = np.zeros((nS, nR, nC, nT), complex)
    Miv = np.zeros_like(Mev)
    Nev = np.zeros((nR, nC))
    Niv = np.zeros((nR, nC))
    run_times = []
    print(f"{len(tasks)} MC runs x {N_SPINS} spins x {len(ms.make_schedule()['dt'])} steps on {N_WORKERS} workers ...")
    ctx = mp.get_context("fork")
    with ctx.Pool(N_WORKERS) as pool:
        for k, (s, ic, ir, mev, miv, nev, niv, rt, stuck) in enumerate(pool.imap_unordered(_task, tasks)):
            Mev[:, ir, ic] += nev * mev        # pool spins across seeds
            Miv[:, ir, ic] += niv * miv
            Nev[ir, ic] += nev
            Niv[ir, ic] += niv
            run_times.append(rt)
            if (k + 1) % 12 == 0:
                print(f"  {k+1}/{len(tasks)} runs, elapsed {time.time()-t_start:.0f} s")
    Mev /= Nev[None, :, :, None]
    Miv /= np.maximum(Niv, 1)[None, :, :, None]
    f_iv = Niv / (Nev + Niv)
    F = f_iv[None, :, :, None]
    TE = ms.TE_MS
    w_ev = np.exp(-TE / ms.T2_EV_MS)
    w_iv = np.exp(-TE / ms.T2_IV_MS)
    tot = np.abs((1 - F) * Mev + F * Miv)
    tot_t2 = np.abs((1 - F) * w_ev * Mev + F * w_iv * Miv)
    np.savez_compressed(
        OUT,
        TE=TE, SO2=SO2, R=R_UM, CBV=CBV,
        signals=tot, signals_norm=tot / tot[..., :1],
        signals_T2=tot_t2, signals_T2_norm=tot_t2 / tot_t2[..., :1],
        ev=np.abs(Mev), iv=np.abs(Miv), M_ev=Mev, M_iv=Miv, f_iv=f_iv,
        axes=np.array("SO2,R,CBV,TE"),
        B0_T=ms.B0_T, gamma_rad_s_T=ms.GAMMA, Hct=ms.HCT, dchi0_cgs_per_Hct=ms.DCHI0_CGS,
        D_um2_ms=ms.D_UM2_MS, dt_max_ms=ms.DT_MAX_MS, t180_ms=ms.T180_MS,
        T2_ev_ms=ms.T2_EV_MS, T2_iv_ms=ms.T2_IV_MS, n_spins_per_seed=N_SPINS,
        seeds=np.array(SEEDS), vessels_per_pct=ms.VESSELS_PER_PCT,
        box_side_in_R=ms.box_side(), rc_near_in_R=ms.RC_NEAR, grid_h_in_R=ms.GRID_H,
        notes=np.array(
            "2D MC, infinite cylinders perpendicular to plane, one random B0 direction per vessel "
            "(cos theta golden-ratio stratified on [0,1], azimuth uniform), periodic voxel, "
            "non-overlapping vessels, impermeable walls, nested CBV geometries. "
            "signals = |(1-f_iv) M_ev + f_iv M_iv| with T2 off (pure susceptibility); "
            "signals_T2 adds exp(-TE/80 ms) to EV and exp(-TE/50 ms) to IV. "
            "M_ev/M_iv are complex T2-free compartment means; rebuild any T2 as "
            "|(1-f) e^{-TE/T2ev} M_ev + f e^{-TE/T2iv} M_iv|."),
    )
    t_tot = time.time() - t_start
    print(f"saved {OUT}  ({os.path.getsize(OUT)/1024:.0f} KB)")
    print(f"geometry build {t_geom:.0f} s | mean run {np.mean(run_times):.1f} s | TOTAL RUNTIME {t_tot:.0f} s ({t_tot/60:.1f} min)")


if __name__ == "__main__":
    main()

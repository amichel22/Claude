"""Playground run: one MC simulation per (approach, radius, seed), stored as signals over an SO2 grid.

Approaches
  3d     : 3D-ANA-MC, the reference (continuous space, analytical field)
  random : 2D, gridded field, each vessel with its own random B0 direction (Berman's best 2D approach)
  3b0    : 2D, same geometry and grid, MrVox-style field averaged over three B0 directions

Usage: python run_play.py [OUT_DIR] [N_WORKERS]   (defaults: out, up to 4 cores)
On ideacluster the same job list maps onto a PBS array (one job per task).
"""
import os
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('NUMBA_NUM_THREADS', '1')

import numpy as np
import multiprocessing

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bss_common as bc

RADII = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0]   # um
SO2_GRID = np.round(np.arange(0.30, 0.951, 0.01), 2)                     # 66 values
SEEDS_3D = [1, 2, 3, 4]        # 1, 2 = "measured" test signals; 3, 4 = 3D dictionary
SEEDS_2D = [11, 12]
NSPINS = 10_000
N_GRID = 512                   # 2D grid: about 8 points across each vessel diameter


def jobs():
    out = [('3d', r, s) for s in SEEDS_3D for r in RADII]
    out += [(a, r, s) for a in ('random', '3b0') for s in SEEDS_2D for r in RADII]
    return out


def run(job, out_dir):
    approach, r_um, seed = job
    path = os.path.join(out_dir, f'{approach}_R{r_um:g}_s{seed}.npz')
    if os.path.exists(path):
        return path, 0.0
    t0 = time.time()
    a = r_um / bc.R0_UM
    if approach == '3d':
        vox = bc.scaled3d(bc.voxel3d(seed), a)
        phi, iv = bc.phase_history(vox, NSPINS, seed=100 + seed, ndims=3)
    else:
        base = bc.voxel2d(seed)
        grid = bc.grids2d(base, N_GRID)[approach]
        vox = bc.grid_voxel(grid, base.size.value * a)
        phi, iv = bc.phase_history(vox, NSPINS, seed=100 + seed, ndims=2)
    tot, ev, ivs = bc.signals_from_history(phi, iv, SO2_GRID)
    np.savez(path, tot=tot.astype(np.float32), ev=ev.astype(np.float32), iv=ivs.astype(np.float32),
             iv_fraction=iv.mean(), so2=SO2_GRID, te=bc.TE, radius_um=r_um, seed=seed, approach=approach)
    return path, time.time() - t0


def _run(args):
    return run(*args)


if __name__ == '__main__':
    out_dir = sys.argv[1] if len(sys.argv) > 1 else 'out'
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else min(4, os.cpu_count() or 1)
    os.makedirs(out_dir, exist_ok=True)
    todo = [(j, out_dir) for j in jobs()]
    # longest jobs first
    todo.sort(key=lambda x: (x[0][0] != '3d', -x[0][1]))
    print(f'Running {len(todo)} simulations on {workers} cores into the folder "{out_dir}".', flush=True)
    print('The 3D runs come first and take about 30-40 s each; the 2D runs at the end take a few seconds.', flush=True)
    print('Runs that are already in the folder are skipped, so you can stop (Ctrl+C) and start again.', flush=True)
    t0 = time.time()
    # 'spawn' is how Mac and Windows start worker processes; use it everywhere so every computer behaves the same
    with multiprocessing.get_context('spawn').Pool(workers) as pool:
        for i, (p, dt) in enumerate(pool.imap_unordered(_run, todo), 1):
            note = 'already done, skipped' if dt == 0 else f'{dt:.0f} s'
            print(f'[{i}/{len(todo)}] {os.path.basename(p)} {note} (elapsed {time.time() - t0:.0f} s)', flush=True)
    print(f'All {len(todo)} runs done in {(time.time() - t0) / 60:.1f} min. Next: python analyze_play.py', flush=True)

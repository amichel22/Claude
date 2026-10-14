# bold-slide-figures

Makes the five simulation figures for slides 8 and 11–14 of the BOLDsωimsuite talk, from scratch.

## Quick start

1. Put this folder in your `Research Work Fan Rotation` folder on the Desktop.
2. Open Jupyter the way you opened Jack's notebook, then open `make_slide_figures.ipynb`.
3. Choose the kernel **BOLDswimsuite (py3.10)** (Kernel → Change kernel).
4. Run the cells from top to bottom with Shift + Enter. Cell 2 takes about 10–15 minutes; everything else is quick.
5. Your figures are in `figs/`. Cell 9 should end with "All 16 checks match".

The full step-by-step guide, with troubleshooting, is the doc "Step by step: run the simulations and make your slide figures".

## Files

| File | What it is |
| --- | --- |
| `make_slide_figures.ipynb` | The notebook: runs the simulations, draws and saves each figure, compares with the reference run |
| `bss_common.py` | Shared settings: the lab's 38 GESFIDE echo times, 3 T, 3% blood volume, Hct 0.40, Δχ0 0.27 ppm (cgs), ADC 1 µm²/ms; voxel builders; the phase-history helper |
| `run_play.py` | One simulation per (recipe, radius, seed); the notebook calls it 96 times. Also runs on its own: `python run_play.py` |
| `analyze_play.py` | The figures as a plain script: `python analyze_play.py` (after `run_play.py`) |
| `compare_to_mine.py` | Compares `figs/summary.json` with `reference_summary.json` |
| `check_install.py` | Jack's install check |
| `reference_summary.json` | The reference numbers |

## What is simulated

Three recipes, each at 12 vessel radii from 1 to 50 µm (a 5 µm geometry scaled), 150 vessels, 3% blood volume, 10,000 spins, T2 off, on the 38-echo GESFIDE timing (180° pulse at 50 ms, spin echo at 100 ms):

| Recipe | Description | Voxels (seeds) |
| --- | --- | --- |
| `3d` | 3D-ANA-MC, Berman's reference: continuous 3D voxel, analytical field, Monte Carlo diffusion | 1, 2 (test signals) and 3, 4 (3D dictionary) |
| `random` | 2D, field on a 512 × 512 grid, each vessel with its own random B0 direction | 11, 12 |
| `3b0` | Same 2D geometry and grid, field = 1/3 (B0 along the vessels) + 2/3 (B0 perpendicular), as in MrVox | 11, 12 |

Each run is stored for 66 SO2 values (30–95%) by rescaling the accumulated phase (Jack's notebook, section 7). The matching uses the normalized inner product over the 38 echoes, jointly over radius and SO2.

## Caveats

A toy version of the project: one blood volume, no T2, total signal (blood + tissue), 2 voxels × 10,000 spins per dictionary entry (Berman used 10 × 40,000).

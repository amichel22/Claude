# Playground runs

Small BOLDsωimsuite experiments on the lab's 38-echo GESFIDE timing. They make the figures on slides 8 and 11–14 of the toolbox talk and the numbers in "Things to try" of the rotation guide.

## What they do

`run_play.py` runs one Monte Carlo simulation per (approach, radius, seed) and stores the 38-echo signal for 66 SO2 values (30–95%), using the phase-history rescaling from Jack's notebook (section 7). Three approaches:

| Name in the code | Recipe | Seeds |
| --- | --- | --- |
| `3d` | 3D-ANA-MC, the reference: continuous 3D voxel, analytical field | 1, 2 (test signals) and 3, 4 (3D dictionary) |
| `random` | 2D, field on a 512 × 512 grid, each vessel with its own random B0 direction | 11, 12 |
| `3b0` | Same 2D geometry and grid, field = 1/3 (B0 along the vessels) + 2/3 (B0 perpendicular), as in MrVox and Jack's notebook | 11, 12 |

Fixed settings: 3 T, 150 vessels, 3% blood volume, Hct 0.40, Δχ0 0.27 ppm (cgs) per unit Hct, ADC 1 µm²/ms, 10,000 spins, T2 decay off, 12 radii from 1 to 50 µm made by scaling a 5 µm geometry. `bss_common.py` holds the timing (echo 35 at 121.360 ms, as in `Echotimes_38.mat`), the physiology and the helpers, copied in spirit from Jack's notebook.

`analyze_play.py` averages seeds, matches the 3D test signals (SO2 50, 60, 70, 80%) against the three dictionaries jointly over radius and SO2 by normalized inner product, repeats it with Rician noise at SNR 100 on the first echo, and draws:

| Figure | Shows |
| --- | --- |
| `fp_radius.png` | Estimated vs true radius at SO2 60% for each dictionary |
| `fp_so2.png` | SO2 error vs radius for each dictionary |
| `boxerman_gesfide.png` | Gradient-echo (29 ms) and spin-echo (100 ms) decay rates vs radius |
| `knobs.png` | GESFIDE signal as SO2 and radius change (3D reference) |
| `fieldmaps.png` | Field maps of one 2D geometry, own direction per vessel vs 3B0 |

`summary.json` has the numbers behind them.

## Run it

Use the Python 3.10 `bold` environment from Jack's setup guide.

```bash
python check_install.py            # should print the field check, S(50 ms) = 0.9334, S(100 ms) = 0.9570, INSTALL OK
MPLBACKEND=Agg python run_play.py out 4      # 96 runs, about 8 minutes on 4 cores; skips runs already in out/
MPLBACKEND=Agg python analyze_play.py out figs
```

On ideacluster, each entry of `run_play.jobs()` is one array task.

## Caveats

This is a toy version of the project: one blood volume, no T2, no noise in the main figures, and 2 voxels × 10,000 spins per entry (Berman used 10 voxels × 40,000 spins). The 3D control dictionary is a few SO2 points low, which is probably voxel-to-voxel noise; more voxels per entry should tell.

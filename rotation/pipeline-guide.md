# What to do with Jack's BOLDswimsuite notebook

Markdown copy of the shared doc of the same name (the doc version has the drawn week-by-week timeline).

For Andres Michel · Fan Lab rotation, Fall 2026 · Mon Oct 5, 2026 (week 2 of 10)

Sources: `BOLDswimsuite_GESFIDE38_walkthrough.ipynb` and its saved outputs, Jack's setup guide (Sep 28), your rotation agreement, Chaussé, Berman & Chen 2025, and Lin, Lu & Fan 2025 (sections 3, 6.2).

## 1. What this notebook is

Jack wrote it to test whether BOLDswimsuite (Python) can replace MRVox (Matlab) for the lab's MR vascular fingerprinting (MRvF) dictionaries. It runs each module (vessels, geometry, spins, sequence, deterministic diffuser) on the echo times of the lab's 38-echo GESFIDE. Then it compares simulation methods with the 3D reference, including an MRVox-like "3B0" voxel, and ends with a dictionary strategy in which SO2, T2 and Hct are nearly free. Outputs are saved, so you can read it before running it.

How it maps to your agreement: sections 1-6 cover "learn the acquisition and the biophysical model". Reproducing CBV, R and SO2 on lab data needs the section 7 dictionary plus matching. The ML vs matching comparison starts once that dictionary exists.

## 2. Before you run anything

**Files, in one folder:** the toolbox source (the folder with `pyproject.toml`), the notebook, and `Echotimes_38.mat` (cell 1 loads it by name). If Jack didn't send the `.mat`, ask him. As a stopgap only, replace the `loadmat` line with the times the notebook states:

```python
TE = np.r_[4.988 + 3.048*np.arange(13),    # part A: 13 FID echoes, 4.988-41.564 ms
           57.328 + 3.048*np.arange(25)]   # B+C: 25 echoes, 57.328-130.48 ms, SE at 100 ms
TE[np.isclose(TE, 121.336)] = 121.360      # the file's echo 35 (see Gotchas)
```

This reproduces the saved schedule (681 steps, 180° before step 260).

**Environment.** The package installs only on Python 3.10 and needs numpy < 2. Make a separate `bold` env; `mrf` is Python 3.11 and will be refused. The notebook's saved kernel is `mrf` (metadata: Python 3.11.4), so after registering `bold` choose **Kernel → Change kernel → BOLDswimsuite (py3.10)**.

Conda (Linux, WSL2, cluster; on a Mac the same commands should work, but Jack only tested Linux and Windows):

```bash
conda create -n bold python=3.10 -y && conda activate bold
cd ~/BOLDswimsuite-main                 # path without spaces
python -m pip install ".[lessons]"      # keep the quotes
python -m ipykernel install --user --name bold --display-name "BOLDswimsuite (py3.10)"
```

uv (the route Jack verified):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
uv venv -p 3.10 ~/envs/bold && source ~/envs/bold/bin/activate
uv pip install ~/BOLDswimsuite-main ipykernel jupyter
python -m ipykernel install --user --name bold --display-name "BOLDswimsuite (py3.10)"
```

Windows: Anaconda Prompt, same conda commands, a folder like `C:\code\`. For offline cluster installs, PyVista and pinned versions, see the setup guide.

## 3. Day 1 checklist

1. Install; `python --version` must say 3.10.x.
2. Run Check 1 from the setup guide (`check_install.py`). The seeded numbers should match; timings won't:
   ```
   field check: 2.035752e-07 vs 2.035752e-07 -> OK
   mini GESFIDE: 32 vessels, S(50 ms) = 0.9334, S(100 ms) = 0.9570
   INSTALL OK
   ```
3. `MPLBACKEND=Agg python examples/3D-ANA-MC_script.py` should print a `ContinuousVoxel3D` with 49 cylinders and walk 600 steps in about 15 s. "Display functionality disabled" only means no PyVista.
4. Run the notebook headless (15-20 min on 2 cores for Jack; yours will differ). `kernel_name=bold` overrides the saved `mrf` kernel:
   ```bash
   jupyter nbconvert --to notebook --execute \
     --ExecutePreprocessor.kernel_name=bold --ExecutePreprocessor.timeout=3600 \
     BOLDswimsuite_GESFIDE38_walkthrough.ipynb --output executed.ipynb
   ```
5. Check `executed.ipynb` against section 4, then work through it interactively, one section per sitting.
6. `python -m pip freeze > bold_env_lock.txt`, and start a Git repo for your own scripts now.

## 4. Section-by-section reading guide

Numbers come from Jack's saved run. Seeded results should match closely on the known-good versions; run times depend on your machine.

**0 · Setup (cells 0-1).** ΔTE = 3.048 ms isn't a multiple of a practical step, so each gap between events is split into equal steps ≤ 0.2 ms and every echo and the 180° land exactly on a step. Also converts SO2 to Δχ and blood T2.
- Check: `38 echoes | A: 13 | B: 15 | C: 10`; 681 steps, dt 0.1890-0.1995 ms, 180° before step 260; irregularities `{121.36: 3.072, 124.384: 3.024}`; Δχ 4.320e-8 cgs (5.429e-7 SI); blood T2 44.6 ms.
- Try: `SO2 = 0.4` gives Δχ 6.48e-8 and blood T2 ≈ 23.7 ms. `DT = 0.1` gives 1325 steps and about twice the run time.

**1 · BOLDvessel.** Vessels know whether a point is inside them and what field they make (Chaussé Eq. 1-2). 2D cylinders stay perpendicular to the plane, and each gets its own random B0 direction.
- Check: EV at r = 2R (2.0358e-7 T) and IV at θ = 90° (−2.7143e-7 T) match Eq. 1 exactly. On the orientation plot, IV goes from +543 nT at θ = 0 through 0 at 54.7° to −271 nT at 90°.
- Try: θ = 0. EV should vanish and IV should equal B0·Δχ_SI/3.

**2a · Continuous voxels.** `from_random` adds random cylinders until the target CBV is reached. Fields are exact analytical sums, so each spin step costs O(N_vessels).
- Check: 3D, 150 cylinders in a 0.408 mm cube, `get_CBV()` 0.0302 vs sampled 0.0307. 2D, 151 cylinders, 0.0302 vs 0.0291. The 2D and 3D EV field histograms overlap.
- Try: `CBV=0.06` with `num_vessels=150` kept. The voxel shrinks by about 1/√2 (to ~0.29 mm) instead of gaining vessels.

**2b · Discrete voxels, FFT fields.** Spins look the field up on a grid (memory ∝ N³). FFT dipole convolution is the only route for custom masks.
- Check: 96³ gives 2.4 points per vessel diameter (Berman recommends ≥ 6). FFT EV RMS error is 56.3 nT cyclic and 53.5 nT padded, against an EV field SD of 75.0 nT. Jack's integer-lag kernel, which removes a half-cell offset in the package's FFT kernel, gives 37.8 nT.
- Try: N = 192 (~4.7 points per diameter, ~113 MB). The error should fall, showing how much of it is under-sampling.

**3 · BOLDspins.** Monte Carlo walk: Gaussian steps with σ = √(2·D·dt), periodic boundaries, steps across walls re-proposed, and a phase of 2π·γ·ΔBz·dt per step.
- Check: MSD after 20 ms is 39.91 µm² (theory 40); σ = 0.632 µm; IV spin fraction 0.0296 (≈ CBV).
- Try: ADC `2e-3` should give an MSD near 80 µm². `IV=False` leaves no IV spins, so the IV signal becomes NaN.

**4 · BOLDsequence: the reference.** 3D-CTN-CYL-ANA-MC with pulses, phase and per-compartment T2; total, EV and IV are averages of Mxy. In the middle panel (T2 off), the no-diffusion signal refocuses to 1.0 at 100 ms. Diffusion attenuates the spin echo, and that gap carries the radius information.
- Check: 31 s for 150 vessels × 10k spins × 681 steps. Fitted R2 = 14.8 s⁻¹ (input 12.5; most of the excess likely comes from IV and diffusion). R2' = 1.81 s⁻¹, against a static-dephasing prediction of 4.36.
- Try: set `SO2 = 0.4` in cell 1 (Δχ and blood T2 follow) and rerun through cell 16. Part A should decay faster and R2' should rise.

**5 · BOLDdeterministic.** 2D diffusion as a convolution of the magnetization grid at each step, which is MRVox's engine. No MC noise, but the grid has to resolve σ ≈ 0.63 µm.
- Check: 1050² grid at 0.597 µm, 15 × 15 kernel, 143 s, against 20 s for 2D-ANA-MC.
- Try: halve `Ndd`. It runs faster, and the curves drift from MC because the kernel is too coarse (Chaussé Fig. 9).

**6 · Head-to-head.** NRMSD over the 38 echoes against the reference, from one voxel each, so a few % is noise.

| Method | Total % | EV % | IV % |
|---|---|---|---|
| 3D-FFT-MC | 2.5 | 1.8 | 24.5 |
| 2D-ANA-MC | 0.5 | 0.6 | 10.6 |
| 2D-ANA-DD | 1.1 | 1.5 | 11.3 |
| 2D-ANA-DD-3B0 (≈ MRVox) | 3.7 | 0.7 | 93.9 |

- Check: 3B0 IV offsets have a median of −0.3 nT, against −74 nT with random B0. The 3D-FFT-MC walk took 1 s, against 31 s.
- Try: rebuild both voxels with another seed. The small EV numbers will move, but the 3B0 IV error should stay large because it's systematic.

**7 · Dictionary strategy.** Records each spin's sign-flipped phase at the 38 echoes once, then rescales it for any SO2, Hct, B0 or T2. R reuses a scaled geometry with a new walk. Each CBV needs its own geometry.
- Check: history and direct simulation agree to 1.7e-14 (SO2 0.6) and 1.8e-14 (0.75). That's round-off, so the trick is exact. 66 SO2 entries take 0.66 s. Radius runs take 23-26 s; fitted R2' rises from about 0.2 to 3.3 s⁻¹ between R = 2 and 20 µm. CBV log-scaling from 3% is off by 0.023 at 1.5% (total drop only 0.095) and 0.012 at 6%, so simulate each CBV.
- Try: an Hct axis from the same `Phi`. Run a second 3% seed too: if the two seeds differ by close to 0.02, much of the CBV error is geometry noise.

**8 · Practical notes.** Read it; the key items are in section 7 below.

## 5. Concepts to have straight

- **GESFIDE.** 90° at 0, 180° at 50 ms, SE at 100 ms. Part A (13 echoes) decays at R2* = R2 + R2'. Part B (15 echoes, up to the SE) refocuses, with a fitted slope of R2 − R2'. Part C (10 echoes) dephases again. So R2 = (slope_A + slope_B)/2 and R2' = (slope_A − slope_B)/2.
- **R2, R2*, R2'.** R2 = 1/T2 is irreversible. R2' comes from static field offsets and can be refocused. R2* = R2 + R2'. Static limit: R2' = CBV·γ·(4/3)π·Δχ_cgs·B0, where CBV and (1−SO2) appear as a product, so expect them to trade off.
- **GE vs SE.** Gradient echoes see all dephasing. The SE refocuses static dephasing around large vessels, but not the dephasing of water that diffuses through the gradients around small vessels (a few µm). The GE/SE contrast encodes R.
- **EV, IV, total.** Means over spins outside vessels, inside them, and all spins (so IV is weighted by ≈ CBV). IV uses one Δχ and has no RBC or flow effects.
- **Δχ units.** The package uses cgs and multiplies by 4π internally. Δχ = 0.27e-6·Hct·(1−SO2) cgs = 3.39e-6·Hct·(1−SO2) SI.
- **Method names.** [2D|3D]-[CTN|GRD]-[CYL|SPH|VAN]-[ANA|FFT]-[MC|DD]. The reference is 3D-CTN-CYL-ANA-MC (shortened to 3D-ANA-MC in the notebook). MRVox is 2D-GRD-CYL-FFT-DD with the 3B0 shortcut.
- **Why 2D random-B0 works.** An infinite cylinder's field depends only on the in-plane position relative to its axis and on θ. A random B0 direction per 2D cylinder reproduces the 3D spread of orientations.
- **Why 3B0 fails for IV.** It uses ⅓ of the parallel field plus ⅔ of the perpendicular field: Δχ/6·(⅓·2 + ⅔·(−1)) = 0. All IV spins sit near 0 nT, so blood barely dephases. With IV in the signal, that biases radius; the notebook cites Delphin et al.'s 14-19% overestimation in vivo.
- **Why SO2 and T2 are free, but R and CBV aren't.** With impermeable walls and ideal pulses, the walk doesn't depend on Δχ and the field is ∝ B0·Δχ, so phase scales linearly with Δχ (SO2, Hct, B0). Spins never change compartment, so T2 is a per-spin factor exp(−TE/T2). R changes vessel size relative to the diffusion distance (~16 µm per axis over 130 ms) and CBV changes the geometry, so both need new walks.

## 6. From notebook to project

**(a) Reproduce the reference.** Move cell 1, the cell 5 voxel, `gesfide_mc` (cell 16) and `phase_history`/`signal_from_history` (cell 25) into a module in your repo. Run it as a script and match CBV 0.0302, R2 14.8, R2' 1.81 and the 1e-14 check.

**(b) Small dictionary on the cluster.** `simulate_node.py TASK` maps the task index to one (R, CBV, seed), runs `phase_history` at a reference Δχ and saves `Phi` (38 × N_spins) and the IV flags to `.npz`. That's about 3 MB per run at 10k spins. Use `run_bss.pbs` from the setup guide (one core per array task) and test one task with `qsub -I` first. A second script expands the histories over SO2, plus Hct and T2 if needed.

Cost: the notebook gives 30-40 s per node at 3% CBV, so 25 R × 15 CBV × 3 seeds ≈ 1,100 runs ≈ 10-12 core-hours. Run time grows with vessel count (12 s at 1.5% CBV and 50 s at 6%, with 8k spins), so a fixed-size voxel going up to 25% CBV costs several times more. Debug with about 8 R × 6 CBV × 2 seeds.

Decide with Jack before the big run:

- EV-only or EV+IV
- Hct fixed at 0.4, or an axis
- blood T2(SO2) model (Lu 2012 is a placeholder) and tissue T2 (fixed 80 ms, or an axis)
- ranges and spacing (Christen 2014: SO2 0-100%, CBV 0.5-25%, R 1-25 µm)
- spins and seeds (notebook: 10k × 1; Berman: 40k spins × 10 voxels)
- 2D or 3D (2D-ANA-MC was within 0.6% on EV), ANA or FFT (ANA for cylinders), time step ≤ 0.2 ms
- voxel sizing: fixed vessel count (flat cost, smaller voxels) or fixed size (cost grows with CBV); keep it well above ~16 µm

**(c) Dictionary matching.**

```python
Dn = D / np.linalg.norm(D, axis=1, keepdims=True)    # (n_entries, 38)
Sn = S / np.linalg.norm(S, axis=1, keepdims=True)    # (n_voxels, 38)
est = params[np.argmax(Sn @ Dn.T, axis=1)]           # SO2, CBV, R; chunk voxels if needed
```

Test on signals that aren't in the dictionary: a held-out seed, and R/CBV values between nodes. Add Rician noise and sweep SNR from 20 to 100:

```python
sigma = clean[:, :1] / snr      # SNR on the first echo; write your definition down
noisy = np.abs(clean + sigma*(rng.standard_normal(clean.shape) + 1j*rng.standard_normal(clean.shape)))
```

The review cites Wheeler et al.: matching needs SNR ≥ ~20. Report bias and RMSE per parameter against SNR.

**(d) ML estimator.** The review suggests an MLP as the simple starting point, with LSTMs (MARVEL) for when temporal order matters. Start with a small MLP: 38 inputs, outputs SO2, CBV and log R scaled to [0, 1]. Train on dictionary signals with fresh Rician noise every epoch (SO2 can be sampled continuously from the histories) and split train/test by geometry seed. Keep PyTorch or scikit-learn in a separate env and pass `.npz` files, so nothing upgrades numpy in `bold`. Compare with matching on:

- accuracy: bias and RMSE per parameter, on and off the grid
- speed: time per 100k voxels
- noise: the same SNR sweep
- out-of-range inputs (R or CBV beyond the training range): matching gives a low best score you can flag; a network gives a plausible number with no warning (review 6.2)

The review cites a network trained on simulated GESFIDE dictionaries that beat matching (Lin et al., ISMRM 2024). Ask Jack what code exists.

**(e) Lab data.** Once (c) and (d) work in simulation, run both on the dataset Jack picks and compare with the lab's current maps, starting with GM/WM ROI values.

## 7. Gotchas

| Trap | Fix |
|---|---|
| `requires a different Python: 3.11.x` | You're in `mrf`/base. `conda activate bold`, reinstall. |
| `No module named 'boldswimsuite'` in Jupyter | Kernel is still `mrf`. Change kernel; `jupyter kernelspec list`. |
| numba errors after another install | numpy went to 2.x. `python -m pip install "numpy<2"`; keep ML tools elsewhere. |
| Wrong numbers, no error | Units: mm (5 µm radius = `diameter=0.010`), ms, ADC mm²/s (`1e-3`), Δχ cgs, fields T. |
| Echoes one step off | `walk()` index j is t = (j+1)·dt; the shipped examples plot it shifted. |
| Pulse index | A pulse at index k acts at t = k·dt, before that step. Section 8's "index 250" is for a uniform 0.2 ms grid; here the 180° is before step 260. |
| Echo 35 = 121.360 ms | 0.024 ms off the regular train (121.336). Use what the scanner used; ask Jack. |
| Blood T2 | Not modeled; you pass `T2IV`. No T1 unless set. |
| `examples/2D-ANA-MC-G-MW_script.py` | Never calls `spins.step()`; vessel-index filter off by one. Don't copy it. |
| Slow cluster jobs | `*_NUM_THREADS=1`, `MPLBACKEND=Agg`, `NUMBA_CACHE_DIR` on scratch. |
| Coarse grids | ≥ 6 points per vessel diameter; DD warns below 5 kernel points and needs one ADC. |

## 8. Questions to bring to Jack

1. Which in-vivo dataset should I use, where does it live, and how do I get access?
2. Is it pre/post ferumoxytol (Christen-style ratio) or contrast-free? That decides whether tissue T2 and B0/B1 must be in the dictionary.
3. Echo 35: 121.360 or 121.336 ms on the scanner? Is the 40-echo notebook in the setup guide a different protocol?
4. The notebook says your manuscript uses an EV-only dictionary, but the MRVox runs used `ExtraVascOnly = 0`. Which should mine match?
5. Which MRVox model file and constants (Δχ0, Hct, T2s, ADC, voxel size, spins, step) made the current dictionary? Where are it and the matching code?
6. Parameter ranges and spacing (linear or log for R and CBV)?
7. Keep Lu 2012 for blood T2(SO2)?
8. What SNR do the lab data have, and how is it defined?
9. Do I have an ideacluster account? Torque or PBS Pro array syntax, which queue, where do outputs go?
10. Did you report the FFT half-cell kernel offset to the authors? Should I avoid the FFT path?
11. Is there code or a trained model from the GESFIDE neural-network work in the review (Lin et al., ISMRM 2024) to compare against?
12. Where should my Git repo live, and when is the end-of-quarter lab meeting?

## 9. Plan for weeks 2-10 (about 15-20 h/week)

| Week | Dates | Goal | Done when |
|---|---|---|---|
| 2 | Oct 5-9 | Install, checks, headless run, sections 0-4 | Numbers match; questions to Jack |
| 3 | Oct 12-16 | Sections 5-7 and the "try" exercises; module in your repo | Reference reproduced from a script |
| 4 | Oct 19-23 | Design decisions with Jack; `simulate_node.py`; test task, small array | Phase histories on disk |
| 5 | Oct 26-30 | Assemble dictionary; matching; SNR sweep | Bias/RMSE vs SNR for SO2, CBV, R |
| 6 | Nov 2-6 | MLP with Rician noise, separate env | First ML vs matching comparison |
| 7 | Nov 9-13 (Nov 11 holiday) | Off-grid, speed, out-of-range tests | Comparison figures |
| 8 | Nov 16-20 | Lab data with both methods vs current maps | ROI table, example maps |
| 9 | Nov 23-25 (Thanksgiving) | Draft 2-page summary; repo README | Draft to Jack and Dr. Fan |
| 10 | Nov 30-Dec 4 | Lab meeting talk; final summary; repo handoff; talk with Dr. Fan about joining | Deliverables done |

If the lab data aren't ready by week 8, the simulation comparison (weeks 5-7) is a complete result on its own. Bring one new figure to each weekly meeting with Dr. Fan.

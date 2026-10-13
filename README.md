# Fan Lab rotation, Fall 2026

Andres Michel. Materials for the journal club on Lin, Lu & Fan 2025, for getting started with the BOLDsωimsuite notebook, and for the toolbox talk and rotation plan.

## Journal club: Lin, Lu & Fan 2025

Lin C-T, Lu H, Fan AP. Technical review of magnetic resonance fingerprinting applications in cerebral physiology. *Magn Reson Med*. 2025. doi:10.1002/mrm.70216 (open access, CC BY).

- **Deck (18 slides, ~16 min):** https://claude.ai/artifact/C8axe25QmaYyTtqDU7op5t. Present from the deck page, or download it as PowerPoint or PDF from there. The speaker notes on each slide hold the script. The link is private until you share it from the deck's Share menu.
- **Script doc:** https://claude.ai/code/artifact/fe8359be-567b-414b-852c-da32422a7799 (same text as `journal-club/script.md`, plus an "If asked" section with short answers to likely questions).

| Path | What it is |
| --- | --- |
| `journal-club/script.md` | Read-aloud script, one section per slide, `[click]` marks the two builds (slides 2 and 16), then "If asked" |
| `journal-club/slides/` | Slide sources: `deck.json` (order, sections) and one HTML file per slide, with the script in each slide's `<aside>` speaker notes. `assets.json` maps the `/_blob/...` image ids to the files below |
| `journal-club/figures/` | Figures 1–8 extracted from the paper PDF, plus the crops used on the slides (Fig. 3A/3B, two rows of Fig. 5, the three cover maps) |
| `journal-club/animations/` | The four animations (`assets/anim/`: GIF for the slides, MP4 copies) and the Python code that generated them (`anim_src/`) |
| `journal-club/slides-static.pdf` | Static PDF backup of the deck, animations at their last frame |

### The animations

All four are simulations written for this talk, for illustration. They are not results from the paper.

| File | Slide | What it shows | How it was made |
| --- | --- | --- | --- |
| `anim_fingerprint` | 3 | Flip angle and TR train, and white matter, gray matter and CSF signals separating into fingerprints | Extended phase graph simulation of an inversion-recovery FISP train, 500 TRs, flip angles 5–70°, TR 11–15 ms, 3 T tissue T1/T2 |
| `anim_matching` | 4 | A simulated 38-echo voxel signal matched against a 756-entry dictionary by normalized inner product; lands on the true SO₂ 60%, CBV 3%, R 5 µm | Dictionary from `mrvf_sim.py` below, test signal from a different geometry seed |
| `anim_spins` | 9 | Water spins around a 3 µm and a 30 µm vessel; the spin echo recovers 0.92 for the large vessel and 0.75 for the small one | One cylinder perpendicular to B0, 3 T, SO₂ 60%, Hct 0.4, D = 1 µm²/ms, 40,000 spins, T2 off |
| `anim_sweep` | 10 | The 38-echo GESFIDE signal as SO₂, then R, then CBV change | 2D Monte Carlo voxel of cylinders, each with its own random B0 direction (the BOLDsωimsuite 2D approach), lab's 38-echo timing, phase-history rescaling for SO₂, T2 off |

The GIFs, MP4 copies and the small dictionary (`mrvf_dict.npz`) are in `journal-club/animations/assets/anim/`. To regenerate them you need Python 3.10+ with numpy, scipy, matplotlib and pillow, and ffmpeg on the path. From `journal-club/animations/anim_src/`, run in this order:

1. `mrvf_make_dict.py` (about 4 minutes on 3 cores) and, optionally, `mrvf_validate.py` (static spins refocus fully at the spin echo; decay grows with CBV and with 1 − SO₂)
2. `a1_render_fingerprint.py`
3. `a2_measured_data.py`, then `a2_render_matching.py`
4. `a3_anim_spins.py` (simulates on the first run, then uses a local cache)
5. `a4_sweep_data.py`, then `a4_render_sweep.py`

`journal-club/slides-static.pdf` is a static backup of the 18 slides, with each animation shown at its last frame. It comes from a local approximation of the deck renderer, so use the deck's own download for the real thing.

## Rotation: what to do with Jack's notebook

`rotation/pipeline-guide.md` is a Markdown copy of the shared doc "What to do with Jack's BOLDsωimsuite notebook" (https://claude.ai/code/artifact/e0be4871-6b10-4fea-9427-abe90dbe6dce). It covers setup (Python 3.10 `bold` environment; the notebook's saved `mrf` kernel will not work), a day-1 checklist, a section-by-section reading guide with the numbers to reproduce, the concepts the notebook assumes, the path from notebook to dictionary, matching and ML comparison, gotchas, questions for Jack, and a plan for weeks 2–10.

## Toolbox talk: BOLDsωimsuite and the rotation plan

Slides for Dr. Fan on what BOLDsωimsuite does, the two papers behind it, what I ran, and how it becomes the rotation project. Papers: Chaussé, Berman & Chen, *Imaging Neuroscience* 2025 (doi:10.1162/imag_a_00519) and Berman, Chaussé, Hartung et al., *Imaging Neuroscience* 2026 (doi:10.1162/IMAG.a.1247).

- **Deck (18 slides, about 18 min):** https://claude.ai/artifact/G5HfX6UUMbQXMXH3u541xn. The speaker notes hold the script. Private until shared from the deck's Share menu.
- **Rotation guide (doc):** https://claude.ai/code/artifact/f628293b-a5d4-4ce6-b50d-1800cd3c4d02. What the toolbox is, the two papers, how to check your notebook run, the cluster, the project in three parts with a timeline to Dec 4, things to try with expected results, a glossary, and questions for Jack and Dr. Fan.

| Path | What it is |
| --- | --- |
| `toolbox-talk/script.md` | Read-aloud script, one section per slide |
| `toolbox-talk/slides/` | Slide sources (`deck.json`, one HTML file per slide); `assets.json` maps the `/_blob/...` ids to the figures below |
| `toolbox-talk/figures/papers/` | Figures cropped from the two papers |
| `toolbox-talk/figures/my-runs/` | Figures from my BOLDsωimsuite runs, plus `summary.json` with the numbers behind them |
| `toolbox-talk/playground/` | The scripts that made them (`check_install.py`, `run_play.py`, `analyze_play.py`); see its README |

The "my runs" figures are a toy version of the project on the lab's 38-echo GESFIDE timing: 3D-ANA-MC reference signals matched against 3D, 2D (own field direction per vessel) and 2D 3B0 (MrVox-style) dictionaries, at one blood volume, T2 off, 2 voxels × 10,000 spins per entry. The 3B0 dictionary got 3 of 12 radii right and shifted SO2 by up to 20 points; the other two missed only by one grid step.

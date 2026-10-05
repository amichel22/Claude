# Fan Lab rotation, Fall 2026

Andres Michel. Materials for the journal club on Lin, Lu & Fan 2025 and for getting started with the BOLDsωimsuite notebook.

## Journal club: Lin, Lu & Fan 2025

Lin C-T, Lu H, Fan AP. Technical review of magnetic resonance fingerprinting applications in cerebral physiology. *Magn Reson Med*. 2025. doi:10.1002/mrm.70216 (open access, CC BY).

- **Deck (18 slides, ~16 min):** https://claude.ai/artifact/C8axe25QmaYyTtqDU7op5t. Present from the deck page, or download it as PowerPoint or PDF from there. The speaker notes on each slide hold the script. The link is private until you share it from the deck's Share menu.
- **Script doc:** https://claude.ai/code/artifact/fe8359be-567b-414b-852c-da32422a7799 (same text as `journal-club/script.md`, plus an "If asked" section with short answers to likely questions).

| Path | What it is |
| --- | --- |
| `journal-club/script.md` | Read-aloud script, one section per slide, `[click]` marks the two builds (slides 2 and 16), then "If asked" |
| `journal-club/slides/` | Slide sources: `deck.json` (order, sections) and one HTML file per slide, with the script in each slide's `<aside>` speaker notes. `assets.json` maps the `/_blob/...` image ids to the files below |
| `journal-club/figures/` | Figures 1–8 extracted from the paper PDF, plus the crops used on the slides (Fig. 3A/3B, two rows of Fig. 5, the three cover maps) |
| `journal-club/animations/` | The four animations (GIF for the slides, MP4 copies) and the Python code that generated them |

### The animations

All four are simulations written for this talk, for illustration. They are not results from the paper.

| File | Slide | What it shows | How it was made |
| --- | --- | --- | --- |
| `anim_fingerprint` | 3 | Flip angle and TR train, and white matter, gray matter and CSF signals separating into fingerprints | Extended phase graph simulation of an inversion-recovery FISP train, 500 TRs, flip angles 5–70°, TR 11–15 ms, 3 T tissue T1/T2 |
| `anim_matching` | 4 | A simulated 38-echo voxel signal matched against a 756-entry dictionary by normalized inner product; lands on the true SO₂ 60%, CBV 3%, R 5 µm | Dictionary from `mrvf_sim.py` below, test signal from a different geometry seed |
| `anim_spins` | 9 | Water spins around a 3 µm and a 30 µm vessel; the spin echo recovers 0.92 for the large vessel and 0.75 for the small one | One cylinder perpendicular to B0, 3 T, SO₂ 60%, Hct 0.4, D = 1 µm²/ms, 40,000 spins, T2 off |
| `anim_sweep` | 10 | The 38-echo GESFIDE signal as SO₂, then R, then CBV change | 2D Monte Carlo voxel of cylinders, each with its own random B0 direction (the BOLDsωimsuite 2D approach), lab's 38-echo timing, phase-history rescaling for SO₂, T2 off |

To regenerate: Python 3.10+ with numpy, scipy, matplotlib, pillow and ffmpeg on the path. From `journal-club/animations/anim_src/`, run `mrvf_make_dict.py` (about 4 minutes on 3 cores), then the `a1_*`, `a2_*`, `a3_*` and `a4_*` render scripts. Outputs go to `journal-club/animations/assets/anim/`. `mrvf_validate.py` checks the simulator (static spins refocus fully at the spin echo, decay grows with CBV and with 1 − SO₂).

## Rotation: what to do with Jack's notebook

`rotation/pipeline-guide.md` is a Markdown copy of the shared doc "What to do with Jack's BOLDsωimsuite notebook" (https://claude.ai/code/artifact/e0be4871-6b10-4fea-9427-abe90dbe6dce). It covers setup (Python 3.10 `bold` environment; the notebook's saved `mrf` kernel will not work), a day-1 checklist, a section-by-section reading guide with the numbers to reproduce, the concepts the notebook assumes, the path from notebook to dictionary, matching and ML comparison, gotchas, questions for Jack, and a plan for weeks 2–10.

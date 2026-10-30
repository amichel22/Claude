"""Blood + tissue vs tissue only: does the 3B0 radius bias come from the blood signal?

Uses the 96 runs in out/ (no new simulations). Same matching as the main analysis,
done twice: once on the total signal (blood + tissue), once on the tissue (EV) signal only,
for the test signals and the dictionaries alike.

Usage: python tissue_only.py [RUN_DIR] [FIG_DIR]   (defaults: out, figs)
"""
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RUN_DIR = sys.argv[1] if len(sys.argv) > 1 else 'out'
FIG_DIR = sys.argv[2] if len(sys.argv) > 2 else 'figs'
os.makedirs(FIG_DIR, exist_ok=True)

INK, MUTED, HAIR = '#1B1B1F', '#5B5F66', '#D9DADC'
RED, BLUE = '#A8322D', '#1F4E79'
plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['Liberation Serif', 'Times New Roman', 'Times', 'DejaVu Serif'],
    'mathtext.fontset': 'stix', 'font.size': 20, 'axes.edgecolor': MUTED, 'axes.linewidth': 1.0,
    'axes.labelcolor': INK, 'xtick.color': MUTED, 'ytick.color': MUTED, 'xtick.labelcolor': INK,
    'ytick.labelcolor': INK, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'axes.labelsize': 21,
    'axes.titlesize': 22, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
})

runs = {}
for p in glob.glob(os.path.join(RUN_DIR, '*.npz')):
    z = np.load(p)
    runs[(str(z['approach']), float(z['radius_um']), int(z['seed']))] = z
if len(runs) < 96:
    sys.exit(f'Found {len(runs)} of 96 runs in "{RUN_DIR}". Run the simulations first.')
RADII = sorted({k[1] for k in runs})
radii = np.array(RADII)
SO2 = runs[('3d', 5.0, 1)]['so2']
i6 = int(np.argmin(np.abs(SO2 - 0.6)))


def mean_signal(approach, seeds, part):
    return np.stack([np.mean([runs[(approach, r, s)][part] for s in seeds], 0) for r in RADII])


def match(D, S):
    Dn = D.reshape(-1, D.shape[-1]); Dn = Dn / np.linalg.norm(Dn, axis=1, keepdims=True)
    Sn = S.reshape(-1, S.shape[-1]); Sn = Sn / np.linalg.norm(Sn, axis=1, keepdims=True)
    best = np.argmax(Sn @ Dn.T, axis=1)
    iR, iS = np.unravel_index(best, D.shape[:2])
    return iR.reshape(S.shape[:-1]), iS.reshape(S.shape[:-1])


RECIPES = {'2D, own field direction': ('random', [11, 12], BLUE, 's'),
           '2D, 3B0 average (MrVox-style)': ('3b0', [11, 12], RED, 'o')}
result = {}
for part in ('tot', 'ev'):
    test = mean_signal('3d', [1, 2], part)[:, i6, :]
    for name, (a, seeds, _, _) in RECIPES.items():
        iR, iS = match(mean_signal(a, seeds, part), test)
        result[(part, name)] = dict(estR=radii[iR], so2err=(SO2[iS] - 0.6) * 100)

fig, ax = plt.subplots(1, 2, figsize=(16, 6.2), sharex=True, sharey=True)
titles = {'tot': 'Blood + tissue signal', 'ev': 'Tissue signal only (no blood)'}
for a, part in zip(ax, ('tot', 'ev')):
    a.plot([0.8, 65], [0.8, 65], color=HAIR, lw=2, zorder=1)
    for name, (_, _, color, marker) in RECIPES.items():
        r = result[(part, name)]
        n_right = int((r['estR'] == radii).sum())
        a.plot(radii, r['estR'], marker, ms=11, mfc='white', mec=color, mew=2.4, ls='none', zorder=3,
               label=f'{name}: {n_right} of 12 right')
    a.set_xscale('log'); a.set_yscale('log')
    a.set_xticks([1, 3, 10, 30]); a.set_xticklabels(['1', '3', '10', '30'])
    a.set_yticks([1, 3, 10, 30]); a.set_yticklabels(['1', '3', '10', '30'])
    a.set_xlim(0.8, 65); a.set_ylim(0.8, 65)
    a.set_title(titles[part], pad=14)
    a.set_xlabel('true radius (µm)')
    a.legend(fontsize=16, loc='upper left')
ax[0].set_ylabel('estimated radius (µm)')
plt.tight_layout(w_pad=3)
fig.savefig(os.path.join(FIG_DIR, 'fp_radius_tissue_only.png'), dpi=110)
plt.close(fig)

summary = {f'{part}|{name}': {'radii_right': int((r['estR'] == radii).sum()),
                               'estR': [float(x) for x in r['estR']],
                               'so2_err_points': [round(float(x), 1) for x in r['so2err']]}
           for (part, name), r in result.items()}
with open(os.path.join(FIG_DIR, 'tissue_only_summary.json'), 'w') as f:
    json.dump(summary, f, indent=1)
for k, v in summary.items():
    print(f'{k:50s} radii right {v["radii_right"]:2d}/12   SO2 error (points) {v["so2_err_points"]}')

"""Analyse the playground runs and draw the figures for the talk.

Usage: python analyze_play.py [RUN_DIR] [FIG_DIR]   (defaults: out, figs)
"""
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bss_common as bc

RUN_DIR = sys.argv[1] if len(sys.argv) > 1 else 'out'
FIG_DIR = sys.argv[2] if len(sys.argv) > 2 else 'figs'
_n_runs = len(glob.glob(os.path.join(RUN_DIR, '*.npz')))
if _n_runs < 96:
    sys.exit(f'Found {_n_runs} of 96 runs in "{RUN_DIR}". Run (or finish) python run_play.py first.')
os.makedirs(FIG_DIR, exist_ok=True)

INK, MUTED, HAIR = '#1B1B1F', '#5B5F66', '#D9DADC'
RED, BLUE, GRAY, GREEN = '#A8322D', '#1F4E79', '#8A8F96', '#3D7A4F'
plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['Liberation Serif', 'Times New Roman', 'Times', 'DejaVu Serif'], 'mathtext.fontset': 'stix',
    'font.size': 20, 'axes.edgecolor': MUTED, 'axes.linewidth': 1.0, 'axes.labelcolor': INK,
    'xtick.color': MUTED, 'ytick.color': MUTED, 'xtick.labelcolor': INK, 'ytick.labelcolor': INK,
    'xtick.labelsize': 18, 'ytick.labelsize': 18, 'axes.labelsize': 21, 'axes.titlesize': 22,
    'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
})

# ---------------- load ----------------
runs = {}
for p in glob.glob(os.path.join(RUN_DIR, '*.npz')):
    z = np.load(p)
    runs[(str(z['approach']), float(z['radius_um']), int(z['seed']))] = z
RADII = sorted({k[1] for k in runs})
SO2 = np.load(next(iter(glob.glob(os.path.join(RUN_DIR, '*.npz')))))['so2']
TE = bc.TE


def mean_signal(approach, seeds, part='tot'):
    """(n_R, n_SO2, 38) signal averaged over seeds."""
    return np.stack([np.mean([runs[(approach, r, s)][part] for s in seeds], 0) for r in RADII])


ref_all = {p: mean_signal('3d', [1, 2, 3, 4], p) for p in ('tot', 'ev', 'iv')}
test = mean_signal('3d', [1, 2])
dicts = {
    '3D, same recipe as the test': mean_signal('3d', [3, 4]),
    '2D, own field direction per vessel': mean_signal('random', [11, 12]),
    '2D, 3B0 average (MrVox-style)': mean_signal('3b0', [11, 12]),
}
colors = {'3D, same recipe as the test': GRAY, '2D, own field direction per vessel': BLUE,
          '2D, 3B0 average (MrVox-style)': RED}


def match(D, S):
    """Normalized inner product matching. D: (nR, nSO2, 38); S: (..., 38) -> (iR, iSO2) per signal."""
    Dn = D.reshape(-1, D.shape[-1])
    Dn = Dn / np.linalg.norm(Dn, axis=1, keepdims=True)
    Sn = S.reshape(-1, S.shape[-1])
    Sn = Sn / np.linalg.norm(Sn, axis=1, keepdims=True)
    best = np.argmax(Sn @ Dn.T, axis=1)
    iR, iS = np.unravel_index(best, D.shape[:2])
    return iR.reshape(S.shape[:-1]), iS.reshape(S.shape[:-1])


SO2_TRUE = [0.5, 0.6, 0.7, 0.8]
iso = [int(np.argmin(np.abs(SO2 - s))) for s in SO2_TRUE]
radii = np.array(RADII)
results = {}
for name, D in dicts.items():
    S = test[:, iso, :]                       # (nR, 4, 38)
    iR, iS = match(D, S)
    estR, estS = radii[iR], SO2[iS]
    results[name] = dict(estR=estR, estSO2=estS)

# noisy version: SNR 100 on the first echo, 200 noise draws, joint matching
rng = np.random.default_rng(0)
SNR = 100
noisy = {}
for name, D in dicts.items():
    S = test[:, iso[1], :]                    # SO2 = 0.6
    sig = S[:, :1] / SNR
    draws = []
    for _ in range(200):
        n = np.abs(S + sig * (rng.standard_normal(S.shape) + 1j * rng.standard_normal(S.shape)))
        draws.append(n)
    draws = np.stack(draws, 1)                # (nR, 200, 38)
    iR, iS = match(D, draws)
    noisy[name] = dict(medR=np.median(radii[iR], 1), q1R=np.percentile(radii[iR], 25, 1),
                       q3R=np.percentile(radii[iR], 75, 1), medS=np.median(SO2[iS], 1))

summary = {'radii_um': RADII, 'so2_true': SO2_TRUE}
for name, r in results.items():
    err = np.abs(r['estR'] - radii[:, None]) / radii[:, None] * 100
    wrong = r['estR'] != radii[:, None]
    summary[name] = {
        'radius_wrong_at_so2_0.6': int(wrong[:, 1].sum()),
        'radius_mean_abs_err_pct_at_so2_0.6': float(err[:, 1].mean()),
        'radius_max_abs_err_pct_at_so2_0.6': float(err[:, 1].max()),
        'radius_wrong_all_so2': int(wrong.sum()),
        'so2_err_points_at_0.6': [round(float(x), 3) for x in (r['estSO2'][:, 1] - 0.6) * 100],
        'so2_mean_abs_err_points_all': float(np.abs(r['estSO2'] - np.array(SO2_TRUE)[None]).mean() * 100),
        'estR_at_so2_0.6': [float(x) for x in r['estR'][:, 1]],
        'noisy_snr100_median_R': [float(x) for x in noisy[name]['medR']],
        'noisy_snr100_median_SO2': [float(x) for x in noisy[name]['medS']],
    }

# Boxerman-style rates at SO2 0.6 (T2 off): GE at 29.4 ms, SE at 100 ms
i6 = iso[1]
def rates(sig):
    return (-np.log(sig[:, i6, bc.GE_IDX]) / (TE[bc.GE_IDX] / 1000),
            -np.log(sig[:, i6, bc.SE_IDX]) / (TE[bc.SE_IDX] / 1000))
box = {'3D reference': {p: rates(ref_all[p]) for p in ('tot', 'ev', 'iv')}}
for lab, a in (('2D, own direction', 'random'), ('2D, 3B0 average', '3b0')):
    box[lab] = {p: rates(mean_signal(a, [11, 12], p)) for p in ('tot', 'ev', 'iv')}
summary['boxerman'] = {k: {p: [[round(float(x), 3) for x in v[0]], [round(float(x), 3) for x in v[1]]]
                           for p, v in d.items()} for k, d in box.items()}
summary['iv_fraction'] = {a: float(np.mean([z['iv_fraction'] for k, z in runs.items() if k[0] == a]))
                          for a in ('3d', 'random', '3b0')}
with open(os.path.join(FIG_DIR, 'summary.json'), 'w') as f:
    json.dump(summary, f, indent=1)

# ---------------- figure 1: radius and SO2 recovered by each dictionary ----------------
fig, ax = plt.subplots(1, 3, figsize=(16, 5.6), sharex=True, sharey=True)
for a, (name, r) in zip(ax, results.items()):
    a.plot([0.8, 65], [0.8, 65], color=HAIR, lw=2, zorder=1)
    a.plot(radii, r['estR'][:, 1], 'o', ms=11, mfc='white', mec=colors[name], mew=2.4, zorder=3)
    a.set_xscale('log'); a.set_yscale('log')
    a.set_xticks([1, 3, 10, 30]); a.set_xticklabels(['1', '3', '10', '30'])
    a.set_yticks([1, 3, 10, 30]); a.set_yticklabels(['1', '3', '10', '30'])
    a.set_xlim(0.8, 65); a.set_ylim(0.8, 65)
    a.set_title(name, color=colors[name] if colors[name] != GRAY else INK, pad=14)
    a.set_xlabel('true radius (µm)')
    n_wrong = summary[name]['radius_wrong_at_so2_0.6']
    a.text(0.04, 0.95, f'{len(radii) - n_wrong} of {len(radii)} radii right', transform=a.transAxes,
           va='top', fontsize=19, color=MUTED)
ax[0].set_ylabel('estimated radius (µm)')
plt.tight_layout(w_pad=2.5)
fig.savefig(os.path.join(FIG_DIR, 'fp_radius.png'), dpi=110)
plt.close(fig)

# ---------------- figure 2: SO2 error vs radius ----------------
fig, a = plt.subplots(figsize=(10, 5.6))
a.axhline(0, color=HAIR, lw=2)
for name, r in results.items():
    e = (r['estSO2'] - np.array(SO2_TRUE)[None]) * 100        # points, (nR, 4)
    a.plot(radii, e[:, 1], 'o-', ms=9, lw=2, color=colors[name], mfc='white', mew=2.2, label=name)
a.set_xscale('log'); a.set_xticks([1, 3, 10, 30]); a.set_xticklabels(['1', '3', '10', '30'])
a.set_xlabel('true radius (µm)'); a.set_ylabel('SO$_2$ error (points)')
a.legend(fontsize=17, loc='best')
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'fp_so2.png'), dpi=110)
plt.close(fig)

# ---------------- figure 3: Boxerman-style plot on GESFIDE timing ----------------
fig, ax = plt.subplots(1, 2, figsize=(16, 5.6))
style = {'3D reference': dict(color=INK, marker='o', ls='-'),
         '2D, own direction': dict(color=BLUE, marker='s', ls='--'),
         '2D, 3B0 average': dict(color=RED, marker='^', ls=':')}
for k, d in box.items():
    ge, se = d['tot']
    ax[0].plot(radii, ge, ms=8, lw=2, mfc='white', mew=2, label=k, **style[k])
    ax[1].plot(radii, se, ms=8, lw=2, mfc='white', mew=2, label=k, **style[k])
for a, t in zip(ax, ['gradient echo, TE 29 ms: ΔR$_2$*', 'spin echo, TE 100 ms: ΔR$_2$']):
    a.set_xscale('log'); a.set_xticks([1, 3, 10, 30]); a.set_xticklabels(['1', '3', '10', '30'])
    a.set_xlabel('vessel radius (µm)'); a.set_title(t, pad=12)
ax[0].set_ylabel('relaxation rate (1/s)')
ax[1].legend(fontsize=17)
plt.tight_layout(w_pad=3)
fig.savefig(os.path.join(FIG_DIR, 'boxerman_gesfide.png'), dpi=110)
plt.close(fig)

# ---------------- figure 4: what the knobs do (3D reference, T2 off) ----------------
def gap(y):
    g = int(np.searchsorted(TE, 50.0))
    return np.insert(TE, g, np.nan), np.insert(y, g, np.nan)

fig, ax = plt.subplots(1, 2, figsize=(16, 5.6), sharey=True)
iR5 = RADII.index(5.0)
so2_show = [0.5, 0.6, 0.7, 0.8, 0.9]
cm = plt.get_cmap('Blues')
for j, s in enumerate(so2_show):
    y = ref_all['tot'][iR5, int(np.argmin(np.abs(SO2 - s)))]
    ax[0].plot(*gap(y), '.-', ms=6, lw=2, color=cm(0.35 + 0.15 * j), label=f'SO$_2$ {int(s * 100)}%')
r_show = [2.0, 5.0, 10.0, 30.0]
cr = plt.get_cmap('Reds')
for j, r in enumerate(r_show):
    y = ref_all['tot'][RADII.index(r), i6]
    ax[1].plot(*gap(y), '.-', ms=6, lw=2, color=cr(0.35 + 0.18 * j), label=f'R {r:g} µm')
for a in ax:
    a.axvline(50, color=HAIR, lw=1.5, ls='--'); a.axvline(100, color=HAIR, lw=1.5, ls=':')
    a.set_xlabel('echo time (ms)')
ax[0].set_ylabel('signal (T2 off)')
ax[0].set_title('Change oxygenation (R 5 µm, CBV 3%)', pad=12)
ax[1].set_title('Change vessel size (SO$_2$ 60%, CBV 3%)', pad=12)
ax[0].legend(fontsize=16, loc='lower left'); ax[1].legend(fontsize=16, loc='lower left')
plt.tight_layout(w_pad=2.5)
fig.savefig(os.path.join(FIG_DIR, 'knobs.png'), dpi=110)
plt.close(fig)

# ---------------- figure 5: field maps, own direction vs 3B0 ----------------
base = bc.voxel2d(11)
g = bc.grids2d(base, 512)
fig, ax = plt.subplots(1, 2, figsize=(14, 6.4))
lim = 60
for a, key, t in zip(ax, ('random', '3b0'), ('each vessel, its own field direction', '3B0: one averaged field for all')):
    idx, dbz, _ = g[key]
    sub = dbz[:200, :200] * 1e9
    im = a.imshow(sub.T, origin='lower', cmap='RdBu_r', vmin=-lim, vmax=lim)
    m = idx[:200, :200].T > 0
    a.contour(m, levels=[0.5], colors=INK, linewidths=0.8)
    a.set_xticks([]); a.set_yticks([])
    for s in a.spines.values():
        s.set_visible(False)
    a.set_title(t, pad=12)
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
cb.set_label('field offset (nT)')
fig.savefig(os.path.join(FIG_DIR, 'fieldmaps.png'), dpi=110, bbox_inches='tight')
plt.close(fig)

# IV field spread
for key in ('random', '3b0'):
    idx, dbz, _ = g[key]
    v = dbz[idx > 0] * 1e9
    summary.setdefault('iv_field_nT_p10_p50_p90', {})[key] = [round(float(x), 1) for x in np.percentile(v, [10, 50, 90])]
with open(os.path.join(FIG_DIR, 'summary.json'), 'w') as f:
    json.dump(summary, f, indent=1)
names = list(results)
right = [len(radii) - summary[n]['radius_wrong_at_so2_0.6'] for n in names]
print(f'Saved 5 figures and summary.json in the folder "{FIG_DIR}":')
print(f'  fieldmaps.png         slide 8   field maps, own direction per vessel vs 3B0')
print(f'  knobs.png             slide 11  GESFIDE signal as SO2 and vessel radius change')
print(f'  boxerman_gesfide.png  slide 12  decay rate vs vessel radius, 3D vs 2D vs 3B0')
print(f'  fp_radius.png         slide 13  radii right at SO2 60%: 3D {right[0]}/12, 2D own direction {right[1]}/12, 3B0 {right[2]}/12')
print(f'  fp_so2.png            slide 14  SO2 error vs radius for the same matches')
print('Next: python compare_to_mine.py')

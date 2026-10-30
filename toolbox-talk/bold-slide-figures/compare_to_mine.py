"""Compare your playground results (figs/summary.json) with the reference run (reference_summary.json).

Usage: python compare_to_mine.py [FIG_DIR]   (default: figs)

The simulations use fixed random seeds, so a working install gives the same
numbers. Tiny differences in the last digit can come from a different processor;
the counts of radii right should match exactly.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main(fig_dir='figs'):
    """Print the comparison table; return True when every check matches."""
    mine_path = os.path.join(fig_dir, 'summary.json')
    if not os.path.exists(mine_path):
        print(f'No {mine_path} yet. Run the figure cells (or python analyze_play.py) first.')
        return False
    yours = json.load(open(mine_path))
    ref = json.load(open(os.path.join(HERE, 'reference_summary.json')))

    DICTS = ['3D, same recipe as the test', '2D, own field direction per vessel', '2D, 3B0 average (MrVox-style)']
    SHORT = {DICTS[0]: '3D', DICTS[1]: '2D own direction', DICTS[2]: '3B0'}
    rows = []

    def check(label, a, b, tol):
        ok = abs(a - b) <= tol
        rows.append((label, a, b, ok))

    for d in DICTS:
        check(f'{SHORT[d]}: radii right at SO2 60% (of 12)',
              12 - yours[d]['radius_wrong_at_so2_0.6'], 12 - ref[d]['radius_wrong_at_so2_0.6'], 0)
    for d in DICTS:
        same = sum(x == y for x, y in zip(yours[d]['estR_at_so2_0.6'], ref[d]['estR_at_so2_0.6']))
        check(f'{SHORT[d]}: estimated radii identical (of 12)', same, 12, 0)
    for d in DICTS:
        worst = max(abs(x - y) for x, y in zip(yours[d]['so2_err_points_at_0.6'], ref[d]['so2_err_points_at_0.6']))
        check(f'{SHORT[d]}: largest SO2 difference (points)', worst, 0, 1)
    radii = ref['radii_um']
    for k in ['3D reference', '2D, own direction', '2D, 3B0 average']:
        for i_r in (radii.index(5.0), radii.index(50.0)):
            check(f'{k}: GE rate at {radii[i_r]:g} um (1/s)',
                  yours['boxerman'][k]['tot'][0][i_r], ref['boxerman'][k]['tot'][0][i_r], 0.05)
    check('3D reference: blood GE rate at 5 um (1/s)',
          yours['boxerman']['3D reference']['iv'][0][radii.index(5.0)],
          ref['boxerman']['3D reference']['iv'][0][radii.index(5.0)], 0.5)

    width = max(len(r[0]) for r in rows)
    print(f'{"check":{width}}  {"yours":>8}  {"mine":>8}')
    for label, a, b, ok in rows:
        print(f'{label:{width}}  {a:8.3g}  {b:8.3g}  {"match" if ok else "DIFFERENT"}')
    n_bad = sum(not r[3] for r in rows)
    if n_bad == 0:
        print(f'\nAll {len(rows)} checks match. Your figures are the ones on the slides.')
    else:
        print(f'\n{n_bad} of {len(rows)} checks differ. See "If something goes wrong" in the guide.')
    return n_bad == 0


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'figs')

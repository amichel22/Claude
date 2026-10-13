"""Shared setup for the BOLDswimsuite playground runs.

GESFIDE timing, physiology and the phase-history helpers follow Jack's
38-echo walkthrough notebook (sections 0, 4 and 7), so the numbers line up
with what you see there.
"""
import numpy as np
from boldswimsuite import BOLDvessel, BOLDgeometry, BOLDspins

# ---------------- 38-echo GESFIDE timing (Echotimes_38.mat, as in the notebook) ----------------
# part A: 13 echoes before the 180 at 50 ms; parts B + C: 25 echoes, spin echo at 100 ms.
# Echo 35 is 121.360 ms in the lab file (0.024 ms off the regular 3.048 ms train); kept as is.
TE = np.concatenate([4.988 + 3.048 * np.arange(13), 57.328 + 3.048 * np.arange(25)])
TE[34] = 121.360
TE = np.round(TE, 3)
T180 = 50.0
SE_IDX = int(np.argmin(np.abs(TE - 100.0)))       # 27
GE_IDX = int(np.argmin(np.abs(TE - 30.0)))        # 8 -> 29.372 ms, Berman's GE TE is 30 ms

DT = 0.2
events = np.unique(np.concatenate([[0.0, T180], TE]))
edges = [0.0]
for a_, b_ in zip(events[:-1], events[1:]):
    n = int(np.ceil((b_ - a_) / DT - 1e-9))
    edges += list(a_ + (b_ - a_) * np.arange(1, n + 1) / n)
edges = np.array(edges)
STEP_DT = np.diff(edges)
T_GRID = edges[1:]
NSTEPS = len(STEP_DT)
ECHO_IDX = np.array([int(np.argmin(np.abs(T_GRID - te))) for te in TE])
assert np.allclose(T_GRID[ECHO_IDX], TE, atol=1e-9)
PULSE_180 = int(np.argmin(np.abs(edges - T180)))

# ---------------- physiology ----------------
B0, ADC, HCT = 3.0, 1e-3, 0.40
DCHI0_CGS = 0.27e-6                  # per unit Hct, cgs (Spees 2001), as in the notebook
SO2_REF = 0.60                       # histories are recorded at this SO2 and rescaled
CBV = 0.03
R0_UM = 5.0                          # geometry is built at 5 um and scaled to other radii
N_VESSELS = 150


def dchi_cgs(so2, hct=HCT):
    return DCHI0_CGS * hct * (1 - so2)


def phase_history(voxel, nspins, seed, ndims):
    """One MC run -> accumulated (sign-flipped) phase of every spin at the 38 echoes + IV flags."""
    cls = BOLDspins.Spins3D if ndims == 3 else BOLDspins.Spins2D
    S = cls(num_spins=nspins, geometry=voxel, dt=DT, IV=True, seed=seed)
    phi = np.zeros(nspins)
    out = np.zeros((len(TE), nspins))
    k = 0
    echo_set = set(ECHO_IDX.tolist())
    for j, dt in enumerate(STEP_DT):
        S.step(dt)
        if j == PULSE_180:
            phi = -phi
        phi += S.phase
        if j in echo_set:
            out[k] = phi
            k += 1
    return out, S.vessel_indices > 0


def signals_from_history(phi, iv, so2_grid):
    """Magnitude signals (T2 off: susceptibility only) for each SO2: arrays (n_so2, 38) for total, EV, IV."""
    tot, ev, iv_s = [], [], []
    for s in so2_grid:
        M = np.exp(1j * (dchi_cgs(s) / dchi_cgs(SO2_REF)) * phi)
        tot.append(np.abs(M.mean(1)))
        ev.append(np.abs(M[:, ~iv].mean(1)))
        iv_s.append(np.abs(M[:, iv].mean(1)))
    return np.array(tot), np.array(ev), np.array(iv_s)


# ---------------- geometries ----------------
def voxel3d(seed):
    return BOLDgeometry.ContinuousVoxel3D.from_random(
        CBV=CBV, B0=B0, ADC=ADC, labels=['v'], weights={'v': 1},
        diameter_distributions={'v': [2 * R0_UM * 1e-3]}, dchis={'v': dchi_cgs(SO2_REF)},
        permeation_probabilities={'v': 0.0}, num_vessels=N_VESSELS, vessel_type='cylinder',
        allow_vessel_intersection=True, seed=seed, progressbar=False)


def scaled3d(vox, a):
    return BOLDgeometry.ContinuousVoxel3D(
        vox.size.value * a, vox.B0, vox.ADC,
        [BOLDvessel.InfiniteCylinder3D(v.diameter * a, v.theta, v.phi, v.origin * a, v.dchi) for v in vox.vessels])


def voxel2d(seed):
    return BOLDgeometry.ContinuousVoxel2D.from_random(
        CBV=CBV, B0=B0, ADC=ADC, labels=['v'], weights={'v': 1},
        diameter_distributions={'v': [2 * R0_UM * 1e-3]}, dchis={'v': dchi_cgs(SO2_REF)},
        num_vessels=N_VESSELS, seed=seed, progressbar=False)


def grids2d(vox, n_grid):
    """Two gridded versions of the same 2D geometry:
    'random': each vessel keeps its own random B0 direction (Berman's 2D approach);
    '3b0'   : MrVox-style field, 1/3 * (B0 along the vessels) + 2/3 * (B0 perpendicular)."""
    def with_b0(theta):
        return BOLDgeometry.ContinuousVoxel2D(
            size=vox.size.value, B0=vox.B0, ADC=ADC,
            vessels=[BOLDvessel.InfiniteCylinder2D(v.diameter, theta, 0.0, v.origin.copy(), v.dchi)
                     for v in vox.vessels])
    g_rand = BOLDgeometry.DiscreteVoxel2D.from_continuous_analytical(n_grid, vox)
    g_par = BOLDgeometry.DiscreteVoxel2D.from_continuous_analytical(n_grid, with_b0(0.0))
    g_perp = BOLDgeometry.DiscreteVoxel2D.from_continuous_analytical(n_grid, with_b0(np.pi / 2))
    dbz_3b0 = g_par.dBz_grid / 3 + 2 * g_perp.dBz_grid / 3
    return {'random': (g_rand.vessel_index_grid, g_rand.dBz_grid, g_rand.permeation_probability_list),
            '3b0': (g_par.vessel_index_grid, dbz_3b0, g_par.permeation_probability_list)}


def grid_voxel(grid, size):
    idx, dbz, perm = grid
    return BOLDgeometry.DiscreteVoxel2D(idx, dbz, ADC=ADC, permeation_probability_list=perm, size=size)

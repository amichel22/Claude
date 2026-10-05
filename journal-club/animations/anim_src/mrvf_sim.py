"""mrvf_sim.py -- 2D Monte Carlo simulator for the 38-echo GESFIDE MRvF signal.

Model (BOLDswimsuite-style "2D trick", Chausse/Berman):
  * Square periodic voxel in the plane; infinite cylinders run perpendicular to
    the plane.  Every cylinder gets its OWN B0 direction: polar angle theta with
    cos(theta) uniform on [0, 1] (equivalent to [-1, 1]: only cos^2 enters) and
    an in-plane azimuth phi_B0 uniform on [0, 2*pi).
  * Field offset (Ogawa 1993, Lorentz-corrected inside):
        r >  R : dBz = B0*dchi_SI/2 * (R/r)^2 * cos(2(phi - phi_B0)) * sin^2(theta)
        r <= R : dBz = B0*dchi_SI/6 * (3 cos^2(theta) - 1)          (uniform, IV)
    The field of all vessels is summed (minimum-image periodic convention).
  * dchi_SI = 4*pi * 0.27e-6 * Hct * (1 - SO2)   (cgs 0.27e-6 per unit Hct, Spees 2001).
  * Water diffuses with D = 1 um^2/ms in both compartments.  Walls are
    impermeable: a step that would cross a wall is re-proposed (up to n_retry
    times, fresh Gaussian draws), otherwise the spin stays where it is.
  * Time grid: every interval between events (0, each echo, the 180 at 50 ms)
    is split into equal steps <= dt_max, so all 38 echoes and the 180 land
    exactly on step boundaries.  Phase uses the trapezoid rule along the path.
  * Phase-history trick: the walk does not depend on dchi, so the accumulated
    (180-sign-flipped) phase of every spin is recorded at the echoes ONCE, at a
    reference state (SO2 = 0, Hct = 0.4, B0 = 3 T).  Any SO2 / Hct / B0 is then
        S = | < w_k * exp(i * s * Phi_k) > |,   s = (1-SO2)*Hct*B0 / (1*0.4*3).
    T2 is applied afterwards per compartment (w = exp(-TE/T2_EV) or exp(-TE/T2_IV)).
  * Everything geometric is done in units of the vessel radius R (R = 1), so one
    geometry serves every radius; only the diffusion step sqrt(2 D dt)/R changes.

Speed: the field is split into a smooth "far" part, pre-computed once per
geometry on a periodic grid (spacing R/4) and bilinearly interpolated, plus an
exact analytic "near" correction for the few vessels within r_c = 4R of a spin
(found through a cell table).  Inside r_c a vessel's far part is replaced by the
C1-continuous polynomial r^2 cos2phi (3/rc^4 - 2 r^2/rc^6), so the grid only
sees smooth fields.  The near correction (exact - polynomial) also gives the
compartment test (r < R) used for the impermeable walls.

Main entry points
  make_schedule(...)                 -> step lengths, echo / 180 step indices
  make_voxel(CBV, seed, ...)         -> VoxelGeometry2D (CBV-targeted, nested in CBV)
  nested_voxels(cbv_list, seed, ...) -> dict CBV -> VoxelGeometry2D (shared prefix)
  phase_history(geom, R_um, ...)     -> dict(Phi (38, Ns) [rad at reference], iv, ...)
  combine_histories([...])           -> pooled history (e.g. several seeds)
  signals_from_history(H, SO2, ...)  -> EV, IV, total signals (complex + magnitude)
  simulate_gesfide(R_um, CBV, SO2)   -> convenience wrapper

Units: lengths um (geometry stored in units of R), time ms, D um^2/ms, B0 T.
"""
from __future__ import annotations

import time
import numpy as np

# ----------------------------------------------------------------- constants
GAMMA_HZ_PER_T = 42.58e6
GAMMA = 2.0 * np.pi * GAMMA_HZ_PER_T          # rad / s / T
B0_T = 3.0
HCT = 0.40
DCHI0_CGS = 0.27e-6                           # per unit Hct, fully deoxygenated blood (cgs)
D_UM2_MS = 1.0
T180_MS = 50.0
DT_MAX_MS = 0.2
T2_EV_MS = 80.0
T2_IV_MS = 50.0                               # blood T2 (simple fixed value)
TE_MS = np.round(np.concatenate([4.988 + 3.048 * np.arange(13),
                                 57.328 + 3.048 * np.arange(25)]), 6)
PART_A = TE_MS < T180_MS
PART_B = (TE_MS > T180_MS) & (TE_MS <= 2 * T180_MS + 1e-9)
PART_C = TE_MS > 2 * T180_MS + 1e-9
SE_INDEX = int(np.argmin(np.abs(TE_MS - 2 * T180_MS)))   # 27 -> TE = 100.000 ms

# geometry discretisation (units of R)
RC_NEAR = 4.0          # exact near-field radius
GRID_H = 0.25          # far-field grid spacing
CELL = 1.0             # cell size for the near-vessel table
VESSELS_PER_PCT = 100  # vessels per 1 % CBV  -> box side L = sqrt(100*pi/0.01) R = 177 R


def dchi_cgs(SO2, Hct=HCT):
    return DCHI0_CGS * Hct * (1.0 - np.asarray(SO2, float))


def dchi_si(SO2, Hct=HCT):
    return 4.0 * np.pi * dchi_cgs(SO2, Hct)


def omega_ref(B0=B0_T, Hct=HCT):
    """gamma*B0*dchi_SI at the reference state SO2 = 0 (rad/ms)."""
    return GAMMA * B0 * float(dchi_si(0.0, Hct)) * 1e-3


OMEGA_REF = omega_ref()   # 1.089 rad/ms


def phase_scale(SO2, Hct=HCT, B0=B0_T):
    """Factor that maps reference-state phases to (SO2, Hct, B0)."""
    return (1.0 - np.asarray(SO2, float)) * (Hct / HCT) * (B0 / B0_T)


def static_dephasing_R2p(CBV, SO2, Hct=HCT, B0=B0_T):
    """Yablonskiy/Haacke long-time static-dephasing R2' (1/s), as written in the notebook:
    CBV * 2*pi*42.58e6 * (4/3)*pi * dchi_cgs * B0."""
    return CBV * GAMMA * (4.0 / 3.0) * np.pi * dchi_cgs(SO2, Hct) * B0


# ------------------------------------------------------------------ schedule
def make_schedule(TE=TE_MS, t180=T180_MS, dt_max=DT_MAX_MS):
    TE = np.asarray(TE, float)
    events = np.unique(np.concatenate([[0.0, t180], TE]))
    edges = [0.0]
    for a, b in zip(events[:-1], events[1:]):
        n = int(np.ceil((b - a) / dt_max - 1e-9))
        edges += list(a + (b - a) * np.arange(1, n + 1) / n)
    edges = np.asarray(edges)
    dt = np.diff(edges)
    t_end = edges[1:]
    echo_step = np.array([int(np.argmin(np.abs(t_end - te))) for te in TE])
    assert np.allclose(t_end[echo_step], TE, atol=1e-9)
    pulse_step = int(np.argmin(np.abs(t_end - t180)))
    assert abs(t_end[pulse_step] - t180) < 1e-9
    return dict(dt=dt, t=t_end, echo_step=echo_step, pulse_step=pulse_step, TE=TE)


# ------------------------------------------------------------------ geometry
def _place_vessels(n, L, rng, min_sep=2.0, max_tries=200000):
    """Random sequential addition of non-overlapping unit disks in a periodic box."""
    xy = np.empty((n, 2))
    k, tries = 0, 0
    while k < n:
        p = rng.random(2) * L
        if k:
            d = xy[:k] - p
            d -= L * np.rint(d / L)
            if np.min(d[:, 0] ** 2 + d[:, 1] ** 2) < min_sep ** 2:
                tries += 1
                if tries > max_tries:
                    raise RuntimeError("could not place vessels (CBV too high?)")
                continue
        xy[k] = p
        k += 1
    return xy


def _far_field_grids(L, xy, A, c2, s2, prefixes, rc=RC_NEAR, h=GRID_H):
    """Smooth 'far' field sum on a periodic Ng x Ng grid, for each vessel-count prefix."""
    Ng = int(np.ceil(L / h))
    h = L / Ng
    g = np.arange(Ng) * h
    prefixes = sorted(set(int(p) for p in prefixes))
    out, G = {}, np.zeros((Ng, Ng))
    rc2, a4, a6 = rc * rc, 3.0 / rc ** 4, 2.0 / rc ** 6
    with np.errstate(divide="ignore", invalid="ignore"):
        for v in range(prefixes[-1]):
            dx = g[:, None] - xy[v, 0]
            dx -= L * np.rint(dx / L)
            dy = g[None, :] - xy[v, 1]
            dy -= L * np.rint(dy / L)
            dx2, dy2 = dx * dx, dy * dy
            r2 = dx2 + dy2
            q = (dx2 - dy2) * c2[v] + 2.0 * dx * dy * s2[v]
            G += A[v] * q * np.where(r2 >= rc2, 1.0 / (r2 * r2), a4 - a6 * r2)
            if v + 1 in prefixes:
                out[v + 1] = G.copy()
    return out, Ng


def _cell_table(L, xy, rc=RC_NEAR, cell=CELL):
    """For every cell, the vessels whose centre lies within rc of the cell (padded with n)."""
    n = len(xy)
    nc = max(3, int(L // cell))
    c = L / nc
    w = int(np.ceil(rc / c)) + 1
    assert 2 * w + 1 < nc, "voxel too small for the cell table"
    o = np.arange(-w, w + 1)
    ox, oy = [a.ravel() for a in np.meshgrid(o, o, indexing="ij")]
    ix0 = np.floor(xy[:, 0] / c).astype(np.int64)
    iy0 = np.floor(xy[:, 1] / c).astype(np.int64)
    IX = ix0[:, None] + ox[None, :]
    IY = iy0[:, None] + oy[None, :]
    dxr = np.maximum(np.maximum(IX * c - xy[:, :1], 0.0), xy[:, :1] - (IX + 1) * c)
    dyr = np.maximum(np.maximum(IY * c - xy[:, 1:], 0.0), xy[:, 1:] - (IY + 1) * c)
    keep = dxr ** 2 + dyr ** 2 < (rc + 1e-9) ** 2
    vid = np.broadcast_to(np.arange(n)[:, None], IX.shape)[keep]
    cid = (IX[keep] % nc) * nc + (IY[keep] % nc)
    order = np.argsort(cid, kind="stable")
    cid, vid = cid[order], vid[order]
    counts = np.bincount(cid, minlength=nc * nc)
    K = max(1, int(counts.max()))
    start = np.concatenate([[0], np.cumsum(counts)[:-1]])
    rank = np.arange(len(cid)) - start[cid]
    table = np.full((nc * nc, K), n, dtype=np.int64)
    table[cid, rank] = vid
    return table, nc, c


class VoxelGeometry2D:
    """Periodic square voxel, side L (units of R); n unit-radius vessels.

    Field values are in units of OMEGA_REF (= gamma*B0*dchi_SI at SO2 = 0):
      EV amplitude A = sin^2(theta)/2 ; IV offset W = (3 cos^2(theta) - 1)/6.
    """

    def __init__(self, L, xy, cos_theta, phi_b0, grid=None, rc=RC_NEAR, h=GRID_H, cell=CELL):
        self.L = float(L)
        self.xy = np.asarray(xy, float)
        self.n = len(self.xy)
        self.cos_theta = np.asarray(cos_theta, float)
        self.phi_b0 = np.asarray(phi_b0, float)
        self.rc = float(rc)
        A = 0.5 * (1.0 - self.cos_theta ** 2)
        W = (3.0 * self.cos_theta ** 2 - 1.0) / 6.0
        c2, s2 = np.cos(2 * self.phi_b0), np.sin(2 * self.phi_b0)
        self.A, self.W, self.c2, self.s2 = A, W, c2, s2
        # padded copies (index n = dummy vessel with NaN position, zero field)
        self._x = np.append(self.xy[:, 0], np.nan)
        self._y = np.append(self.xy[:, 1], np.nan)
        self._A, self._W = np.append(A, 0.0), np.append(W, 0.0)
        self._c2, self._s2 = np.append(c2, 0.0), np.append(s2, 0.0)
        if grid is None:
            grids, Ng = _far_field_grids(self.L, self.xy, A, c2, s2, [self.n], rc, h)
            grid = grids[self.n]
        self.G = np.ascontiguousarray(grid)
        self.Ng = self.G.shape[0]
        self.h = self.L / self.Ng
        self._Gf = self.G.ravel()
        self.table, self.nc, self.cell = _cell_table(self.L, self.xy, rc, cell)
        self.K = self.table.shape[1]

    @property
    def cbv(self):
        return self.n * np.pi / self.L ** 2

    # ---- field and compartment at arbitrary points (units of R) ----
    def _interp(self, P):
        fx, fy = P[:, 0] / self.h, P[:, 1] / self.h
        i0, j0 = np.floor(fx).astype(np.int64), np.floor(fy).astype(np.int64)
        tx, ty = fx - i0, fy - j0
        Ng = self.Ng
        i0 %= Ng
        j0 %= Ng
        i1, j1 = (i0 + 1) % Ng, (j0 + 1) % Ng
        Gf = self._Gf
        return ((1 - tx) * ((1 - ty) * Gf[i0 * Ng + j0] + ty * Gf[i0 * Ng + j1])
                + tx * ((1 - ty) * Gf[i1 * Ng + j0] + ty * Gf[i1 * Ng + j1]))

    def evaluate(self, P):
        """Return (field [units of OMEGA_REF], vessel index or -1) at points P (n, 2)."""
        L, nc = self.L, self.nc
        ix = (P[:, 0] // self.cell).astype(np.int64) % nc
        iy = (P[:, 1] // self.cell).astype(np.int64) % nc
        cand = self.table[ix * nc + iy]
        with np.errstate(divide="ignore", invalid="ignore"):
            dx = P[:, :1] - self._x[cand]
            dx -= L * np.rint(dx / L)
            dy = P[:, 1:] - self._y[cand]
            dy -= L * np.rint(dy / L)
            dx2, dy2 = dx * dx, dy * dy
            r2 = dx2 + dy2
            ins = r2 < 1.0
            near = r2 < self.rc * self.rc
            q = self._A[cand] * ((dx2 - dy2) * self._c2[cand] + 2.0 * dx * dy * self._s2[cand])
            exact = np.where(ins, self._W[cand], q / (r2 * r2))
            inner = q * (3.0 / self.rc ** 4 - 2.0 * r2 / self.rc ** 6)
            corr = np.where(near, exact - inner, 0.0).sum(1)
        any_in = ins.any(1)
        idx = np.where(any_in, cand[np.arange(len(P)), ins.argmax(1)], -1)
        return self._interp(P) + corr, idx

    def field_direct(self, P):
        """Brute-force analytic sum over all vessels (minimum image) -- for validation."""
        f = np.zeros(len(P))
        for v in range(self.n):
            dx = P[:, 0] - self.xy[v, 0]
            dx -= self.L * np.rint(dx / self.L)
            dy = P[:, 1] - self.xy[v, 1]
            dy -= self.L * np.rint(dy / self.L)
            r2 = dx * dx + dy * dy
            with np.errstate(divide="ignore", invalid="ignore"):
                ev = self.A[v] * ((dx * dx - dy * dy) * self.c2[v] + 2 * dx * dy * self.s2[v]) / (r2 * r2)
            f += np.where(r2 < 1.0, self.W[v], ev)
        return f


def _orientations(n, seed, orientation="golden"):
    rng = np.random.default_rng([int(seed), 7])
    u0 = rng.random()
    if orientation == "golden":       # low-discrepancy cos(theta), still uniform on [0, 1]
        cos_t = (u0 + 0.6180339887498949 * np.arange(n)) % 1.0
    else:
        cos_t = rng.random(n)
    phi = np.random.default_rng([int(seed), 8]).random(n) * 2 * np.pi
    return cos_t, phi


def box_side(vessels_per_pct=VESSELS_PER_PCT):
    """Box side (units of R) that makes n vessels give CBV = n / vessels_per_pct  %."""
    return np.sqrt(vessels_per_pct * np.pi / 0.01)


def nested_voxels(cbv_list, seed=0, vessels_per_pct=VESSELS_PER_PCT, orientation="golden",
                  verbose=False):
    """CBV-targeted random voxels sharing one box and one vessel sequence: the voxel for a
    given CBV holds the first n = CBV/1% * vessels_per_pct vessels (non-overlapping, random
    sequential addition).  Returns {CBV: VoxelGeometry2D}."""
    L = box_side(vessels_per_pct)
    ns = {c: int(round(c / 0.01 * vessels_per_pct)) for c in cbv_list}
    nmax = max(ns.values())
    t0 = time.time()
    xy = _place_vessels(nmax, L, np.random.default_rng([int(seed), 5]))
    cos_t, phi = _orientations(nmax, seed, orientation)
    A = 0.5 * (1 - cos_t ** 2)
    grids, _ = _far_field_grids(L, xy, A, np.cos(2 * phi), np.sin(2 * phi), list(ns.values()))
    out = {c: VoxelGeometry2D(L, xy[:n], cos_t[:n], phi[:n], grid=grids[n]) for c, n in ns.items()}
    if verbose:
        print(f"  seed {seed}: {nmax} vessels, L = {L:.1f} R, grids for {len(ns)} CBVs in {time.time()-t0:.1f} s")
    return out


def make_voxel(CBV, seed=0, **kw):
    return nested_voxels([CBV], seed, **kw)[CBV]


# ------------------------------------------------------------- Monte Carlo
def phase_history(geom, R_um, n_spins=20000, seed=0, D=D_UM2_MS, dt_max=DT_MAX_MS,
                  TE=TE_MS, t180=T180_MS, static=False, n_retry=8, record_steps=False):
    """Random walk of n_spins in `geom` (scaled to vessel radius R_um).

    Returns dict with
      Phi (nTE, Ns): accumulated phase (rad) at each echo, reference state (SO2=0,
                     Hct=0.4, B0=3 T), sign flipped by the 180 at t180;
      iv (Ns,) bool, vessel (Ns,) int, n_rejected (total final rejections), runtime.
    The same `seed` gives the same initial positions and main Gaussian draws for any
    R_um and any nested CBV voxel (common random numbers -> smooth sweeps).
    """
    t0 = time.time()
    sch = make_schedule(TE, t180, dt_max)
    rng = np.random.default_rng([int(seed), 1])
    rng_retry = np.random.default_rng([int(seed), 2])
    L = geom.L
    P = rng.random((n_spins, 2)) * L
    f, own = geom.evaluate(P)
    Phi = np.zeros(n_spins)
    out = np.zeros((len(sch["TE"]), n_spins))
    echo_of_step = {int(s): k for k, s in enumerate(sch["echo_step"])}
    n_stuck = 0
    trace = [] if record_steps else None
    for j, dt in enumerate(sch["dt"]):
        if static:
            fn = f
        else:
            sig = np.sqrt(2.0 * D * dt) / R_um
            Q = P + sig * rng.standard_normal((n_spins, 2))
            Q %= L
            fn, idx = geom.evaluate(Q)
            bad = idx != own
            tries = 0
            while tries < n_retry and bad.any():
                b = np.flatnonzero(bad)
                Qb = (P[b] + sig * rng_retry.standard_normal((len(b), 2))) % L
                fb, ib = geom.evaluate(Qb)
                ok = ib == own[b]
                g = b[ok]
                Q[g], fn[g] = Qb[ok], fb[ok]
                bad[g] = False
                tries += 1
            if bad.any():
                n_stuck += int(bad.sum())
                Q[bad], fn[bad] = P[bad], f[bad]
            P = Q
        Phi += (0.5 * dt * OMEGA_REF) * (f + fn)
        f = fn
        if j == sch["pulse_step"]:
            Phi = -Phi
        k = echo_of_step.get(j)
        if k is not None:
            out[k] = Phi
        if record_steps:
            trace.append(Phi.copy())
    res = dict(Phi=out, iv=own >= 0, vessel=own, R_um=float(R_um), cbv=geom.cbv,
               n_stuck=n_stuck, runtime=time.time() - t0, TE=np.asarray(sch["TE"]))
    if record_steps:
        res["Phi_steps"] = np.asarray(trace)
        res["t_steps"] = sch["t"]
    return res


def combine_histories(hs):
    """Pool spins of several runs (e.g. geometry seeds) into one history."""
    return dict(Phi=np.concatenate([h["Phi"] for h in hs], axis=1),
                iv=np.concatenate([h["iv"] for h in hs]),
                R_um=hs[0]["R_um"], cbv=float(np.mean([h["cbv"] for h in hs])), TE=hs[0]["TE"])


def signals_from_history(H, SO2, Hct=HCT, B0=B0_T, T2_ev=None, T2_iv=None):
    """EV, IV and total GESFIDE signals for one or many SO2 values.

    T2_ev / T2_iv: None (off), a number in ms, or (T2_iv only) a callable T2_iv(SO2).
    Returns dict of arrays shaped (nSO2, nTE): total, ev, iv (magnitudes), and the
    complex compartment means M_ev, M_iv (T2-free), plus f_iv (IV spin fraction).
    total = (1 - f_iv) * w_ev * M_ev + f_iv * w_iv * M_iv  (spin-density weighted).
    """
    SO2 = np.atleast_1d(np.asarray(SO2, float))
    TE = np.asarray(H["TE"], float)
    Phi, iv = H["Phi"], np.asarray(H["iv"], bool)
    n_iv, n = int(iv.sum()), len(iv)
    f_iv = n_iv / n
    s = phase_scale(SO2, Hct, B0)
    M_ev = np.zeros((len(SO2), len(TE)), complex)
    M_iv = np.zeros_like(M_ev)
    Pev, Piv = Phi[:, ~iv], Phi[:, iv]
    for i, si in enumerate(s):
        M_ev[i] = np.exp(1j * si * Pev).mean(1)
        if n_iv:
            M_iv[i] = np.exp(1j * si * Piv).mean(1)
    w_ev = np.exp(-TE / T2_ev)[None, :] if T2_ev else np.ones((1, len(TE)))
    if T2_iv is None:
        w_iv = np.ones((1, len(TE)))
    elif callable(T2_iv):
        w_iv = np.exp(-TE[None, :] / np.array([T2_iv(x) for x in SO2])[:, None])
    else:
        w_iv = np.exp(-TE / T2_iv)[None, :]
    tot = (1 - f_iv) * w_ev * M_ev + f_iv * w_iv * M_iv
    return dict(total=np.abs(tot), ev=np.abs(w_ev * M_ev), iv=np.abs(w_iv * M_iv),
                M_ev=M_ev, M_iv=M_iv, M_total=tot, f_iv=f_iv, SO2=SO2, TE=TE)


def T2_blood_lu2012(SO2, Hct=HCT):
    """Optional SO2-dependent blood T2 (ms), Lu et al. 2012 3T bovine-blood calibration."""
    a = -13.5 + 80.2 * Hct - 75.9 * Hct ** 2
    b = (-0.5 * Hct + 3.4 * Hct ** 2) * (1 - SO2)
    c = 247.4 * Hct * (1 - SO2) ** 2
    return 1000.0 / (a + b + c)


def simulate_gesfide(R_um, CBV, SO2, n_spins=20000, seeds=(0, 1), T2_ev=None, T2_iv=None,
                     static=False, D=D_UM2_MS, dt_max=DT_MAX_MS):
    """Convenience: build voxels for each seed, run, pool, return signals_from_history()."""
    hs = [phase_history(make_voxel(CBV, s), R_um, n_spins, seed=s, static=static, D=D, dt_max=dt_max)
          for s in seeds]
    return signals_from_history(combine_histories(hs), SO2, T2_ev=T2_ev, T2_iv=T2_iv)


# ------------------------------------------------------------------ fitting
def fit_rate(TE, S, mask):
    """Mono-exponential rate (1/s) of S over TE[mask] (ms)."""
    return -np.polyfit(np.asarray(TE)[mask] / 1000.0, np.log(np.asarray(S)[mask]), 1)[0]


def gesfide_rates(S, TE=TE_MS):
    """Ni 2015 readout: part A decays at R2 + R2', part B (rephasing) at R2 - R2'."""
    a = fit_rate(TE, S, TE < T180_MS)
    b = fit_rate(TE, S, (TE > T180_MS) & (TE <= 2 * T180_MS + 1e-9))
    return (a + b) / 2, (a - b) / 2

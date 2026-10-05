"""A3 spin-dephasing Monte Carlo (single infinite cylinder perpendicular to B0).

2D cross-section, vessel centred at the origin, B0 along +y (in-plane, theta = 90 deg).
Extravascular field offset (SI):  dBz = B0 * dchi/2 * (R/r)^2 * cos(2*phi),
phi measured from the B0 direction (+y)  ->  cos(2 phi) = (y^2 - x^2) / r^2.
Angular frequency offset  dw = gamma * dBz.

Spins: start uniformly outside the vessel in a square [-L, L)^2, L = fov_factor * R.
Boundary: periodic (the single-cylinder field is even in x and in y, so it is
continuous across the wrap).  Vessel wall: impermeable; a step that would end inside
the vessel is re-proposed with a fresh random step (up to 100 tries, else the spin
stays put for that step).
Diffusion: Gaussian steps, sigma = sqrt(2 D dt) per axis (true physical step length).
Phase: phi_GE(t) = int dw dt (trapezoid).  Spin echo with 180 deg at t180:
phi_SE(t) = phi_GE(t) - 2 phi_GE(t180) for t > t180 (sign flip of the phase at t180).
"""
import numpy as np

GAMMA = 2 * np.pi * 42.58e6          # rad/s/T
B0 = 3.0                             # T


def dchi_si(so2, hct=0.4, dchi0=0.27e-6):
    """SI susceptibility difference blood - tissue: 4*pi*dchi0_cgs*Hct*(1-SO2)."""
    return 4 * np.pi * dchi0 * hct * (1 - so2)


def dw_wall(so2, hct=0.4):
    """Angular frequency offset at the wall (r = R, phi = 0), rad/ms."""
    return GAMMA * B0 * dchi_si(so2, hct) / 2 / 1e3


def field(x, y, R, dw0):
    r2 = x * x + y * y
    return dw0 * R * R * (y * y - x * x) / (r2 * r2)


def simulate(R, so2=0.6, D=1.0, dt=0.05, T=100.0, t180=50.0, n=40000,
             fov_factor=5.0, n_keep=300, seed=1, static=False):
    rng = np.random.default_rng(seed)
    L = fov_factor * R
    dw0 = dw_wall(so2)
    # uniform start outside the vessel
    pos = np.empty((0, 2))
    while len(pos) < n:
        p = rng.uniform(-L, L, size=(2 * n, 2))
        p = p[(p ** 2).sum(1) > R * R]
        pos = np.vstack([pos, p])
    pos = pos[:n].copy()
    nsteps = int(round(T / dt))
    i180 = int(round(t180 / dt))
    t = np.arange(nsteps + 1) * dt
    sig = np.sqrt(2 * D * dt)

    phase = np.zeros(n)
    m_ge = np.empty(nsteps + 1, dtype=complex)
    m_se = np.empty(nsteps + 1, dtype=complex)
    keep_pos = np.empty((nsteps + 1, n_keep, 2))
    keep_phi = np.empty((nsteps + 1, n_keep))     # GE (unflipped) phase
    m_ge[0] = m_se[0] = 1.0
    keep_pos[0] = pos[:n_keep]
    keep_phi[0] = 0.0
    w_old = field(pos[:, 0], pos[:, 1], R, dw0)
    phi180 = None
    n_stuck = 0
    for k in range(1, nsteps + 1):
        if not static:
            step = rng.normal(0, sig, size=(n, 2))
            new = pos + step
            bad = (new ** 2).sum(1) < R * R
            tries = 0
            while bad.any() and tries < 100:
                idx = np.nonzero(bad)[0]
                new[idx] = pos[idx] + rng.normal(0, sig, size=(len(idx), 2))
                bad[idx] = (new[idx] ** 2).sum(1) < R * R
                tries += 1
            if bad.any():
                n_stuck += bad.sum()
                new[bad] = pos[bad]
            # periodic wrap
            new = (new + L) % (2 * L) - L
            pos = new
        w_new = field(pos[:, 0], pos[:, 1], R, dw0)
        phase += 0.5 * (w_old + w_new) * dt
        w_old = w_new
        if k == i180:
            phi180 = phase.copy()
        m_ge[k] = np.exp(1j * phase).mean()
        if phi180 is None:
            m_se[k] = m_ge[k]
        else:
            m_se[k] = np.exp(1j * (phase - 2 * phi180)).mean()
        keep_pos[k] = pos[:n_keep]
        keep_phi[k] = phase[:n_keep]
    return dict(t=t, m_ge=m_ge, m_se=m_se, keep_pos=keep_pos, keep_phi=keep_phi,
                phi180_keep=phi180[:n_keep], R=R, L=L, dw0=dw0, so2=so2, D=D, dt=dt,
                n=n, i180=i180, n_stuck=n_stuck, fov_factor=fov_factor)


if __name__ == "__main__":
    import sys
    for so2 in [0.6, 0.4, 0.2]:
        for ff in [5.0, 4.0]:
            line = f"SO2={so2:.1f} dw0={dw_wall(so2):.3f} rad/ms  L={ff}R  zeta={np.pi/(2*ff)**2:.3f} |"
            for R in [3.0, 30.0]:
                for st in [False, True]:
                    s = simulate(R, so2=so2, fov_factor=ff, n=20000, static=st, n_keep=10)
                    i = lambda tt: int(round(tt / s['dt']))
                    line += (f" R={R:g}{'S' if st else 'D'}: GE10={abs(s['m_ge'][i(10)]):.2f}"
                             f" GE50={abs(s['m_ge'][i(50)]):.2f} SE100={abs(s['m_se'][-1]):.3f}"
                             f" GE100={abs(s['m_ge'][-1]):.2f};")
            print(line, flush=True)

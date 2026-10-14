import sys, time, numpy as np, numba, scipy, matplotlib
print('python', sys.version.split()[0], '| numpy', np.__version__, '| numba', numba.__version__,
      '| scipy', scipy.__version__, '| matplotlib', matplotlib.__version__)
from boldswimsuite import BOLDvessel, BOLDgeometry, BOLDspins, BOLDsequence
t = time.time()
c = BOLDvessel.InfiniteCylinder3D(diameter=0.01, theta=np.pi/2, phi=0.0, origin=np.zeros(3), dchi=4.32e-8)
got = c.dBz_EV(np.array([[0., 0., 0.01]]), 3.0)[0]
exp = 3.0 * (4*np.pi*4.32e-8) / 2 * 0.25
print(f'field check: {got:.6e} vs {exp:.6e} -> {"OK" if np.isclose(got, exp) else "FAIL"} (numba compile {time.time()-t:.1f} s)')
vox = BOLDgeometry.ContinuousVoxel3D.from_random(CBV=0.02, B0=3.0, ADC=1e-3, labels=['v'], weights={'v': 1},
      diameter_distributions={'v': [0.006]}, dchis={'v': 3e-8}, num_vessels=30, seed=0, progressbar=False)
sp = BOLDspins.Spins3D(num_spins=2000, geometry=vox, dt=0.2, seed=0)
seq = BOLDsequence.SpinSequence(sp, [0, 250], [np.pi/2, np.pi], ['y', 'x'])
t = time.time(); tot, ev, iv = seq.walk(dt=0.2, num_steps=500, progressbar=False)
print(f'mini GESFIDE: {len(vox.vessels)} vessels, S(50 ms) = {tot[249]:.4f}, S(100 ms) = {tot[499]:.4f}, walk {time.time()-t:.1f} s')
print('INSTALL OK')

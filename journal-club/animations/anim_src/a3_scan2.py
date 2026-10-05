import numpy as np
from a3_spins_sim import simulate, dw_wall
for so2 in [0.6, 0.5, 0.4]:
    for R in [3.0, 30.0]:
        s = simulate(R, so2=so2, fov_factor=5.0, n=40000, n_keep=10, seed=2)
        i = lambda tt: int(round(tt / s['dt']))
        ge = ' '.join(f"{abs(s['m_ge'][i(tt)]):.3f}" for tt in [5,10,20,30,50,75,100])
        se = ' '.join(f"{abs(s['m_se'][i(tt)]):.3f}" for tt in [60,75,90,100])
        print(f"SO2={so2} R={R:g}: GE(5,10,20,30,50,75,100)= {ge} | SE(60,75,90,100)= {se} stuck={s['n_stuck']}", flush=True)

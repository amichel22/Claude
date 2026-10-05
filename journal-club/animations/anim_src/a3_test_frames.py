import sys, os
import a3_anim_spins as A
sims = A.get_sims()
fig, H = A.build_figure(sims)
sched = A.frame_schedule()
print("n frames", len(sched), "secs", len(sched)/A.FPS + 2)
idx = [int(a) for a in sys.argv[1:]]
for i in idx:
    k, fl, tag = sched[i]
    A.update(fig, H, sims, k, fl, tag)
    fig.savefig(os.path.join(A.HERE, f"_a3_test_{i:04d}.png"), dpi=A.DPI)
    print(i, k*A.DT, fl, tag)

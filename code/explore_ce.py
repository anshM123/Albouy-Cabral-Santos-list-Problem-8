"""Step 1: validate the 1-D reduction against brute force; step 2: random search for the
maximal number of critical points in the circle-ellipse case."""
import sys
import time
import numpy as np
from ce_core import unit_roots, count_branch, nmax_over_s, brute_force_2d

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 1)

print("== validation: polynomial roots vs branch sign changes vs 2-D Newton ==")
for t in range(8):
    e = rng.uniform(0.05, 0.95); sig = rng.uniform(0.02, 1.0); w = rng.uniform(0, np.pi / 2)
    s = np.exp(rng.uniform(np.log(0.2), np.log(5)))
    ang, gap, _ = unit_roots(e, sig, w, s)
    nm, np_ = count_branch(e, sig, w, s)
    bf = brute_force_2d(e, sig, w, s, grid=40)
    print(f"e={e:.3f} sig={sig:.3f} w={w:.3f} s={s:.3f}: poly={len(ang)}  branches={nm}+{np_}={nm+np_}  2D={len(bf)}")

print("== random search over (e, sig, w), max over s ==")
t0 = time.time()
hist = {}
best = []
for t in range(4000):
    e = rng.uniform(0.0, 0.995); sig = rng.uniform(0.0, 1.0); w = rng.uniform(0, np.pi / 2)
    n, lv = nmax_over_s(e, sig, w, m=6000)
    hist[n] = hist.get(n, 0) + 1
    if n >= 10:
        best.append((n, e, sig, w, lv))
print("histogram of max_s N:", dict(sorted(hist.items())), f"({time.time()-t0:.0f}s)")
print("configurations with N>=10:", len(best))
for b in best[:10]:
    print("  N=%d e=%.4f sig=%.4f w=%.4f s=%.5f" % b)

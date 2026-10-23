"""Step 3: how close can the circle-ellipse case get to 12 critical points?
For configurations with 10 real critical points, the two remaining roots of z^6 Pol form a
pair (z0, 1/conj z0) off the unit circle; delta = |log|z0|| > 0.  Minimise delta while
keeping 10 roots on the circle.  delta -> 0 in the interior would signal that 12 is reachable."""
import sys
import numpy as np
from scipy.optimize import minimize
from ce_core import unit_roots, nmax_over_s

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 7)


def delta_of(x):
    e, sig, w, ls = x
    if not (0.0 < e < 0.9999 and 0.0 <= sig <= 1.0):
        return 10.0
    ang, gap, roots = unit_roots(e, sig, w, np.exp(ls), tol=1e-7)
    if len(ang) < 10:
        return 5.0 + (10 - len(ang))
    return gap


starts = []
while len(starts) < 60:
    e = rng.uniform(0.6, 0.995); sig = rng.uniform(0.5, 1.0); w = rng.uniform(0, 0.9)
    n, lv = nmax_over_s(e, sig, w, m=6000)
    if n >= 10:
        starts.append((e, sig, w, np.log(lv)))

res = []
for x0 in starts:
    r = minimize(delta_of, np.array(x0), method="Nelder-Mead",
                 options=dict(xatol=1e-10, fatol=1e-12, maxiter=4000))
    res.append((r.fun, *r.x))
res.sort()
print("smallest delta reached (delta, e, sig, w, s):")
for d, e, sig, w, ls in res[:8]:
    print(f"  delta={d:.3e}  e={e:.6f} sig={sig:.6f} w={w:.6f} s={np.exp(ls):.6f}")
print("largest of the local minima:", f"{res[-1][0]:.3e}")

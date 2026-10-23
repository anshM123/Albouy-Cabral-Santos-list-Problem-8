"""AUDIT a8: independent search for configurations where Pol is (nearly) real-rooted (12 roots on the circle).

Objective on x = (eta, l, w, t) (a-chart, e = sin eta, S = t/(1-t)):
  gap(x) = min |log|z|| over the roots of Q(z) = z^6 Pol that are NOT on the circle, provided >= 10 roots are on it.
The working notes claimed min gap >= 0.836 (60 Nelder-Mead runs, e in [0.6, 0.995]).  Here: 400 starts spread over the whole
domain (including near the faces l = 1, w = 0, w = pi/2, t -> 0, t -> 1), Nelder-Mead in an unconstrained
reparametrisation (sigmoid), plus a report restricted to the bulk domain e in [0.2, 0.95].
usage: python a8_gap_search.py NSTART SEED EMIN EMAX
"""
import sys
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, '.')
from nbport import phipsi_coeffs

NS = int(sys.argv[1]) if len(sys.argv) > 1 else 200
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 8)
EMIN = float(sys.argv[3]) if len(sys.argv) > 3 else 0.2
EMAX = float(sys.argv[4]) if len(sys.argv) > 4 else 0.95
lo = np.array([np.arcsin(EMIN), 0.0, 0.0, 0.0]); hi = np.array([np.arcsin(EMAX), 1.0, np.pi / 2, 1.0])


def roots(x):
    (aP, sP), (aS, sS) = phipsi_coeffs(0, x[0], x[1], x[2])
    t = x[3]
    a = (1 - t) * aP - t * aS; s = (1 - t) * sP - t * sS
    c = np.zeros(13, complex); c[6] = a[0]
    for k in range(1, 7):
        c[6 + k] = (a[k] - 1j * s[k]) / 2; c[6 - k] = (a[k] + 1j * s[k]) / 2
    mx = np.abs(c).max()
    top = 12
    while top > 0 and abs(c[top]) < 1e-14 * mx:
        top -= 1
    low = 0
    while low < top and abs(c[low]) < 1e-14 * mx:
        low += 1
    return np.roots(c[low:top + 1][::-1]), (12 - top) + low


def gap(x):
    z, lost = roots(x)
    lz = np.sort(np.abs(np.log(np.abs(z))))
    n_on = int(np.sum(lz < 1e-7))
    if n_on >= 10 and lost == 0:
        return float(lz[n_on]) if n_on < len(lz) else 0.0
    return 5.0 + (10 - min(n_on, 10))


def to_x(y):
    return lo + (hi - lo) / (1 + np.exp(-y))


def to_y(x):
    f = np.clip((x - lo) / (hi - lo), 1e-9, 1 - 1e-9)
    return np.log(f / (1 - f))


starts = []
tries = 0
while len(starts) < NS and tries < 200000:
    tries += 1
    x = lo + (hi - lo) * rng.uniform(size=4)
    # bias some starts towards faces
    k = rng.integers(0, 6)
    if k == 0: x[1] = 1 - 10 ** rng.uniform(-6, -1)
    if k == 1: x[2] = 10 ** rng.uniform(-6, -1)
    if k == 2: x[2] = np.pi / 2 - 10 ** rng.uniform(-6, -1)
    if k == 3: x[3] = 1 - 10 ** rng.uniform(-6, -1)
    if k == 4: x[3] = 10 ** rng.uniform(-6, -1)
    if gap(x) < 5:
        starts.append(x)
print(f"{len(starts)} starts with >= 10 real roots (from {tries} random draws)")
res = []
for x0 in starts:
    r = minimize(lambda y: gap(to_x(y)), to_y(x0), method='Nelder-Mead',
                 options=dict(xatol=1e-9, fatol=1e-12, maxiter=3000))
    res.append((r.fun, to_x(r.x)))
res.sort(key=lambda a: a[0])
print(f"domain e in [{EMIN}, {EMAX}]: smallest gaps found (gap, e, l, w, t):")
for g, x in res[:10]:
    print("   gap=%.4f  e=%.5f l=%.6f w=%.6f t=%.6f" % (g, np.sin(x[0]), x[1], x[2], x[3]))
print("median of local minima: %.4f" % np.median([g for g, _ in res]))

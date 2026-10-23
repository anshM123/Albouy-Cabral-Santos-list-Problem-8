"""AUDIT a7b: root structure of Pol in the second failure cluster of costmap (e~0.4358, l~0.99, w~0.261, t~0.197)."""
import json, sys
import numpy as np
from a3_reduction import crit_2d
sys.path.insert(0, '.')
from nbport import phipsi_coeffs

def roots_all(chart, x0, l, w, t):
    (aP, sP), (aS, sS) = phipsi_coeffs(chart, x0, l, w)
    a = (1 - t) * aP - t * aS; s = (1 - t) * sP - t * sS
    c = np.zeros(13, complex); c[6] = a[0]
    for k in range(1, 7):
        c[6 + k] = (a[k] - 1j * s[k]) / 2; c[6 - k] = (a[k] + 1j * s[k]) / 2
    z = np.roots(c[::-1])
    return z

d = json.load(open('../code/costmap.json'))
F = [f for r in d['results'] for f in r['fails']]
B = np.array([f[0] for f in F])
mid = np.stack([(B[:, 2*i] + B[:, 2*i+1]) / 2 for i in range(4)], 1)
sel = mid[:, 1] > 0.98
print("cluster-2 leaves:", sel.sum())
hull = [(B[sel, 2*i].min(), B[sel, 2*i+1].max()) for i in range(4)]
print("hull: eta", hull[0], "e", np.sin(hull[0]), "l", hull[1], "w", hull[2], "t", hull[3])
rng = np.random.default_rng(1)
mx = 0; hist = {}
best = None
for it in range(20000):
    x = [rng.uniform(*hull[i]) for i in range(4)]
    z = roots_all(0, *x)
    lz = np.abs(np.log(np.abs(z)))
    n = int(np.sum(lz < 1e-6)); gap = np.sort(lz[lz >= 1e-6])
    hist[n] = hist.get(n, 0) + 1
    if n > mx or (n == mx and best is not None and gap.size and gap[0] < best[1]):
        mx = n; best = (x, gap[0] if gap.size else np.inf)
print("real-root count histogram over 20000 samples of the hull:", dict(sorted(hist.items())))
x = best[0]
print("example point:", np.round(x, 6), "e=%.6f" % np.sin(x[0]), "min off-circle |log|z|| =", best[1])
z = roots_all(0, *x)
order = np.argsort(np.angle(z))
for zz in z[order]:
    print("   root angle %+.6f  |log|z|| %.3e" % (np.angle(zz), abs(np.log(abs(zz)))))
e = np.sin(x[0]); l = x[1]; w = x[2]; S = x[3] / (1 - x[3])
I = np.arcsin(l)
cps = crit_2d(e, I, w, np.sqrt(S), n=72)
print("2-D Newton critical points at that point:", len(cps))

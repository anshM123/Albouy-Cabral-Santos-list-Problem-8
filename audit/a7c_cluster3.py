"""AUDIT a7c: root structure at the bulk_a3 failure cluster (e = 0.2, l ~ 0.2, w ~ 0.193, t ~ 0.509)."""
import json, sys
import numpy as np
sys.path.insert(0, '.')
from nbport import phipsi_coeffs
from a3_reduction import crit_2d
recs = [json.loads(x) for x in open('../code/bulk_a3.json.fails.jsonl') if x.strip()]
F = [f for r in recs for f in r['fails']]
B = np.array([f[0] for f in F])
hull = [(B[:, 2*i].min(), B[:, 2*i+1].max()) for i in range(4)]
print("hull eta", hull[0], "-> e", np.sin(hull[0]), " l", hull[1], " w", hull[2], " t", hull[3])
rng = np.random.default_rng(3)
def roots(x):
    (aP, sP), (aS, sS) = phipsi_coeffs(0, x[0], x[1], x[2]); t = x[3]
    a = (1-t)*aP - t*aS; s = (1-t)*sP - t*sS
    c = np.zeros(13, complex); c[6] = a[0]
    for k in range(1, 7): c[6+k] = (a[k]-1j*s[k])/2; c[6-k] = (a[k]+1j*s[k])/2
    return np.roots(c[::-1])
hist = {}; best = None
for it in range(20000):
    x = [rng.uniform(*hull[i]) for i in range(4)]
    z = roots(x); lz = np.abs(np.log(np.abs(z)))
    n = int(np.sum(lz < 1e-6)); hist[n] = hist.get(n, 0) + 1
    on = np.sort(np.angle(z[lz < 1e-6]))
    sep = np.min(np.diff(np.concatenate([on, on[:1] + 2*np.pi]))) if len(on) > 1 else np.inf
    if best is None or sep < best[1]:
        best = (x, sep, n, np.sort(lz[lz >= 1e-6])[:2])
print("count histogram:", dict(sorted(hist.items())))
x, sep, n, offs = best
print("tightest point: x=", np.round(x, 7), "e=%.7f" % np.sin(x[0]), "#real", n, "min root separation %.2e" % sep, "off-circle |log|z||", offs)
z = roots(x); lz = np.abs(np.log(np.abs(z)))
for zz in sorted(z, key=np.angle):
    print("   angle %+.6f |log|z|| %.2e" % (np.angle(zz), abs(np.log(abs(zz)))))
S = x[3]/(1-x[3]); print("2-D Newton count:", len(crit_2d(np.sin(x[0]), np.arcsin(x[1]), x[2], np.sqrt(S), n=80)))

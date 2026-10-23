"""Empirical root geometry in the C3 region {e in [0.99,1], l in [0.95,1], w in [1.4137, pi/2], mu in [0, 0.0505]}:
real roots of Q(th) (th = tan((f - pi)/2)) split by |th|; weighted radius delta = max(Qz, sqrt(kappa), sqrt(eps), mu^(1/3))."""
import numpy as np, math, sys
import sympy as sp
th, Pz, Qz, e, mu = sp.symbols('th Pz Qz e mu', real=True)
cf, sf = sp.symbols('cf sf', real=True)
zeta = Pz * cf + Qz * sf; dzeta = -Pz * sf + Qz * cf
B = -2 * zeta * dzeta + 2 * e * (sf - Qz * zeta); D = 1 + e * cf; k = 2 * e * (1 + e)
P = sp.expand(mu**2 * k**2 * sf**2 * (1 - zeta**2) - D**2 * B**2)
den = 1 + th**2
Q = sp.Poly(sp.expand(sp.cancel(P.subs({cf: -(1 - th**2) / den, sf: -2 * th / den}) * den**6)), th)
cfs = [sp.lambdify((Pz, Qz, e, mu), c, 'math') for c in Q.all_coeffs()]
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
rows = []
for it in range(20000):
    ee = 1 - 0.01 * rng.random() ** 2; l = 1 - 0.05 * rng.random() ** 2; w = math.pi / 2 - (math.pi / 2 - 1.4137) * rng.random() ** 2
    m = 0.0505 * rng.random() ** 3
    pz = l * math.sin(w); qz = l * math.cos(w)
    c = np.array([f(pz, qz, ee, m) for f in cfs], dtype=float)
    r = np.roots(c)
    re = np.sort(np.abs(r[np.abs(r.imag) < 1e-7 * (1 + np.abs(r))].real))
    nreal = len(re) + (12 - (len(c) - 1 - np.argmax(np.abs(c) > 0)) if False else 0)
    delta = max(qz, math.sqrt(max(1 - pz, 0)), math.sqrt(1 - ee), m ** (1 / 3))
    rows.append((delta, re, len(r)))
# for each sample: largest |th| among roots with |th| < 1 (near-pi group) relative to delta; smallest |th| among |th| >= 1
near_ratio = []; far_min = []; cnt = {}
for delta, re, nr in rows:
    nearr = re[re < 1.0]; farr = re[re >= 1.0]
    if len(nearr):
        near_ratio.append(nearr.max() / delta)
    if len(farr):
        far_min.append(farr.min())
    cnt[len(re)] = cnt.get(len(re), 0) + 1
print('real-root counts (finite th):', dict(sorted(cnt.items())))
print('near group max |th|/delta: quantiles', np.quantile(near_ratio, [0.5, 0.9, 0.99, 1.0]))
print('far group min |th|:', np.quantile(far_min, [0, 0.01, 0.1, 0.5]) if far_min else None)
print('delta range', min(r[0] for r in rows), max(r[0] for r in rows))
# split counts: near |th| <= 3 delta, mid 3 delta < |th| < 1.11, far |th| >= 1.11 (plus roots at infinity = 12 - deg)
nn = {}; ff = {}; mid = 0; joint = {}
for delta, re, nr in rows:
    a = int(np.sum(re <= 3 * delta)); b = int(np.sum((re > 3 * delta) & (re < 1.11))); c = int(np.sum(re >= 1.11)) + (12 - nr)
    nn[a] = nn.get(a, 0) + 1; ff[c] = ff.get(c, 0) + 1; mid += b; joint[(a, c)] = joint.get((a, c), 0) + 1
print('near counts', dict(sorted(nn.items())), 'far(+inf) counts', dict(sorted(ff.items())), 'mid-region roots', mid)
print('joint (near, far):', dict(sorted(joint.items())))

"""Numerical test of Lemmas N1, N2 and A (CORNER_E1.md) on random corner parameters."""
import sys
import numpy as np

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 4)


def F(f, e, l, w, mu, sign):
    p, q = l * np.sin(w), l * np.cos(w)
    zeta = p * np.cos(f) + q * np.sin(f)
    zp = -p * np.sin(f) + q * np.cos(f)
    sig = np.sqrt(np.maximum(1 - zeta ** 2, 0))
    D = 1 + e * np.cos(f)
    B = -2 * zeta * zp + 2 * e * (np.sin(f) - q * zeta)
    return 2 * mu * e * (1 + e) * np.sin(f) * sig + sign * D * B      # sign=-1: F_-, +1: F_+


def nzeros(v):
    return int(np.sum(np.sign(v[:-1]) != np.sign(v[1:])))


phi0 = 0.3
worst = dict(N1=np.inf, N2=-np.inf)
viol = dict(N1=0, N2=0, A_minus=0, A_plus=0, total=0)
for it in range(int(sys.argv[2]) if len(sys.argv) > 2 else 3000):
    eps = np.exp(rng.uniform(np.log(1e-8), np.log(1e-2)))
    mu = np.exp(rng.uniform(np.log(1e-8), np.log(1e-2)))
    l = rng.uniform(0, 1); w = rng.uniform(0, np.pi / 2)
    if 1 - (l * np.sin(w)) ** 2 < 0.05:          # c >= c0 (C3 excluded)
        continue
    e = 1 - eps
    # aphelion arc
    f = np.pi + np.linspace(-phi0, phi0, 200001)
    h = f[1] - f[0]
    Fm = F(f, e, l, w, mu, -1); Fp = F(f, e, l, w, mu, +1)
    d3 = np.diff(Fm, 3) / h ** 3
    d1 = np.diff(Fp) / h
    worst['N1'] = min(worst['N1'], d3.min()); worst['N2'] = max(worst['N2'], d1.max())
    viol['N1'] += d3.min() <= 0
    viol['N2'] += d1.max() >= 0
    # away from aphelion
    g = np.linspace(-np.pi + phi0, np.pi - phi0, 400001)
    zm = nzeros(F(g, e, l, w, mu, -1)); zp_ = nzeros(F(g, e, l, w, mu, +1))
    viol['A_minus'] += zm > 3; viol['A_plus'] += zp_ > 3
    viol['total'] += (zm + nzeros(Fm)) + (zp_ + nzeros(Fp)) > 10
print("violations:", viol)
print("min F_-''' on aphelion arc: %.3g ; max F_+' on aphelion arc: %.3g" % (worst['N1'], worst['N2']))

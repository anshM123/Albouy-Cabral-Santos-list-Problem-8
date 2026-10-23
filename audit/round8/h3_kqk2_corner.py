"""AUDIT-8 h3: the KQK2 corner where level r2 left 2,146 failed leaves (all with certified total 12):
e in [0.99998, 1], l in [0.99922, 1], w in [0.156, 0.314], t in [0.99999992, 1] (mu <= 2.8e-4).
Counts N_mult of P = mu^2 k^2 sin^2 f sigma^2 - D^2 B^2 (100 digits, own formula) at the limit e = 1, mu = 0 (non-physical)
and at nearby physical points, to see whether a uniform bound <= 11 can exist on boxes touching e = 1, t = 1."""
import mpmath as mp
mp.mp.dps = 100


def nmult(e, l, w, mu):
    e, l, w, mu = [mp.mpf(v) for v in (e, l, w, mu)]
    pz, qz = l * mp.sin(w), l * mp.cos(w)
    N = 32; vals = []
    for n in range(N):
        f = 2 * mp.pi * n / N
        zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
        B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta); D = 1 + e * mp.cos(f); k = 2 * e * (1 + e)
        vals.append(mu ** 2 * k ** 2 * mp.sin(f) ** 2 * (1 - zeta ** 2) - D ** 2 * B ** 2)
    c = [mp.fsum(vals[n] * mp.expj(-kk * 2 * mp.pi * n / N) for n in range(N)) / N for kk in range(-6, 7)]
    sc = max(abs(v) for v in c); lo = 0
    while abs(c[lo]) < mp.mpf(10) ** -80 * sc:
        lo += 1
    hi = 12
    while abs(c[hi]) < mp.mpf(10) ** -80 * sc:
        hi -= 1
    z = mp.polyroots(c[lo:hi + 1][::-1], maxsteps=1000, extraprec=1500)
    on = sorted([mp.arg(r) % (2 * mp.pi) for r in z if abs(mp.log(abs(r))) < mp.mpf(10) ** -20])
    return len(on), [mp.nstr(x, 6) for x in on]


for w in (0.16, 0.2, 0.25, 0.3, 0.05, 0.4):
    n, roots = nmult(1, 1, w, 0)
    print(f'limit e = 1, mu = 0, l = 1, w = {w}: N_mult = {n}   roots f = {roots}')
for (e, l, w, mu) in ((1 - 1e-5, 0.9995, 0.2, 1e-4), (1 - 2e-5, 1.0, 0.25, 2e-4), (1 - 1e-6, 0.9999, 0.3, 1e-5), (1 - 1e-5, 1.0, 0.16, 5e-5)):
    n, roots = nmult(e, l, w, mu)
    print(f'physical e = {e}, l = {l}, w = {w}, mu = {mu}: N_mult = {n}')

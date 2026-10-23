"""AUDIT-3 c8: KA leaves (away arc |f - pi| >= PHI_A for the C3 cell, bound 2): 100-digit count of the real roots of P
(R-normalised, mu = q/R) on the away arc at sampled points of every leaf."""
import json, math, random
import mpmath as mp
mp.mp.dps = 100
random.seed(8)
D = json.load(open('../../code/runs/KA.json')); H = D['header']; PHI_A = mp.mpf(H['phi_a'])
def P(e, pz, qz, mu, f):
    zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
    B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta); D_ = 1 + e * mp.cos(f); k = 2 * e * (1 + e)
    return mu ** 2 * k ** 2 * mp.sin(f) ** 2 * (1 - zeta ** 2) - D_ ** 2 * B ** 2
def roots_away(e, pz, qz, mu):
    N = 24
    vals = [P(e, pz, qz, mu, 2 * mp.pi * n / N) for n in range(N)]
    c = [mp.fsum(vals[n] * mp.exp(-1j * k * 2 * mp.pi * n / N) for n in range(N)) / N for k in range(-6, 7)]
    sc = max(abs(x) for x in c)
    while abs(c[-1]) < mp.mpf(10) ** -80 * sc: c.pop()
    lo = 0
    while abs(c[lo]) < mp.mpf(10) ** -80 * sc: lo += 1
    z = mp.polyroots(c[lo:][::-1], maxsteps=400, extraprec=300)
    on = [r for r in z if abs(mp.log(abs(r))) < mp.mpf(10) ** -20]
    return sum(1 for r in on if abs(mp.pi - abs(mp.arg(r))) >= PHI_A - mp.mpf(10) ** -12 or abs(mp.arg(r)) <= mp.pi - PHI_A), len(on)
viol = 0; n = 0; mx = 0
for r in D['results']:
    bx = r['box']
    assert r['tree'] in 'abcdefghijk' and len(r['tree']) == 1
    bound = ord(r['tree']) - ord('a')
    for k in range(3):
        e = mp.mpf(random.uniform(bx[0], min(bx[1], 1 - 1e-12))); l = mp.mpf(random.uniform(bx[2], bx[3]))
        w = mp.mpf(random.uniform(bx[4], bx[5])); mu = mp.mpf(random.uniform(max(bx[6], 1e-7), bx[7]))
        a, tot = roots_away(e, l * mp.sin(w), l * mp.cos(w), mu)
        n += 1; mx = max(mx, a)
        if a > bound: viol += 1; print('VIOLATION', bx, a, bound)
print(f"KA: {len(D['results'])} leaves x 3 points = {n}; violations {viol}; max real roots on the away arc {mx} (bound 2)")

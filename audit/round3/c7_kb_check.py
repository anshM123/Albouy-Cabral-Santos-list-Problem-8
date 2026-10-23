"""AUDIT-3 c7: empirical soundness of the KB band engine (kb_run.certify) and of corner2.check_plus_sub.
Random boxes in the band e in [0.99, 1], l in [0.95, 1], w in [3pi/20, 9pi/20], mu in [0, 0.0505]; when certified with
total T, sampled points (e < 1, mu > 0) must have #real roots of P on the circle (100 digits, companion polynomial in
z = e^{if}, multiplicity counted via near-circle clusters) <= T.  When check_plus_sub accepts a box, the number of
sign changes of F_+ on |phi| <= 0.1 at sampled points must be <= 1.   (.venv python; run from audit/round3)"""
import os, sys, math, random
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.join(HERE, '..', '..', 'code')
sys.path.insert(0, CODE); sys.path.insert(0, os.path.join(CODE, 'c3'))
os.chdir(CODE)                                   # kb_run inserts 'c3' relative to cwd (read only)
import kb_run
kb_run.init()
import corner2
mp.mp.dps = 100
random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 9)
NB = int(sys.argv[1]) if len(sys.argv) > 1 else 30


def P(e, pz, qz, mu, f):
    zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
    B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta); D = 1 + e * mp.cos(f); k = 2 * e * (1 + e)
    return mu ** 2 * k ** 2 * mp.sin(f) ** 2 * (1 - zeta ** 2) - D ** 2 * B ** 2


def nreal(e, pz, qz, mu):
    N = 24
    vals = [P(e, pz, qz, mu, 2 * mp.pi * n / N) for n in range(N)]
    c = [mp.fsum(vals[n] * mp.exp(-1j * k * 2 * mp.pi * n / N) for n in range(N)) / N for k in range(-6, 7)]
    sc = max(abs(x) for x in c)
    while abs(c[-1]) < mp.mpf(10) ** -80 * sc:
        c.pop()
    lo = 0
    while abs(c[lo]) < mp.mpf(10) ** -80 * sc:
        lo += 1
    z = mp.polyroots(c[lo:][::-1], maxsteps=400, extraprec=300)
    return sum(1 for r in z if abs(mp.log(abs(r))) < mp.mpf(10) ** -20)


def Fplus_changes(e, pz, qz, mu):
    fs = [mp.pi - mp.mpf('0.1') + mp.mpf('0.2') * k / 4000 for k in range(4001)]
    vals = []
    for f in fs:
        zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
        sig = mp.sqrt(1 - zeta ** 2); B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta); D = 1 + e * mp.cos(f)
        vals.append(mu * 2 * e * (1 + e) * mp.sin(f) * sig + D * B)
    return sum(1 for a, b in zip(vals, vals[1:]) if mp.sign(a) != mp.sign(b))


st = dict(boxes=0, cert=0, viol=0, n2_acc=0, n2_viol=0); hist = {}
W0, W1 = 3 * math.pi / 20, 9 * math.pi / 20
for ib in range(NB):
    c = [random.uniform(0.99, 1.0), random.uniform(0.95, 1.0), random.uniform(W0, W1), random.uniform(0, 0.0505)]
    r = np.array([0.01, 0.05, W1 - W0, 0.0505]) * 0.02 * np.exp(np.array([random.uniform(-3, 0) for _ in range(4)]))
    lims = [(0.99, 1.0), (0.95, 1.0), (W0, W1), (0.0, 0.0505)]
    bx = []
    for i in range(4):
        bx += [max(c[i] - r[i], lims[i][0]), min(c[i] + r[i], lims[i][1])]
    ok, tot = kb_run.certify(tuple(bx))
    st['boxes'] += 1
    pz0 = kb_run.dn(bx[2] * kb_run.dn(math.sin(bx[4]))); pz1 = kb_run.up(bx[3] * kb_run.up(math.sin(bx[5])))
    qz0 = kb_run.dn(bx[2] * kb_run.dn(math.cos(bx[5]))); qz1 = kb_run.up(bx[3] * kb_run.up(math.cos(bx[4])))
    n2 = corner2.check_plus_sub(bx[0], bx[1], pz0, pz1, qz0, qz1, 0.1 + 1e-12, 32)
    for k in range(4):
        e = mp.mpf(random.uniform(bx[0], min(bx[1], 1 - 1e-12))); l = mp.mpf(random.uniform(bx[2], bx[3]))
        w = mp.mpf(random.uniform(bx[4], bx[5])); mu = mp.mpf(random.uniform(max(bx[6], 1e-7), bx[7]))
        pz, qz = l * mp.sin(w), l * mp.cos(w)
        if ok:
            n = nreal(e, pz, qz, mu)
            if n > tot:
                st['viol'] += 1; print('VIOLATION', bx, (e, l, w, mu), n, tot)
        if n2:
            ch_ = Fplus_changes(e, pz, qz, mu)
            if ch_ > 1:
                st['n2_viol'] += 1; print('N2 VIOLATION', bx, ch_)
    if ok:
        st['cert'] += 1; hist[int(tot)] = hist.get(int(tot), 0) + 1
    st['n2_acc'] += bool(n2)
print(st, 'certified totals', dict(sorted(hist.items())))

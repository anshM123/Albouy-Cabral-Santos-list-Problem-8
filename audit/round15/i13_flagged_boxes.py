"""AUDIT-15 i13 (as AUDIT-12 i11): count leaves of the new files that g1's strict criterion flagged (recorded bound b,
today's first-success total > b, <= 11).  For the flagged boxes with b <= 7 and every 6th of the others:
(a) 60-digit real-root counts (my own P) at the centre and 8 random points: never above b?
(b) today's engine, every option (perihelion count on/off, both split angles, regular and deep): smallest total."""
import os, sys, math, random, re
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.normpath(os.path.join(HERE, '..', '..', 'code'))
sys.path.insert(0, CODE); sys.path.insert(0, os.path.join(CODE, 'c3'))
src = os.path.join(HERE, 'out', 'G_g1.out')
os.chdir(CODE)
mp.mp.dps = 60
random.seed(7)
ALL = []
for line in open(src):
    m = re.match(r'\s+COUNT RECERT PROBLEM \(([^)]*)\) (\d+) True (\d+)', line)
    if m:
        ALL.append((int(m.group(2)), int(m.group(3)), tuple(float(v) for v in m.group(1).split(','))))
BOXES = [b for b in ALL if b[0] <= 7] + [b for i, b in enumerate(ALL) if b[0] > 7 and i % 6 == 0]
print(f'flagged by g1: {len(ALL)}; checked here: {len(BOXES)}', flush=True)


def P(e, pz, qz, mu, f):
    zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
    B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta); D = 1 + e * mp.cos(f); k = 2 * e * (1 + e)
    return mu ** 2 * k ** 2 * mp.sin(f) ** 2 * (1 - zeta ** 2) - D ** 2 * B ** 2


def nreal(e, l, w, t):
    pz, qz, mu = l * mp.sin(w), l * mp.cos(w), mp.sqrt((1 - t) / t)
    N = 32
    vals = [P(e, pz, qz, mu, 2 * mp.pi * n / N) for n in range(N)]
    c = [mp.fsum(vals[n] * mp.expj(-kk * 2 * mp.pi * n / N) for n in range(N)) / N for kk in range(-6, 7)]
    co = c[::-1]
    sc = max(abs(v) for v in co)
    while abs(co[0]) < mp.mpf(10) ** -50 * sc:
        co = co[1:]
    z = mp.polyroots(co, maxsteps=3000, extraprec=1500)
    return sum(1 for r in z if abs(mp.log(abs(r))) < mp.mpf(10) ** -25)


import kq_run
kq_run.init([1.0, 1.0, math.pi / 2, 1.0]); kq_run.DEEP_ALWAYS = True
under = 0; notrep = 0
for b, tnow, bx in BOXES:
    pts = [[0.5 * (bx[2 * i] + bx[2 * i + 1]) for i in range(4)]] + [[random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(4)] for _ in range(8)]
    ns = []
    for x in pts:
        x[0] = min(x[0], 1 - 1e-12); x[3] = min(x[3], 1 - 1e-15)
        ns.append(nreal(*[mp.mpf(v) for v in x]))
    best = 99
    for peri in (None, kq_run.PHI_P):
        for phi0 in sorted({kq_run.phi_split(bx), kq_run.PHI_0}):
            for rounds, cap in ((kq_run.ROUNDS1, kq_run.CAP1), (kq_run.ROUNDS, kq_run.CAP)):
                ok, tot = kq_run.count_total(bx, rounds, cap, phi0, peri)
                if ok:
                    best = min(best, int(tot))
    under += max(ns) > b; notrep += best > b
    print(f'   recorded {b}, g1 first-success {tnow}: 60-digit counts max {max(ns)} (<= recorded: {max(ns) <= b}); '
          f'smallest certified total over all options {best} (== recorded: {best == b})', flush=True)
print(f'checked {len(BOXES)}: 60-digit count above the recorded bound: {under}; recorded bound not reproduced by any '
      f'option of today\'s engine: {notrep}', flush=True)

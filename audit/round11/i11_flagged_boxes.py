"""AUDIT-12 i11: the KQ count leaves that g1's strict criterion flagged (recorded bound b, today's kq_run total > b, <= 11).
(a) 100-digit real-root counts (my own P, FFT + companion roots) at the centre and 8 random points of each box: the
    recorded bound must not be exceeded (an old-engine under-count would be a soundness red flag);
(b) today's engine with every option (split angle phi_split / PHI_0, perihelion count on / off, deep): the minimum."""
import os, sys, math, random
import mpmath as mp
CODE = r'<repo>\code'
sys.path.insert(0, CODE); sys.path.insert(0, os.path.join(CODE, 'c3')); os.chdir(CODE)
mp.mp.dps = 60
random.seed(5)
BOXES = [('KQ1', 9, (0.962425, 0.97495, 0.7000000000000001, 0.75, 0.7068583470577035, 0.7853981633974483, 0.898, 0.9477500000000001)),
         ('KQ1.r1', 8, (0.987475, 1.0, 0.9, 0.95, 0.9424777960769379, 0.9817477042468103, 0.9228750000000001, 0.9477500000000001)),
         ('KQ1.r1', 10, (0.987475, 1.0, 0.95, 0.975, 0.8639379797371931, 0.9424777960769379, 0.9228750000000001, 0.9477500000000001)),
         ('KQB.r3', 8, (0.99, 0.995, 0.95, 0.955, 0.6872233929727672, 0.6970408700152353, 0.9993043750000001, 0.9993242500000001)),
         ('KQK1.r2', 10, (0.9924999999999999, 0.9937499999999999, 0.7445312499999999, 0.7453124999999999, 0.019634954084936207, 0.022089323345553233, 0.999987578125, 0.9999900625)),
         ('KQK1.r4', 9, (0.9999781250001094, 0.9999812500000937, 0.9375, 0.95, 0.8344855486097889, 0.8443030256522569, 0.999999837574185, 0.9999998400000273))]


def P(e, pz, qz, mu, f):
    zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
    B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta); D = 1 + e * mp.cos(f); k = 2 * e * (1 + e)
    return mu ** 2 * k ** 2 * mp.sin(f) ** 2 * (1 - zeta ** 2) - D ** 2 * B ** 2


def nreal(e, l, w, t):
    pz, qz, mu = l * mp.sin(w), l * mp.cos(w), mp.sqrt((1 - t) / t)
    N = 32
    vals = [P(e, pz, qz, mu, 2 * mp.pi * n / N) for n in range(N)]
    c = [mp.fsum(vals[n] * mp.expj(-kk * 2 * mp.pi * n / N) for n in range(N)) / N for kk in range(-6, 7)]
    z = mp.polyroots(c[::-1], maxsteps=2000, extraprec=1500)
    return sum(1 for r in z if abs(mp.log(abs(r))) < mp.mpf(10) ** -25)


import kq_run
kq_run.init([1.0, 1.0, math.pi / 2, 1.0]); kq_run.DEEP_ALWAYS = True
for name, b, bx in BOXES:
    pts = [[0.5 * (bx[2 * i] + bx[2 * i + 1]) for i in range(4)]] + [[random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(4)] for _ in range(8)]
    ns = []
    for x in pts:
        x[0] = min(x[0], 1 - 1e-12); x[3] = min(x[3], 1 - 1e-15)
        ns.append(nreal(*[mp.mpf(v) for v in x]))
    tots = []
    for peri in (None, kq_run.PHI_P):
        for phi0 in sorted({kq_run.phi_split(bx), kq_run.PHI_0}):
            for rounds, cap in ((kq_run.ROUNDS1, kq_run.CAP1), (kq_run.ROUNDS, kq_run.CAP)):
                ok, tot = kq_run.count_total(bx, rounds, cap, phi0, peri)
                tots.append((peri, phi0, rounds, bool(ok), int(min(tot, 99))))
    best = min(t[4] for t in tots if t[3]) if any(t[3] for t in tots) else None
    print(f'{name}: recorded {b}; 100-digit real roots at 9 points {ns} (max {max(ns)} <= recorded: {max(ns) <= b}); '
          f'today\'s engine, smallest certified total over all options: {best}; all: {tots}', flush=True)

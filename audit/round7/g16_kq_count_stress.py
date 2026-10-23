"""AUDIT-7 g16: stress test of the KQ COUNT certificate (kq_run.certify: away count + KBQe near count; since 08:40 the
away count uses nb_cert2.away_count_t2) on random boxes of the five KQ domains (widths 1e-5..1e-1 of the KQ scale,
biased towards the faces).  Every certified box is checked at 100 digits (MY OWN P(f)) at its centre, 4 corners and 3
random points: N_mult(P) on the whole circle must be <= the certified total (<= 11).
usage: python g16_kq_count_stress.py NBOX SEED"""
import os, sys, math, random
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.normpath(os.path.join(HERE, '..', '..', 'code'))
NBOX, SEED = int(sys.argv[1]), int(sys.argv[2])
sys.path.insert(0, CODE); sys.path.insert(0, os.path.join(CODE, 'c3')); sys.path.insert(0, HERE)
os.chdir(CODE)
import kq_run
kq_run.init()
mp.mp.dps = 100
random.seed(SEED)
HP = math.pi / 2
DOMS = {'KQ1': (0.9499, 1.0, 0, 1, 0, HP, 0.5, 0.9975), 'KQ2': (0.9499, 0.99, 0, 1, 0, HP, 0.9975, 1.0),
        'KQK1': (0.99, 1.0, 0, 0.95, 0, HP, 0.997456, 1.0), 'KQK2': (0.99, 1.0, 0.95, 1.0, 0, 0.47123889803846897, 0.997456, 1.0),
        'KQB': (0.99, 1.0, 0.95, 1.0, 0.47123889803846897, 1.413716694115407, 0.997456, 1.0)}
SC = (1.0, 1.0, HP, 1.0)


def coeffs(e, l, w, t):
    e, l, w, t = [mp.mpf(v) for v in (e, l, w, t)]
    pz, qz = l * mp.sin(w), l * mp.cos(w)
    mu = mp.sqrt((1 - t) / t)
    N = 24; vals = []
    for n in range(N):
        f = 2 * mp.pi * n / N
        zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
        B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta); D = 1 + e * mp.cos(f); k = 2 * e * (1 + e)
        vals.append(mu ** 2 * k ** 2 * mp.sin(f) ** 2 * (1 - zeta ** 2) - D ** 2 * B ** 2)
    return [mp.fsum(vals[n] * mp.expj(-kk * 2 * mp.pi * n / N) for n in range(N)) / N for kk in range(-6, 7)]


def inside_count(c, rho):
    sc = max(abs(v) for v in c)
    lo = 0
    while abs(c[lo]) < mp.mpf(10) ** -80 * sc:
        lo += 1
    hi = 12
    while abs(c[hi]) < mp.mpf(10) ** -80 * sc:
        hi -= 1
    z = list(mp.polyroots(c[lo:hi + 1][::-1], maxsteps=800, extraprec=600)) if hi > lo else []
    marg = min([abs(abs(r) - rho) / rho for r in z] + [mp.inf])
    return lo + sum(1 for r in z if abs(r) < rho), marg


def rand_box(dom):
    bx = []
    for i in range(4):
        lo, hi = dom[2 * i], dom[2 * i + 1]
        wdt = min(SC[i] * 10 ** random.uniform(-5, -1), hi - lo)
        u = random.random()
        if u < 0.2:
            a = lo                          # lower face
        elif u < 0.4:
            a = hi - wdt                    # upper face (e = 1, t = 1, l = 1, w = pi/2 ...)
        else:
            a = random.uniform(lo, hi - wdt)
        bx += [a, a + wdt]
    return tuple(bx)


tot = {k: [0, 0, 0] for k in DOMS}          # boxes, certified, violations
minmarg = mp.inf
names = list(DOMS)
for it in range(NBOX):
    nm = names[it % len(names)]
    bx = rand_box(DOMS[nm])
    ok, T = kq_run.certify(bx)
    tot[nm][0] += 1
    if not ok:
        continue
    tot[nm][1] += 1
    pts = [[0.5 * (bx[2 * i] + bx[2 * i + 1]) for i in range(4)]]
    pts += [[bx[2 * i + (1 if (m >> i) & 1 else 0)] for i in range(4)] for m in random.sample(range(16), 4)]
    pts += [[random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(4)] for _ in range(3)]
    for x in pts:
        x[0] = min(x[0], 1 - 1e-12); x[3] = min(max(x[3], 1e-12), 1 - 1e-15)      # physical e < 1, 0 < t < 1
        c = coeffs(*x)
        sc = max(abs(v) for v in c); lo = 0
        while abs(c[lo]) < mp.mpf(10) ** -80 * sc:
            lo += 1
        hi = 12
        while abs(c[hi]) < mp.mpf(10) ** -80 * sc:
            hi -= 1
        z = list(mp.polyroots(c[lo:hi + 1][::-1], maxsteps=800, extraprec=600)) if hi > lo else []
        n = sum(1 for r in z if abs(mp.log(abs(r))) < mp.mpf(10) ** -30)
        minmarg = min(minmarg, mp.mpf(T - n))
        if n > T:
            tot[nm][2] += 1
            print('VIOLATION', nm, bx, x, 'certified total', T, 'N_mult', n, flush=True)
for nm in names:
    print(f"{nm}: random boxes {tot[nm][0]}, count-certified {tot[nm][1]}, violations {tot[nm][2]}", flush=True)
print('min (certified total - N_mult) over all checked points', mp.nstr(minmarg, 3), '; away count engine:', kq_run.away.__name__ if hasattr(kq_run.away, '__name__') else kq_run.away)

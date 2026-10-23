"""AUDIT-11 h21: the corner K of the face construction, with MY OWN P at 100 digits (tau = tan((f - pi)/2)).
 (a) the authors' finding: at sigma = 0.01877, E = 0.25, Pz = 0.9863, Qz = 0.1646 the FIP polynomial
     G(m) = Q(3 sigma; e = 1 - sigma^2 E, mu = sigma^2 m) changes sign for m in (0.995, 1.0]; and the zeros of P there.
 (b) case A points: (l, w) in K_lw, s = max(sqrt eps, sqrt mu) <= SF2 = 0.010005: near |tau| <= 3s (<= 5),
     gap 3s < |tau| < T_X2 = 3 SF2 (0), outer |tau| >= T_X2 (<= 6), total <= 10.
 (c) case B points: (l, w) in K_lw, SF2 < s <= 0.02: total <= 10 (whole-circle run KQS).
 (d) points just outside K_lw (l in [0.95, 1], w in [1.20, 1.29) or l in [0.85, 0.88), w in [1.29, 1.45)), s <= 0.02:
     the original construction needs gap 3s < |tau| < 0.06 empty, near <= 5, outer <= 6.
usage: python h21_corner_K.py NPTS SEED"""
import sys, math, random
import mpmath as mp
mp.mp.dps = 100
NPTS, SEED = int(sys.argv[1]), int(sys.argv[2]); random.seed(SEED)
SF2 = mp.mpf('0.010005'); W9 = 1.413716694115407


def P(e, pz, qz, mu, f):
    zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
    B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta); D = 1 + e * mp.cos(f); k = 2 * e * (1 + e)
    return mu ** 2 * k ** 2 * mp.sin(f) ** 2 * (1 - zeta ** 2) - D ** 2 * B ** 2


def taus(e, pz, qz, mu):
    N = 32
    vals = [P(e, pz, qz, mu, 2 * mp.pi * n / N) for n in range(N)]
    c = [mp.fsum(vals[n] * mp.expj(-kk * 2 * mp.pi * n / N) for n in range(N)) / N for kk in range(-6, 7)]
    sc = max(abs(v) for v in c); lo = 0
    while abs(c[lo]) < mp.mpf(10) ** -90 * sc:
        lo += 1
    hi = 12
    while abs(c[hi]) < mp.mpf(10) ** -90 * sc:
        hi -= 1
    z = mp.polyroots(c[lo:hi + 1][::-1], maxsteps=2000, extraprec=3000)
    return [mp.tan((mp.arg(r) % (2 * mp.pi) - mp.pi) / 2) for r in z if abs(mp.log(abs(r))) < mp.mpf(10) ** -12]


# (a)
sg, E_, pz, qz = mp.mpf('0.01877'), mp.mpf('0.25'), mp.mpf('0.9863'), mp.mpf('0.1646')
th = 3 * sg
Q = lambda m: P(1 - sg ** 2 * E_, pz, qz, sg ** 2 * m, mp.pi + 2 * mp.atan(th)) * (1 + th ** 2) ** 6
print(f'(a) G(m = 0.995) = {mp.nstr(Q(mp.mpf("0.995")), 6)}, G(m = 1.0) = {mp.nstr(Q(mp.mpf(1)), 6)}  (sign change: '
      f'{mp.sign(Q(mp.mpf("0.995"))) != mp.sign(Q(mp.mpf(1)))})')
m_ = mp.mpf('0.9975'); e_ = 1 - sg ** 2 * E_; mu_ = sg ** 2 * m_; s_ = mp.sqrt(max(sg ** 2 * E_, mu_))
tz = sorted(abs(t) for t in taus(e_, pz, qz, mu_))
print(f'    zeros at m = 0.9975: |tau| = {[mp.nstr(t, 5) for t in tz]}; s = {mp.nstr(s_, 5)}, gap = [3s, 0.06] = '
      f'[{mp.nstr(3 * s_, 5)}, 0.06]: zeros in the gap {sum(1 for t in tz if 3 * s_ < t < mp.mpf("0.06"))}')


def rand_eps_mu(smax, smin=mp.mpf(0)):
    while True:
        eps = mp.mpf(10) ** mp.mpf(random.uniform(-12, float(mp.log10(smax ** 2))))
        mu = mp.mpf(10) ** mp.mpf(random.uniform(-12, float(mp.log10(smax ** 2))))
        s = mp.sqrt(max(eps, mu))
        if smin < s <= smax:
            return eps, mu, s


def run(name, npts, lw, smax, smin, tx, near_t, gap_t, outer_t, whole_only=False):
    mx = [0, 0, 0, 0]; bad = 0; odd = 0
    for _ in range(npts):
        l, w = lw()
        eps, mu, s = rand_eps_mu(smax, smin)
        t = taus(1 - eps, mp.mpf(l) * mp.sin(w), mp.mpf(l) * mp.cos(w), mu)
        a = [abs(x) for x in t]
        cnt = [sum(1 for x in a if x <= 3 * s), sum(1 for x in a if 3 * s < x < tx), sum(1 for x in a if x >= tx), len(a)]
        mx = [max(u, v) for u, v in zip(mx, cnt)]; odd += cnt[3] % 2
        if (not whole_only and (cnt[0] > near_t or cnt[1] > gap_t or cnt[2] > outer_t)) or cnt[3] > 10:
            bad += 1; print('   exceeds a target:', name, l, w, float(eps), float(mu), cnt, flush=True)
    print(f'{name}: {npts} points; max near {mx[0]}, max gap {mx[1]}, max outer {mx[2]}, max total {mx[3]}; odd totals {odd}; '
          f'points exceeding a target {bad}', flush=True)


def in_K_lw():
    while True:                                   # bias to the corner l ~ 1, w ~ 1.3 .. 1.41 (outside the C3 cell)
        l = 1 - 10 ** random.uniform(-6, math.log10(0.12)); w = random.uniform(1.29, 1.5707963)
        if not (l >= 0.95 and w >= W9):
            return l, w


def outside_K_lw():
    if random.random() < 0.6:
        return 1 - 10 ** random.uniform(-6, math.log10(0.05)), random.uniform(1.20, 1.29)
    return random.uniform(0.85, 0.88), random.uniform(1.29, 1.45)


run('(b) case A, (l, w) in K_lw, s <= SF2', NPTS, in_K_lw, SF2, mp.mpf(0), 3 * SF2, 5, 0, 6)
run('(c) case B, (l, w) in K_lw, SF2 < s <= 0.02', NPTS // 2, in_K_lw, mp.mpf('0.02'), SF2, 3 * SF2, 99, 99, 99, whole_only=True)
run('(d) just outside K_lw, s <= 0.02 (original construction, T_X = 0.06)', NPTS, outside_K_lw, mp.mpf('0.02'), mp.mpf(0),
    mp.mpf('0.06'), 5, 0, 6)

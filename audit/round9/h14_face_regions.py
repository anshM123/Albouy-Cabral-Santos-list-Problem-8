"""AUDIT-9 h14: the three regions of the face construction at random physical points of the face neighbourhood
(eps, mu <= SF^2 = 4e-4, log-uniform down to 1e-12; l in [0, 1], w in [0, pi/2] outside the C3 cell), counted with MY OWN
P at 100 digits (zeros on the unit circle with multiplicity, tau = tan((f - pi)/2)):
  near  |tau| <= 3 s, s = max(sqrt eps, sqrt mu)   (FE/FM certify <= 5)
  gap   3 s < |tau| < T_X = 0.06                   (FIP/FIM certify 0)
  outer |tau| >= 0.06                              (KQF certifies <= 6)
Also the total (even, <= 10)."""
import random, math
import mpmath as mp
mp.mp.dps = 100
random.seed(14)
W9 = 1.413716694115407


def zeros(e, pz, qz, mu):
    N = 32; vals = []
    for n in range(N):
        f = 2 * mp.pi * n / N
        zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
        B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta); D = 1 + e * mp.cos(f); k = 2 * e * (1 + e)
        vals.append(mu ** 2 * k ** 2 * mp.sin(f) ** 2 * (1 - zeta ** 2) - D ** 2 * B ** 2)
    c = [mp.fsum(vals[n] * mp.expj(-kk * 2 * mp.pi * n / N) for n in range(N)) / N for kk in range(-6, 7)]
    sc = max(abs(v) for v in c); lo = 0
    while abs(c[lo]) < mp.mpf(10) ** -90 * sc:
        lo += 1
    hi = 12
    while abs(c[hi]) < mp.mpf(10) ** -90 * sc:
        hi -= 1
    z = mp.polyroots(c[lo:hi + 1][::-1], maxsteps=2000, extraprec=3000)
    return [r for r in z if abs(mp.log(abs(r))) < mp.mpf(10) ** -12]


mx = {'near': 0, 'gap': 0, 'outer': 0, 'total': 0}; odd = 0; n = 0; viol = 0
while n < 150:
    l = random.random() ** 0.3; w = random.uniform(0, math.pi / 2)
    if l >= 0.95 and w >= W9:
        continue
    eps = mp.mpf(10) ** mp.mpf(random.uniform(-12, math.log10(4e-4))); mu = mp.mpf(10) ** mp.mpf(random.uniform(-12, math.log10(4e-4)))
    s = mp.sqrt(max(eps, mu))
    z = zeros(1 - eps, mp.mpf(l) * mp.sin(w), mp.mpf(l) * mp.cos(w), mu)
    taus = [abs(mp.tan((mp.arg(r) % (2 * mp.pi) - mp.pi) / 2)) for r in z]
    cnt = {'near': sum(1 for t in taus if t <= 3 * s), 'gap': sum(1 for t in taus if 3 * s < t < mp.mpf('0.06')),
           'outer': sum(1 for t in taus if t >= mp.mpf('0.06')), 'total': len(taus)}
    for k in mx:
        mx[k] = max(mx[k], cnt[k])
    odd += cnt['total'] % 2
    if cnt['near'] > 5 or cnt['gap'] > 0 or cnt['outer'] > 6 or cnt['total'] > 10:
        viol += 1; print('  exceeds a certified target:', float(eps), float(mu), l, w, cnt, flush=True)
    n += 1
print(f'{n} random face-neighbourhood points: max near {mx["near"]} (target 5), max gap {mx["gap"]} (target 0), '
      f'max outer {mx["outer"]} (target 6), max total {mx["total"]}; odd totals {odd}; points exceeding a target {viol}')

"""AUDIT-11 h22 (light): the corner K, with my own P (as in h21; matches the physics, AUDIT_9 h9).
 (A) zeros of B near aphelion at e = 1 (the doc's 'B zero at phi0 = 2c^2/(3 Pz Qz)');
 (B) the FIP / FIM gap thresholds: sigma*(l, w) = smallest sigma at which P(tau = +-3 sigma; eps = sigma^2 E,
     mu = sigma^2 m) has a sign change over (E, m) in [0, 1] x [0, R1] (R1 = 1.001), i.e. the gap chart can no longer be
     zero-free; leading order: negative for m < 9c (c^2 = 1 - Pz^2);
 (C) zero positions for s <= SF2 at the worst corner of K outside the C3 cell (item 4 of the brief)."""
import math
import mpmath as mp
mp.mp.dps = 40
W9 = mp.mpf(1.413716694115407); R1 = mp.mpf('1.001'); SF2 = mp.mpf(0.010005); TX2 = mp.mpf(0.030015)


def P(e, pz, qz, mu, f):
    zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
    B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta); D = 1 + e * mp.cos(f); k = 2 * e * (1 + e)
    return mu ** 2 * k ** 2 * mp.sin(f) ** 2 * (1 - zeta ** 2) - D ** 2 * B ** 2


def Bf(e, pz, qz, f):
    zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
    return -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta)


print('(A) zeros of B(pi + phi) at e = 1, 0 < |phi| < 0.6 (tau = tan(phi/2)); doc: tau0 = c^2/(3 Pz Qz)')
for lab, pz, qz in (('l = 1, w = 9pi/20', mp.sin(W9), mp.cos(W9)), ('Pz = 0.9863, Qz = 0.1646', mp.mpf('0.9863'), mp.mpf('0.1646')),
                    ('l = 1, w = 1.35', mp.sin(1.35), mp.cos(1.35)), ('l = 0.97, w = 1.40', 0.97 * mp.sin(1.4), 0.97 * mp.cos(1.4))):
    pz, qz = mp.mpf(pz), mp.mpf(qz)
    g = [mp.mpf(-0.6) + mp.mpf(1.2) * k / 2400 for k in range(2401)]
    zs = []
    for a, b in zip(g, g[1:]):
        if abs(a) < 1e-9 or abs(b) < 1e-9:
            continue
        if mp.sign(Bf(1, pz, qz, mp.pi + a)) != mp.sign(Bf(1, pz, qz, mp.pi + b)):
            zs.append(mp.findroot(lambda x: Bf(1, pz, qz, mp.pi + x), (a, b), solver='bisect'))
    c2 = 1 - pz ** 2
    print(f'   {lab}: tau of the B zeros {[mp.nstr(mp.tan(z / 2), 5) for z in zs]}; doc tau0 = {mp.nstr(c2 / (3 * pz * qz), 4)};'
          f' pi/2 - w -> tau {mp.nstr(mp.tan((mp.pi / 2 - mp.atan2(pz, qz)) / 2), 5)}')


def gap_sign_change(pz, qz, sg, side):
    """True if P(tau = side 3 sg) takes both signs (or a zero) on the (E, m) grid"""
    f = mp.pi + 2 * mp.atan(side * 3 * sg)
    sg_ = set()
    for E in (0, mp.mpf('0.25'), mp.mpf('0.5'), 1):
        for m in (0, mp.mpf('0.5'), mp.mpf('0.9'), 1, R1):
            v = P(1 - sg ** 2 * E, pz, qz, sg ** 2 * m, f)
            sg_.add(int(mp.sign(v)))
    return len(sg_ - {0}) > 1 or 0 in sg_ or 1 in sg_      # leading-order sign is negative


print('(B) gap thresholds sigma* (first sigma in steps of 2.5e-4 with a nonnegative value on the (E, m) grid); need > SF2 = 0.010005 in K')
rows = []
for l in (mp.mpf(1), mp.mpf('0.999'), mp.mpf('0.99'), mp.mpf('0.97'), mp.mpf('0.9499'), mp.mpf('0.92')):
    for w in (mp.mpf('1.29'), mp.mpf('1.33'), mp.mpf('1.37'), mp.mpf('1.40'), W9 - mp.mpf('1e-6'), W9, mp.mpf('1.5')):
        incell = l >= mp.mpf('0.95') and w >= W9
        if incell:
            continue
        pz, qz = l * mp.sin(w), l * mp.cos(w)
        inK = pz >= mp.mpf('0.875') and mp.mpf('0.125') <= qz <= mp.mpf('0.25')
        res = []
        for side in (1, -1):
            s_star = None
            for k in range(1, 161):
                sg = mp.mpf('0.00025') * k
                if gap_sign_change(pz, qz, sg, side):
                    s_star = sg; break
            res.append(s_star)
        c = mp.sqrt(1 - pz ** 2)
        rows.append((l, w, inK, c, res))
        print(f'   l = {mp.nstr(l, 5)}, w = {mp.nstr(w, 8)}, in K {inK}, c = {mp.nstr(c, 5)}: sigma*_FIP = '
              f'{mp.nstr(res[0], 4) if res[0] else "> 0.04"}, sigma*_FIM = {mp.nstr(res[1], 4) if res[1] else "> 0.04"}; '
              f'estimate (c - 1/9)/3 = {mp.nstr((c - R1 / 9) / 3, 4)}', flush=True)
mn = min((r[4][0] for r in rows if r[2] and r[4][0] is not None), default=None)
print('   min sigma*_FIP over the grid points in K (outside the C3 cell):', mn)


def taus(e, pz, qz, mu):
    N = 32
    vals = [P(e, pz, qz, mu, 2 * mp.pi * n / N) for n in range(N)]
    c = [mp.fsum(vals[n] * mp.expj(-kk * 2 * mp.pi * n / N) for n in range(N)) / N for kk in range(-6, 7)]
    z = mp.polyroots(c[::-1], maxsteps=2000, extraprec=1500)
    return sorted(mp.tan((mp.arg(r) % (2 * mp.pi) - mp.pi) / 2) for r in z if abs(mp.log(abs(r))) < mp.mpf(10) ** -10)


mp.mp.dps = 60
print('(C) zeros (signed tau, |tau| < 0.25) at the worst corner of K outside the C3 cell, s <= SF2')
for l, w in ((mp.mpf(1), W9 - mp.mpf('1e-6')), (mp.mpf('0.99999'), W9 - mp.mpf('0.003')), (mp.mpf(1), mp.mpf('1.38'))):
    pz, qz = l * mp.sin(w), l * mp.cos(w)
    for s in (SF2, mp.mpf('0.008'), mp.mpf('0.005')):
        for E, m in ((0, 1), (mp.mpf('1e-6'), 1), (1, 1), (1, 0)):
            t = taus(1 - s ** 2 * E, pz, qz, s ** 2 * m)
            near = [x for x in t if abs(x) <= 3 * s]; gap = [x for x in t if 3 * s < abs(x) < TX2]
            outer = [x for x in t if abs(x) >= TX2]
            ymax = max((abs(x) / s for x in near), default=0)
            beyond = min((abs(x) for x in t if abs(x) > 3 * s), default=None)
            print(f'   l = {mp.nstr(l, 6)}, w = {mp.nstr(w, 7)}, s = {mp.nstr(s, 6)}, E = {mp.nstr(E, 3)}, m = {m}: total {len(t)},'
                  f' near {len(near)} (max |y| {mp.nstr(ymax, 4)}), gap {len(gap)}, outer {len(outer)}; smallest |tau| beyond 3s: '
                  f'{mp.nstr(beyond, 5)}; |tau| < 0.25: {[mp.nstr(x, 4) for x in t if abs(x) < 0.25]}', flush=True)

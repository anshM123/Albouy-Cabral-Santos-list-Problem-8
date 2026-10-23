"""Independent spot-check of disc-winding leaves (code '8', KQF files): at the centre and random points of sampled
leaves, map to physical (e, l, w, s), build the original a-chart polynomial Pol(u) at 60 digits (verify_mp.pol_count's
construction), and count its zeros on the unit circle z = e^{iu} (with multiplicity) whose true anomaly f
(tan(f/2) = sqrt((1+e)/(1-e)) tan(u/2)) lies in the outer region |th| >= T_X, th = tan((f - pi)/2), i.e.
|tan(f/2)| <= 1/T_X.  Checks that this outer count is <= 6.  Independent of the chart polynomials and of the engines;
a check at sample points, not a proof.
usage: python disc_mp_check.py FILE.json [--n 100] [--pts 3] [--seed 0]"""
import argparse
import math
import random
import mpmath as mp
from runio import load_json
from verify_mp import phys
from splitrule import split_box

mp.mp.dps = 60


def outer_count(e, l, w, s, tx):
    e, l, w, s = mp.mpf(e), mp.mpf(l), mp.mpf(w), mp.mpf(s)
    b = mp.sqrt(1 - e * e)
    p, q = l * mp.sin(w), l * mp.cos(w)
    n = 64
    vals = []
    for j in range(n):
        u = 2 * mp.pi * j / n
        r = 1 - e * mp.cos(u)
        Z = p * (mp.cos(u) - e) + q * b * mp.sin(u)
        Zp = -p * mp.sin(u) + q * b * mp.cos(u)
        A = 2 * e * mp.sin(u) * r
        rho2 = r * r - Z * Z
        d = A - 2 * Z * Zp
        vals.append(A * A * rho2 - s * s * d * d)
    cs = [mp.fsum(vals[j] * mp.expj(-2 * mp.pi * j * k / n) for j in range(n)) / n for k in range(-6, 7)]
    big = max(abs(c) for c in cs)
    while cs and abs(cs[-1]) < mp.mpf(10) ** (-45) * big:
        cs.pop()
    lead = 0
    while abs(cs[lead]) < mp.mpf(10) ** (-45) * big:
        lead += 1
    cs = cs[lead:]
    roots = mp.polyroots(cs[::-1], maxsteps=500, extraprec=600) if len(cs) > 1 else []
    k = mp.sqrt((1 + e) / (1 - e))
    out = 0
    for z in roots:
        if abs(mp.log(abs(z))) < mp.mpf(10) ** (-12):
            u = mp.arg(z)
            tf = k * mp.tan(u / 2)
            if abs(tf) <= 1 / mp.mpf(tx):
                out += 1
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('file'); ap.add_argument('--n', type=int, default=100); ap.add_argument('--pts', type=int, default=3)
    ap.add_argument('--seed', type=int, default=0)
    A = ap.parse_args()
    D = load_json(A.file)
    H = D['header']
    assert H['engine'] == 'KQF'
    tx = H['T_X']
    leaves = []
    for r in D['results']:
        stack = [tuple(r['box'])]
        for ch in r['tree']:
            bx = stack.pop()
            if ch == '+':
                if H.get('split_rule') is not None:
                    stack.extend(split_box(bx, H['scale'], H['split_rule']))
                else:
                    wd = [(bx[2 * i + 1] - bx[2 * i]) / H['scale'][i] for i in range(4)]
                    i = max(range(4), key=lambda k: wd[k])
                    mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
                    a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
                    stack.extend([tuple(a), tuple(c)])
            elif ch == '8':
                leaves.append(bx)
    random.seed(A.seed)
    sample = random.sample(leaves, min(A.n, len(leaves)))
    bad = 0; checked = 0; worst = 0
    for bx in sample:
        for k in range(A.pts + 1):
            pt = [0.5 * (bx[2 * i] + bx[2 * i + 1]) if k == 0 else random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(4)]
            ph = phys(dict(engine='KQ'), pt)
            if ph is None:
                continue
            N = outer_count(*ph, tx)
            checked += 1; worst = max(worst, N)
            if N > 6:
                bad += 1
                print('VIOLATION', bx, pt, ph, N, flush=True)
    print(f'{A.file}: {len(leaves)} disc leaves, sampled {len(sample)}, points checked {checked}, max outer count {worst}, '
          f'violations {bad}', flush=True)

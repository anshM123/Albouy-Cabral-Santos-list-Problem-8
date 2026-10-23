"""Independent spot-check of recorded leaves (bb_run / c3_run / kb_run JSON files).

For a random sample of certified leaves: map the box centre and random points of the box to physical parameters
(e, l, w, s) with s = R/a, build the ORIGINAL a-chart polynomial Pol(u) = A^2 rho^2 - s^2 ((rho^2)')^2 at 60 digits
(mpmath), count its roots on the unit circle (z = e^{iu}) with multiplicity, and check N <= recorded bound <= 10
(for C3 near charts, whose bound is a partial count, check the total <= 10).  This independent of every generated
chart polynomial and of the certificate engines; it checks the claim at sample points, it is not a proof.
usage: python verify_mp.py FILE.json [--n 200] [--pts 4] [--seed 0]"""
import argparse
import json
from runio import load_json
from splitrule import split_box
import math
import random
import mpmath as mp

mp.mp.dps = 60


def pol_count(e, l, w, s):
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
    while cs and abs(cs[-1]) < mp.mpf(10) ** (-45) * max(abs(c) for c in cs):
        cs.pop()                                  # degree drop: roots at infinity, off the circle
    lead = 0
    while abs(cs[lead]) < mp.mpf(10) ** (-45) * max(abs(c) for c in cs):
        lead += 1                                 # roots at 0, off the circle
    cs = cs[lead:]
    if len(cs) < 2:
        return 0
    roots = mp.polyroots(cs[::-1], maxsteps=500, extraprec=600)
    return sum(1 for z in roots if abs(mp.log(abs(z))) < mp.mpf(10) ** (-12))


def phys(H, bx_pt):
    """box point -> (e, l, w, s) physical, or None if outside the physical domain."""
    eng = H.get('engine') or H.get('chart')
    x0, x1, x2, x3 = bx_pt
    if eng == 'a':
        e = math.sin(x0); l = x1; w = x2; t = x3
        return (e, l, w, math.sqrt(t / (1 - t))) if 0 < t < 1 else None
    if eng == 'f':
        e = x0; l = x1; w = x2; t = x3
        return (e, l, w, math.sqrt(t / (1 - t)) * (1 - e)) if 0 < t < 1 and e < 1 else None
    if eng == 'KQ':
        e = x0; l = x1; w = x2; t = x3
        return (e, l, w, math.sqrt(t / (1 - t)) * (1 - e)) if 0 < t < 1 and e < 1 else None
    if eng in ('K', 'KA', 'KB'):
        e, l, w, mu = x0, x1, x2, x3
        return (e, l, w, (1 - e) / mu) if mu > 0 and e < 1 else None
    if eng == 'C1':
        e = x0 * math.cos(x1); L = x0 * math.sin(x1); t = x3
        return (e, math.sqrt(L), x2, math.sqrt(t / (1 - t))) if 0 < t < 1 and e > 0 else None
    if eng == 'C1p':
        e = x0; t = e * e * x3 / (1 - x3)
        return (e, x1, x2, math.sqrt(t / (1 - t))) if 0 < t < 1 and e > 0 else None
    if eng == 'C1pp':
        e = x0 * x1 * x1; t = x0 * x0 * x3 / (1 - x3)
        return (e, x1, x2, math.sqrt(t / (1 - t))) if 0 < t < 1 and e > 0 else None
    if eng == 'B1':
        e = x0; p = e * x1; q = e * x2; tau = -2.0 + 4.0 * x3; t = 0.5 + e * tau
        return (e, math.hypot(p, q), math.atan2(p, q), math.sqrt(t / (1 - t))) if e > 0 else None
    if eng == 'B2':
        th, eh, w = x0, x1, x2; e = th * eh; L = e * th; tau = 0.6825 * (2 * x3 - 1); t = 0.5 + th * tau
        return (e, math.sqrt(L), w, math.sqrt(t / (1 - t))) if e > 0 else None
    if eng in ('B3p', 'B3m'):
        sg = 1 if eng == 'B3p' else -1
        eta, eh, thh, w = x0, x1, x2, x3; e = eta * eh; L = e * eta * thh; t = 0.5 + sg * eta
        return (e, math.sqrt(L), w, math.sqrt(t / (1 - t))) if e > 0 else None
    if eng.startswith('L2') and eng[2] in 'EPAT' and len(eng) <= 4 and eng not in ('L2aE',):
        s0, s1, s2, s3 = 0.01, 0.25, 0.25, 0.5
        nu = x0
        if eng == 'L2E':
            E_, P_, A_, T_ = 1.0, x1, x2, -1.05 + 2.1 * x3
        elif eng == 'L2P':
            E_, P_, A_, T_ = x1, 1.0, x2, -1.05 + 2.1 * x3
        elif eng in ('L2Ap', 'L2Am'):
            E_, P_, A_, T_ = x1, x2, (1.0 if eng == 'L2Ap' else -1.0), -1.05 + 2.1 * x3
        else:
            E_, P_, A_, T_ = x1, x2, x3, (1.0 if eng == 'L2Tp' else -1.0)
        e = s0 * E_ * nu; ph = s1 * P_ * nu; qh = 1 + s2 * A_ * nu; tauh = s3 * T_ * nu
        p, q = e * ph, e * qh; t = 0.5 + e * tauh
        return (e, math.hypot(p, q), math.atan2(p, q), math.sqrt(t / (1 - t))) if e > 0 else None
    # C3 charts
    nu = x0
    if eng == 'Ca':
        qz, kap, eps, mu = nu, nu * nu * x1, nu * nu * x2, nu ** 3 * x3
    elif eng == 'CK':
        qz, kap, eps, mu = nu * x1, nu * nu, nu * nu * x2, nu ** 3 * x3
    elif eng == 'CE':
        qz, kap, eps, mu = nu * x1, nu * nu * x2, nu * nu, nu ** 3 * x3
    elif eng == 'Cm':
        qz, kap, eps, mu = nu * x1, nu * nu * x2, nu * nu * x3, nu ** 3
    else:
        return None
    pz = 1 - kap; l = math.hypot(pz, qz)
    if l > 1 or eps <= 0 or mu <= 0:
        return None
    e = 1 - eps
    return (e, l, math.atan2(pz, qz), eps / mu)


def leaves_of(D):
    H = D['header']
    scale = H['scale']
    out = []
    for r in D['results']:
        stack = [tuple(r['box'])]
        for ch in r['tree']:
            bx = stack.pop()
            if ch == '+' and H.get('split_rule') is not None:     # kq_run --facesplit / --splitrule (splitrule.py)
                stack.extend(split_box(bx, scale, H['split_rule']))
            elif ch == '+':
                w = [(bx[2 * i + 1] - bx[2 * i]) / scale[i] for i in range(4)]
                i = max(range(4), key=lambda k: w[k])
                mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
                a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
                stack.extend([tuple(a), tuple(c)])
            elif ch in 'abcdefghijk':
                out.append((bx, ord(ch) - ord('a')))
            elif ch in 'ABCDEFGHIJK':
                out.append((bx, ord(ch) - ord('A')))
            elif ch in 'lL':
                out.append((bx, 10))          # total 11 certified; the even true count is <= 10
            elif ch == '7':
                out.append((bx, 10))          # tau-winding certificate (tau_wind.py): <= 10
            elif ch in '123456':
                out.append((bx, 12 - 2 * (ord(ch) - ord('0'))))
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('file'); ap.add_argument('--n', type=int, default=200); ap.add_argument('--pts', type=int, default=4)
    ap.add_argument('--seed', type=int, default=0)
    A = ap.parse_args()
    D = load_json(A.file)
    H = D['header']
    eng = H.get('engine') or H.get('chart')
    partial = eng in ('Ca', 'CK', 'CE', 'Cm', 'KA')    # partial-arc bound only: check total <= 10
    L = leaves_of(D)
    random.seed(A.seed)
    sample = random.sample(L, min(A.n, len(L)))
    bad = 0; checked = 0; worst = 0
    for bx, bound in sample:
        for k in range(A.pts + 1):
            pt = [0.5 * (bx[2 * i] + bx[2 * i + 1]) if k == 0 else random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(4)]
            ph = phys(H, pt)
            if ph is None:
                continue
            N = pol_count(*ph)
            checked += 1
            worst = max(worst, N)
            lim = 10 if partial else min(bound, 10)
            if N > lim:
                bad += 1
                print('VIOLATION', eng, bx, 'bound', bound, 'point', pt, 'phys', ph, 'count', N, flush=True)
    print(f'{A.file}: engine {eng}, {len(L)} leaves, sampled {len(sample)}, physical points checked {checked}, '
          f'max count {worst}, violations {bad}', flush=True)

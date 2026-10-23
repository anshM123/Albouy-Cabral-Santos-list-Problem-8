"""AUDIT-3 c2: own replay of the c3_run.py trees and 50-digit spot checks of C3 leaves.

Replay: same split rule (header scale); codes '+' split, '!' fail, 'u'/'o' pruned, 'a'+c certified (bound c),
'A'+c zone leaf (outer count c); checks that the tree is consumed exactly, fails/left consistent, volume conserved.

Spot check (independent of gen_c3 / nbcoef): for sampled leaves and sampled box points, map to the physical
(e, Pz, Qz, mu) with MY OWN implementation of the chart maps of C3_CORE.md, evaluate
  P(f) = mu^2 k^2 sin^2 f sigma^2 - D^2 B^2,  f = pi + 2 atan(tau),  Q(tau) = P (1 + tau^2)^6
at 60 digits, rebuild the chart polynomial in the chart's root variable by interpolation (15 Chebyshev nodes, degree 12
fit, residual checked), and count its real roots (with multiplicity) on the certified root interval:
  level-1 Ca/CK/CE/Cm: y = tau/nu in [-3, 3]                <= bound (zone leaves: |y| in [0.25, 3] <= c, inner <= 4)
  PhiP/PhiM:           root m in [0, R1^3] at tau = +-3 nu   == 0
  L2*E / L2*m:         Y in [-R2, R2], tau = nu Y2 eta Y      <= bound
  L2*Yp / L2*Ym:       root mp in [0, R2^2] at y = +-Y2 eta   == 0
usage: python c2_c3_replay_spot.py NLEAF NPT SEED RUN.json [RUN.json ...]
"""
import sys
import json
import random
import numpy as np
import mpmath as mp

mp.mp.dps = 100
NLEAF = int(sys.argv[1]); NPT = int(sys.argv[2]); random.seed(int(sys.argv[3])); files = sys.argv[4:]
Y0, Y2, E2, M2 = mp.mpf(3), mp.mpf('0.25'), mp.mpf('0.01'), mp.mpf('0.01')
R1 = R2 = mp.mpf('1.001')


def split(bx, S):
    w = [(bx[2 * i + 1] - bx[2 * i]) / S[i] for i in range(4)]
    i = int(np.argmax(w))
    mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)


def vol(b):
    return float(np.prod([b[2 * i + 1] - b[2 * i] for i in range(4)]))


def P_phys(e, pz, qz, mu, f):
    zeta = pz * mp.cos(f) + qz * mp.sin(f)
    dz = -pz * mp.sin(f) + qz * mp.cos(f)
    B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta)
    D = 1 + e * mp.cos(f)
    k = 2 * e * (1 + e)
    return mu ** 2 * k ** 2 * mp.sin(f) ** 2 * (1 - zeta ** 2) - D ** 2 * B ** 2


def Q_phys(e, pz, qz, mu, tau):
    return P_phys(e, pz, qz, mu, mp.pi + 2 * mp.atan(tau)) * (1 + tau ** 2) ** 6


def chart_fn(ch, x):
    """returns (function r -> chart value, root interval, kind) for box point x (my own maps)."""
    x = [mp.mpf(v) for v in x]
    nu = x[0]
    if ch in ('Ca', 'CK', 'CE', 'Cm'):
        if ch == 'Ca':
            qz, kap, ep, mu = nu, nu ** 2 * x[1], nu ** 2 * x[2], nu ** 3 * x[3]
        elif ch == 'CK':
            qz, kap, ep, mu = nu * x[1], nu ** 2, nu ** 2 * x[2], nu ** 3 * x[3]
        elif ch == 'CE':
            qz, kap, ep, mu = nu * x[1], nu ** 2 * x[2], nu ** 2, nu ** 3 * x[3]
        else:
            qz, kap, ep, mu = nu * x[1], nu ** 2 * x[2], nu ** 2 * x[3], nu ** 3
        return (lambda r: Q_phys(1 - ep, 1 - kap, qz, mu, nu * r) / nu ** 10), (-Y0, Y0), 'count'
    if ch in ('PhiP', 'PhiM'):
        sg = 1 if ch == 'PhiP' else -1
        a_, K_, E_ = x[1], x[2], x[3]
        return (lambda r: Q_phys(1 - nu ** 2 * E_, 1 - nu ** 2 * K_, nu * a_, nu ** 3 * r, sg * Y0 * nu) / nu ** 10), \
            (mp.mpf(0), R1 ** 3), 'nozero'
    # level 2
    rv, eta, w4 = x[1], x[2], x[3]
    if ch.startswith('L2a'):
        qz, kap = nu, nu ** 2 * rv
    else:
        qz, kap = nu * rv, nu ** 2
    kind = ch[3:]
    if kind == 'E':
        f = lambda r: Q_phys(1 - nu ** 2 * E2 * eta ** 2, 1 - kap, qz, nu ** 3 * M2 * eta ** 2 * w4, nu * Y2 * eta * r)
        return (lambda r: f(r) / (nu ** 10 * eta ** 6)), (-R2, R2), 'count'
    if kind == 'm':
        f = lambda r: Q_phys(1 - nu ** 2 * E2 * eta ** 2 * w4, 1 - kap, qz, nu ** 3 * M2 * eta ** 2, nu * Y2 * eta * r)
        return (lambda r: f(r) / (nu ** 10 * eta ** 6)), (-R2, R2), 'count'
    sg = 1 if kind == 'Yp' else -1
    f = lambda r: Q_phys(1 - nu ** 2 * E2 * eta ** 2 * w4, 1 - kap, qz, nu ** 3 * M2 * eta ** 2 * r, sg * nu * Y2 * eta)
    return (lambda r: f(r) / (nu ** 10 * eta ** 6)), (mp.mpf(0), R2 ** 2), 'nozero'


def real_roots_in(fn, lo, hi, deg=12):
    """interpolate fn (a polynomial of degree <= deg in r) on Chebyshev nodes of a wider interval; return real roots
    in [lo, hi] (with multiplicity) and the max interpolation residual at check points."""
    c = (lo + hi) / 2; h = (hi - lo) / 2 * mp.mpf('1.5') + mp.mpf('0.5')
    n = deg + 1
    nodes = [c + h * mp.cos(mp.pi * (2 * i + 1) / (2 * n)) for i in range(n)]
    vals = [fn(t) for t in nodes]
    V = mp.matrix([[t ** j for j in range(n)] for t in nodes])
    coef = mp.lu_solve(V, mp.matrix(vals))
    cl = [coef[j] for j in range(n)]
    # residual check at extra points
    res = max(abs(fn(t) - sum(cl[j] * t ** j for j in range(n))) for t in [c + h * mp.mpf(k) / 7 for k in (-6, -3, 1, 5)])
    scale = max(abs(v) for v in vals) + mp.mpf(10) ** -50
    while len(cl) > 1 and abs(cl[-1]) < mp.mpf(10) ** -40 * scale:
        cl.pop()
    if len(cl) <= 1:
        return [], res / scale, cl
    roots = mp.polyroots(cl[::-1], maxsteps=300, extraprec=200)
    rr = [r for r in roots if abs(mp.im(r)) < mp.mpf(10) ** -12 * (1 + abs(r))]
    inside = [mp.re(r) for r in rr if lo - mp.mpf(10) ** -12 <= mp.re(r) <= hi + mp.mpf(10) ** -12]
    return inside, res / scale, cl


for fn in files:
    D = json.load(open(fn)); H = D['header']; ch = H['chart']; S = H['scale']
    leaves = []; bad = 0; nf = 0; nl = 0
    for r in D['results']:
        root = tuple(r['box']); stack = [root]; v = 0.0; fi = 0; cnt = 0
        fl = [tuple(f[0]) for f in r['fails']]
        for code in r['tree']:
            bx = stack.pop()
            if code == '+':
                stack.extend(split(bx, S))
            elif code == '!':
                if fi >= len(fl) or fl[fi] != bx:
                    bad += 1
                fi += 1; v += vol(bx)
            else:
                cnt += 1; v += vol(bx)
                if code in 'uo':
                    leaves.append((bx, code, None))
                elif 'a' <= code <= 'k':
                    leaves.append((bx, 'a', ord(code) - ord('a')))
                elif 'A' <= code <= 'E':
                    leaves.append((bx, 'Z', ord(code) - ord('A')))
                else:
                    bad += 1
        if r['done'] and stack:
            bad += 1
        if (not r['done']) and stack != [tuple(b) for b in r['left']]:
            bad += 1
        v += sum(vol(b) for b in r['left'])
        if abs(v - vol(root)) > 1e-9 * vol(root) or cnt != r['nleaf'] or fi != len(fl):
            bad += 1
        nf += len(fl); nl += cnt
    kinds = {}
    for _, k, b in leaves:
        key = k if b is None else f'{k}{b}'
        kinds[key] = kinds.get(key, 0) + 1
    print(f"{fn.split('/')[-1]}: records {len(D['results'])}, leaves {nl}, fails {nf}, replay problems {bad}; leaf kinds {dict(sorted(kinds.items()))}",
          flush=True)
    # spot checks
    cert = [l for l in leaves if l[1] in ('a', 'Z')]
    sample = random.sample(cert, min(NLEAF, len(cert)))
    viol = 0; worst_res = mp.mpf(0); checked = 0; maxcnt = 0
    for bx, kind, bound in sample:
        for k in range(NPT):
            x = [bx[2 * i] + (bx[2 * i + 1] - bx[2 * i]) * random.random() for i in range(4)]
            if x[0] < 1e-4:
                x[0] = bx[0] + (bx[1] - bx[0]) * 0.5 + 1e-4       # keep nu away from 0 (blow-up face)
                if x[0] > bx[1]:
                    continue
            if ch.startswith('L2') and x[2] < 1e-3:
                if bx[5] < 1e-3:
                    continue
                x[2] = max(bx[4], 1e-3)
            f, (lo, hi), mode = chart_fn(ch, x)
            if kind == 'Z':
                ins, res, _ = real_roots_in(f, -Y0, Y0)
                outer = [r for r in ins if abs(r) >= Y2 - mp.mpf(10) ** -12]
                inner = [r for r in ins if abs(r) <= Y2 + mp.mpf(10) ** -12]
                n = len(outer)
                if len(outer) > bound or len(inner) > 4:
                    viol += 1; print('   VIOLATION (zone)', ch, bx, x, len(outer), len(inner), bound)
            else:
                ins, res, _ = real_roots_in(f, lo, hi)
                n = len(ins)
                if (mode == 'count' and n > bound) or (mode == 'nozero' and n > 0):
                    viol += 1; print('   VIOLATION', ch, bx, x, n, bound, [mp.nstr(r, 8) for r in ins])
            maxcnt = max(maxcnt, n)
            worst_res = max(worst_res, res); checked += 1
    print(f"    spot check: {len(sample)} leaves x {NPT} pts = {checked} points; violations {viol}; max real-root count on "
          f"the certified interval {maxcnt}; max interpolation residual (rel) {mp.nstr(worst_res, 3)}", flush=True)

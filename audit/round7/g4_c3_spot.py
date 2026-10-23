"""AUDIT-7 g4 (final pass; from AUDIT-3 c2): own replay of the c3_run.py trees, 100-digit spot checks of C3 leaves,
and re-certification of sampled leaves with c3_run.certify (same kind, bound <= recorded).
Changes against c2: zone codes 'A'..'F' (ZONE_OUTER = 5 since round 6), gzip-aware loading, records.jsonl input
(completed records only, chart from the file name), re-certification.

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
usage: python g4_c3_spot.py NLEAF NPT SEED NRECERT RUN.json[.gz|.records.jsonl] [...]
"""
import os
import sys
import json
import gzip
import random
import numpy as np
import mpmath as mp

mp.mp.dps = 100
NLEAF = int(sys.argv[1]); NPT = int(sys.argv[2]); random.seed(int(sys.argv[3])); NREC = int(sys.argv[4])
files = [os.path.abspath(f) for f in sys.argv[5:]]
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
C3DIR = os.path.normpath(os.path.join(HERE, '..', '..', 'code', 'c3'))
sys.path.insert(0, os.path.dirname(C3DIR)); sys.path.insert(0, C3DIR)
ZONE_MAX = 'F'


def load_any(fn):
    if fn.endswith('.records.jsonl'):
        recs = []
        for line in open(fn):
            try:
                r = json.loads(line)
            except Exception:
                continue
            if r.get('done'):
                recs.append(r)
        import c3_run
        ch = os.path.basename(fn).split('.')[0].split('_')[0]
        d = c3_run.DOM[ch]
        sc_ = [d[2 * i + 1] - d[2 * i] for i in range(4)]
        if os.environ.get('AUDIT_SCALE'):              # AUDIT-13: continuation levels run with their own --scale
            sc_ = [float(v) for v in os.environ['AUDIT_SCALE'].split(',')]
        else:
            print('   WARNING partial file without header: leaves rebuilt with the base domain widths as scale '
                  '(set AUDIT_SCALE for a continuation level)', flush=True)
        return dict(header=dict(chart=ch, scale=sc_), results=recs)
    p = fn if os.path.exists(fn) else fn + '.gz'
    return json.load(gzip.open(p, 'rt') if p.endswith('.gz') else open(p))
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
    if ch in ('FE', 'FM', 'FIP', 'FIM'):          # round 9 face charts (weights th:1, eps:2, mu:2), G = Q / s^6
        sv, rat, pzv, qzv = x
        if ch == 'FE':
            return (lambda r: Q_phys(1 - sv ** 2, pzv, qzv, sv ** 2 * rat, sv * r) / sv ** 6), (-mp.mpf(3), mp.mpf(3)), 'count'
        if ch == 'FM':
            return (lambda r: Q_phys(1 - sv ** 2 * rat, pzv, qzv, sv ** 2, sv * r) / sv ** 6), (-mp.mpf(3), mp.mpf(3)), 'count'
        sg = 3 if ch == 'FIP' else -3
        return (lambda r: Q_phys(1 - sv ** 2 * rat, pzv, qzv, sv ** 2 * r, sg * sv) / sv ** 6), (mp.mpf(0), R1), 'nozero'
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
    if not (os.path.exists(fn) or os.path.exists(fn + '.gz')):
        print(os.path.basename(fn), 'not present', flush=True); continue
    D = load_any(fn); H = D['header']; ch = H['chart']; S = H['scale']
    if H.get('split_rule'):                   # AUDIT-12: g4 replays the plain c3_run bisection only
        print('   PROBLEM: header split_rule %r not supported by g4 (leaf boxes would be wrong)' % H['split_rule'], flush=True)
    leaves = []; bad = 0; nf = 0; nl = 0
    for r in D['results']:
        root = tuple(r['box']); stack = [root]; v = 0.0; fi = 0; cnt = 0
        twr = r.get('tw', []); ti = 0                 # round 13: tau-winding contours of this record
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
                elif 'A' <= code <= ZONE_MAX and ch in ('Ca', 'CK'):
                    leaves.append((bx, 'Z', ord(code) - ord('A')))
                elif code == '7' and ch in ('Ca', 'CK', 'CE', 'Cm'):     # round 13: tau-winding, whole circle <= 10
                    if ti < len(twr) and len(twr[ti]) == 4 and twr[ti][0] == 0 and twr[ti][2] > twr[ti][3] > 0:
                        leaves.append((bx, 'T', tuple(twr[ti])))
                    else:
                        bad += 1
                    ti += 1
                else:
                    bad += 1
        if r['done'] and stack:
            bad += 1
        if (not r['done']) and stack != [tuple(b) for b in r['left']]:
            bad += 1
        v += sum(vol(b) for b in r['left'])
        if ti != len(twr):
            bad += 1
        if abs(v - vol(root)) > 1e-9 * vol(root) or cnt != r['nleaf'] or fi != len(fl):
            bad += 1
        nf += len(fl); nl += cnt
    kinds = {}
    for _, k, b in leaves:
        key = k if b is None else ('T' if k == 'T' else f'{k}{b}')
        kinds[key] = kinds.get(key, 0) + 1
    print(f"{os.path.basename(fn)}: records {len(D['results'])}, leaves {nl}, fails {nf}, replay problems {bad}; leaf kinds {dict(sorted(kinds.items()))}",
          flush=True)
    # spot checks
    cert = [l for l in leaves if l[1] in ('a', 'Z')]
    sample = random.sample(cert, min(NLEAF, len(cert)))
    viol = 0; worst_res = mp.mpf(0); checked = 0; maxcnt = 0
    for bx, kind, bound in sample:
        for k in range(NPT):
            x = [bx[2 * i] + (bx[2 * i + 1] - bx[2 * i]) * random.random() for i in range(4)]
            if x[0] < (1e-7 if ch.startswith('F') else 1e-4):
                x[0] = bx[0] + (bx[1] - bx[0]) * 0.5 + (1e-7 if ch.startswith('F') else 1e-4)       # keep nu away from 0 (blow-up face)
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
    # round 13: tau-winding leaves of the C3 charts
    tws = [l for l in leaves if l[1] == 'T']
    if tws:
        import tau_wind
        tsamp = random.sample(tws, min(NREC if NREC else NLEAF, len(tws)))
        tbad = sum(1 for bx, _, C_ in tsamp if not tau_wind.check_tw_chart(ch, bx, C_[1], C_[2], C_[3]))
        tviol = 0; tmax = 0
        for bx, _, C_ in random.sample(tws, min(NLEAF, len(tws))):
            c0 = mp.mpc(C_[1], C_[2]); cnts = set()
            for k in range(NPT):
                x = [mp.mpf(bx[2 * i] + (bx[2 * i + 1] - bx[2 * i]) * random.random()) for i in range(4)]
                if x[0] <= 0:
                    continue
                nu = x[0]
                if ch == 'Ca':
                    qz, kap, ep, mu = nu, nu ** 2 * x[1], nu ** 2 * x[2], nu ** 3 * x[3]
                elif ch == 'CK':
                    qz, kap, ep, mu = nu * x[1], nu ** 2, nu ** 2 * x[2], nu ** 3 * x[3]
                elif ch == 'CE':
                    qz, kap, ep, mu = nu * x[1], nu ** 2 * x[2], nu ** 2, nu ** 3 * x[3]
                else:
                    qz, kap, ep, mu = nu * x[1], nu ** 2 * x[2], nu ** 2 * x[3], nu ** 3
                vals = [P_phys(1 - ep, 1 - kap, qz, mu, 2 * mp.pi * n / 32) for n in range(32)]
                cf = [mp.fsum(vals[n] * mp.expj(-kk * 2 * mp.pi * n / 32) for n in range(32)) / 32 for kk in range(-6, 7)]
                co = cf[::-1]
                while abs(co[0]) < mp.mpf(10) ** -50 * max(abs(t_) for t_ in co):
                    co = co[1:]
                zz = mp.polyroots(co, maxsteps=3000, extraprec=1500)
                Np = sum(1 for z in zz if abs(abs(z) - 1) < mp.mpf(10) ** -20)
                ys = [(-1j * (z + 1) / (z - 1)) / nu for z in zz if abs(z - 1) > mp.mpf(10) ** -40]
                ins = sum(1 for yv in ys if abs(yv - c0) < C_[3])
                cnts.add(ins); tmax = max(tmax, Np)
                if Np > 10 or ins < 1:
                    tviol += 1; print('   TW SPOT VIOLATION', ch, bx, [float(v_) for v_ in x], Np, ins, C_, flush=True)
            if len(cnts) > 1:
                tviol += 1; print('   TW SPOT VIOLATION (inside count not constant)', ch, bx, C_, cnts, flush=True)
        print(f"    tau-winding leaves {len(tws)}: re-certified {len(tsamp)} with tau_wind.check_tw_chart (problems {tbad}); "
              f"100-digit spot: violations {tviol}, max N {tmax}", flush=True)
    # re-certification with the production code (c3_run.certify): same kind, bound <= recorded
    if NREC:
        import c3_run
        c3_run.init(ch)
        if hasattr(c3_run, 'DEEP_ALWAYS'):
            c3_run.DEEP_ALWAYS = True                # round 8: two-tier counts; allow the deep pass
        pool = [l for l in leaves if l[1] != 'T']
        rs = random.sample(pool, min(NREC, len(pool)))
        rbad = 0
        for bx, kind, bound in rs:
            k2, c2 = c3_run.certify(bx)
            if kind in ('u', 'o'):
                okk = (k2 == kind)
            elif kind == 'a':
                okk = (k2 == 'a' and c2 <= bound)
            else:
                okk = (k2 == 'Z' and c2 <= bound)
            if not okk:
                rbad += 1; print('   RECERT PROBLEM', ch, bx, kind, bound, k2, c2, flush=True)
        print(f"    re-certified {len(rs)} leaves (incl. pruned) with c3_run.certify: problems {rbad}", flush=True)

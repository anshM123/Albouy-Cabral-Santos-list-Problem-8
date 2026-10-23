"""AUDIT-7 g3 (final pass): 100-digit spot checks of sampled leaves of every whole-circle family and of KA, with MY OWN
chart maps (audit round2/b3 for B*/L2*, round6/f3 for D2*, round4/d2 for a/C1*) and MY OWN closed forms:
  a-chart based charts (a, C1, C1p, C1pp, B*, L2*, D2*): f_t(u) = (1-t) A^2 rho^2 - t ((rho^2)')^2 at the mapped
      (e, p, q, t); the chart polynomial is f_t divided by a positive normalisation (same zeros);
  f-chart and KQ: t P(f) = (1-t) k^2 sin^2 f sigma^2 - t D^2 B^2  (= the f-chart polynomial, AUDIT_6);
  KA: zeros of P(f) on the away arc |f - pi| >= phi_a (header), bound <= 2.
At NPT points (centre, random) of NLEAF sampled leaves per family (winding leaves, total-11 leaves and the rest
sampled separately): N_mult = zeros on |z| = 1 with multiplicity (from the degree-12 self-inversive polynomial) must be
<= min(bound, 10) and even; for winding leaves #{|z| < e^L} must equal W (and the margin of the nearest zero to the
test circle is reported).  Points with nu = 0 / e = 0 (blow-up faces) are moved inside the box.
usage: python g3_spot_all.py NLEAF NPT SEED NAME [NAME ...]   (NAME = bulk, C1, ..., F0, KQ1, ..., KA)"""
import os, sys, json, gzip, math, random
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.normpath(os.path.join(HERE, '..', '..', 'code'))
NLEAF, NPT, SEED = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
names = sys.argv[4:]
mp.mp.dps = 100
random.seed(SEED)
sys.path.insert(0, CODE)
import splitrule
S0, S1, S2, S3 = mp.mpf('0.01'), mp.mpf('0.25'), mp.mpf('0.25'), mp.mpf('0.5')
SD = mp.mpf(5) / 1000


def load(name):
    for p in (os.path.join(CODE, 'runs', name + '.json'), os.path.join(CODE, 'runs', name + '.json.gz')):
        if os.path.exists(p):
            return json.load(gzip.open(p, 'rt') if p.endswith('.gz') else open(p))
    p = os.path.join(CODE, 'runs', name + '.json.records.jsonl')        # partial run: completed records only
    if os.path.exists(p):
        recs = []
        for line in open(p):
            try:
                r = json.loads(line)
            except Exception:
                continue
            if r.get('done'):
                recs.append(r)
        return dict(partial=True, header=None, results=recs)
    return None


def a_map(eng, x):
    """box point -> (e, p, q, t) for the a-chart closed forms (t may exceed 1 in C1p/C1pp/D2: positive factor kept)"""
    x = [mp.mpf(v) for v in x]
    if eng == 'a':
        eta, l, w, t = x; return mp.sin(eta), l * mp.sin(w), l * mp.cos(w), t
    if eng == 'C1':
        eps, th, w, t = x; e = eps * mp.cos(th); l = mp.sqrt(eps * mp.sin(th)); return e, l * mp.sin(w), l * mp.cos(w), t
    if eng == 'C1p':
        eps, l, w, tau = x; return eps, l * mp.sin(w), l * mp.cos(w), eps ** 2 * tau / (1 - tau)
    if eng == 'C1pp':
        epp, l, w, tau = x; e = epp * l * l; return e, l * mp.sin(w), l * mp.cos(w), epp ** 2 * tau / (1 - tau)
    if eng == 'B1':
        e, ph, qh, x3 = x; return e, e * ph, e * qh, mp.mpf(1) / 2 + e * (-2 + 4 * x3)
    if eng == 'B2':
        th, eh, w, x3 = x; e = th * eh; l = mp.sqrt(e * th)
        return e, l * mp.sin(w), l * mp.cos(w), mp.mpf(1) / 2 + th * (-mp.mpf('0.6825') + 2 * mp.mpf('0.6825') * x3)
    if eng in ('B3p', 'B3m'):
        sg = 1 if eng == 'B3p' else -1
        eta, eh, thh, w = x; e = eta * eh; l = mp.sqrt(e * eta * thh)
        return e, l * mp.sin(w), l * mp.cos(w), mp.mpf(1) / 2 + sg * eta
    if eng.startswith('L2'):
        nu = x[0]
        if eng == 'L2E':
            E, P, A, T = 1, x[1], x[2], -mp.mpf('1.05') + mp.mpf('2.1') * x[3]
        elif eng == 'L2P':
            E, P, A, T = x[1], 1, x[2], -mp.mpf('1.05') + mp.mpf('2.1') * x[3]
        elif eng in ('L2Ap', 'L2Am'):
            E, P, A, T = x[1], x[2], (1 if eng == 'L2Ap' else -1), -mp.mpf('1.05') + mp.mpf('2.1') * x[3]
        else:
            E, P, A, T = x[1], x[2], x[3], (1 if eng == 'L2Tp' else -1)
        e = S0 * E * nu; ph = S1 * P * nu; qh = 1 + S2 * A * nu; tauh = S3 * T * nu
        return e, e * ph, e * qh, mp.mpf(1) / 2 + e * tauh
    if eng.startswith('D2'):
        if eng in ('D2Tp', 'D2Tm'):
            n, E, P, K = x; sg = 1 if eng == 'D2Tp' else -1
            e, p, k, tp = SD * E * n, SD * P * n, SD * K * n, sg * SD * n
        else:
            n, r1, r2, x3 = x
            T = (mp.mpf('1.04') + mp.mpf('0.24') * x3) if eng == 'D2E2' else (-mp.mpf('1.05') + mp.mpf('2.1') * x3)
            if eng in ('D2E', 'D2E2'):
                e, p, k = SD * n, SD * r1 * n, SD * r2 * n
            elif eng == 'D2P':
                e, p, k = SD * r1 * n, SD * n, SD * r2 * n
            else:
                e, p, k = SD * r1 * n, SD * r2 * n, SD * n
            tp = SD * T * n
        tau = mp.mpf(1) / 2 + tp
        return e, p, 1 - k, e * e * tau / (1 - tau)
    raise ValueError(eng)


def coeffs_a(e, p, q, t):
    b = mp.sqrt(1 - e * e); N = 32; vals = []
    for n in range(N):
        u = 2 * mp.pi * n / N
        r = 1 - e * mp.cos(u); Z = p * (mp.cos(u) - e) + q * b * mp.sin(u); Zp = -p * mp.sin(u) + q * b * mp.cos(u)
        A = 2 * e * mp.sin(u) * r; rho2 = r * r - Z * Z; d = A - 2 * Z * Zp
        vals.append((1 - t) * A * A * rho2 - t * d * d)
    return [mp.fsum(vals[n] * mp.expj(-k * 2 * mp.pi * n / N) for n in range(N)) / N for k in range(-6, 7)]


def tP(e, l, w, t, f):
    pz, qz = l * mp.sin(w), l * mp.cos(w)
    zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
    B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta); D = 1 + e * mp.cos(f); k = 2 * e * (1 + e)
    return (1 - t) * k ** 2 * mp.sin(f) ** 2 * (1 - zeta ** 2) - t * D ** 2 * B ** 2


def coeffs_f(e, l, w, t):
    N = 32
    vals = [tP(e, l, w, t, 2 * mp.pi * n / N) for n in range(N)]
    return [mp.fsum(vals[n] * mp.expj(-k * 2 * mp.pi * n / N) for n in range(N)) / N for k in range(-6, 7)]


def zeros(c):
    sc = max(abs(v) for v in c)
    lo = 0
    while abs(c[lo]) < mp.mpf(10) ** -80 * sc:
        lo += 1
    hi = 12
    while abs(c[hi]) < mp.mpf(10) ** -80 * sc:
        hi -= 1
    z = list(mp.polyroots(c[lo:hi + 1][::-1], maxsteps=800, extraprec=600)) if hi > lo else []
    return z + [mp.mpf(0)] * lo


SPLIT_RULE = None                                    # header split_rule of the file being replayed (round 9)


def split(bx, S):
    if SPLIT_RULE:
        return splitrule.split_box(bx, S, SPLIT_RULE)
    w = [(bx[2 * i + 1] - bx[2 * i]) / S[i] for i in range(4)]
    i = max(range(4), key=lambda j: w[j]); mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)


def leaves(D, S, n):
    """reservoir samples (size n) of winding leaves, total-11 leaves and other leaves, plus the class sizes"""
    res = {'w': [], 'e': [], 'r': []}; cnt = {'w': 0, 'e': 0, 'r': 0}

    def keep(c, item):
        cnt[c] += 1
        if len(res[c]) < n:
            res[c].append(item)
        else:
            j = random.randrange(cnt[c])
            if j < n:
                res[c][j] = item
    for r in D['results']:
        stack = [tuple(r['box'])]; wl = r.get('wL', []); wi = 0
        for ch in r['tree']:
            bx = stack.pop()
            if ch == '+':
                stack.extend(split(bx, S))
            elif '1' <= ch <= '6':
                keep('w', (bx, 12 - 2 * int(ch), int(ch), wl[wi])); wi += 1
            elif 'a' <= ch <= 'l' or 'A' <= ch <= 'L':
                b = ord(ch) - (ord('a') if ch >= 'a' else ord('A'))
                keep('e' if b == 11 else 'r', (bx, b, None, None))
    return res, cnt


ENG = {'bulk': 'a', 'F0': 'f', 'KQ1': 'KQ', 'KQ2': 'KQ', 'KQK1': 'KQ', 'KQK2': 'KQ', 'KQB': 'KQ', 'KA': 'KA'}
for name in names:
    D = load(name)
    if D is None:
        print(name, 'not present'); continue
    eng = ((D.get('header') or {}).get('engine')) or ENG.get(name, name)     # continuation files (NAME.r1) too
    if D.get('header'):
        S = D['header']['scale']; H = D['header']
        SPLIT_RULE = H.get('split_rule')
    else:                                     # partial records: scale as in the drivers
        S = [1.0, 1.0, math.pi / 2, 1.0] if eng in ('f', 'KQ') else ([math.pi / 2, 1.0, math.pi / 2, 1.0] if eng == 'a' else None)
        H = {}
        SPLIT_RULE = None
    res, cnt = leaves(D, S, NLEAF)
    wl, el, rl = res['w'], res['e'], res['r']
    sample = random.sample(wl, min(NLEAF // 3, len(wl))) + random.sample(el, min(NLEAF // 3, len(el)))
    sample += random.sample(rl, min(NLEAF - len(sample), len(rl)))
    viol = 0; odd = 0; npt = 0; nphys = 0; maxN = 0; minmarg = mp.inf; nw = 0
    phi_a = mp.mpf(H.get('phi_a', 0)) if eng == 'KA' else None
    for bx, bound, W, Lw in sample:
        for k in range(NPT):
            x = [0.5 * (bx[2 * i] + bx[2 * i + 1]) if k == 0 else random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(4)]
            if eng not in ('a', 'f', 'KQ', 'KA') and x[0] <= 0:
                x[0] = 0.5 * (bx[0] + bx[1])
            if eng in ('L2P', 'L2Ap', 'L2Am', 'L2Tp', 'L2Tm', 'D2P', 'D2K', 'D2Tp', 'D2Tm') and x[1] <= 0:   # (D2E, D2E2: e = SD nu)
                x[1] = 0.5 * (bx[2] + bx[3])          # E ratio > 0 so that e > 0
            if eng in ('C1p', 'C1pp') and x[3] >= 1:
                x[3] = 0.5 * (bx[6] + bx[7])          # tau < 1 (t = e^2 tau / (1 - tau))
            if eng in ('f', 'KQ', 'KQF'):
                x[0] = min(x[0], 1.0); x[3] = min(max(x[3], 0.0), 1.0)
                z = zeros(coeffs_f(*[mp.mpf(v) for v in x]))
                phys = x[0] < 1 and 0 < x[3] < 1 and x[1] <= 1
            elif eng == 'KA':
                e, l, w, mu = [mp.mpf(v) for v in x]
                e = min(e, 1 - mp.mpf(10) ** -12); t = 1 / (1 + mu * mu)
                z = zeros(coeffs_f(e, l, w, t))
                phys = True
            else:
                e, p, q, t = a_map(eng, x)
                if not (0 < e <= 1):
                    continue
                z = zeros(coeffs_a(e, p, q, t))
                phys = e < 1 and 0 < t < 1 and p * p + q * q <= 1
            on = [r for r in z if r != 0 and abs(mp.log(abs(r))) < mp.mpf(10) ** -30]
            npt += 1; nphys += phys
            if eng == 'KA':
                Na = sum(1 for r in on if abs(mp.pi - abs(mp.arg(r))) >= phi_a - mp.mpf(10) ** -12)
                maxN = max(maxN, Na)
                if Na > bound:
                    viol += 1; print('   KA VIOLATION', bx, x, Na, bound, flush=True)
                continue
            if eng == 'KQF':                          # outer count of the face construction: |f - pi| >= 2 atan(T_X)
                No = sum(1 for r in on if mp.pi - abs(mp.arg(r)) >= 2 * mp.atan(mp.mpf(0.06)) - mp.mpf(10) ** -12)
                maxN = max(maxN, No)
                if No > bound:
                    viol += 1; print('   KQF VIOLATION', bx, x, No, bound, flush=True)
                continue
            Nm = len(on); maxN = max(maxN, Nm)
            if Nm % 2:
                odd += 1
            if Nm > min(bound, 10):
                viol += 1; print('   VIOLATION', name, bx, x, Nm, bound, flush=True)
            if W is not None:
                rho = mp.e ** mp.mpf(Lw); nw += 1
                ins = sum(1 for r in z if abs(r) < rho)
                minmarg = min([minmarg] + [abs(abs(r) - rho) / rho for r in z])
                if ins != W:
                    viol += 1; print('   WINDING VIOLATION', name, bx, x, ins, W, Lw, flush=True)
    print(f"{name} [{eng}]{' (partial)' if D.get('partial') else ''}: leaves {sum(cnt.values())} (winding {cnt['w']}, total-11 {cnt['e']}); "
          f"sampled {len(sample)} leaves, {npt} points ({nphys} physical, {nw} winding); violations {viol}; odd N_mult {odd}; "
          f"max {'away count' if eng == 'KA' else 'N_mult'} {maxN}; min winding margin {mp.nstr(minmarg, 3)}", flush=True)

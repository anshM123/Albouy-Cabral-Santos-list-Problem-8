"""AUDIT-5 e2: D2 runs (D2E, D2P, D2K, D2Tp, D2Tm): own replay, re-certification of EVERY leaf with the right engine
(regular -> nbeng_D2*.certify_box must reproduce the recorded bound; winding -> fixed wind_cert.wind_box with the stored
L must reproduce W), pruned 'u' leaves checked unphysical at sampled points, and 100-digit spot checks with MY OWN
physics: real roots of f_t (C1p function up to a positive factor) <= bound, and #{|z| < e^L} == W for winding leaves.
usage: python e2_d2_recert.py NSPOT SEED RUN.json [...]   (.venv python)"""
import os, sys, json, math, random, importlib, itertools
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.join(HERE, '..', '..', 'code'); sys.path.insert(0, CODE)
import wind_cert
mp.mp.dps = 100
NSPOT = int(sys.argv[1]); random.seed(int(sys.argv[2])); files = sys.argv[3:]
SD = 5e-3


def split(bx, S):
    w = [(bx[2 * i + 1] - bx[2 * i]) / S[i] for i in range(4)]
    i = int(np.argmax(w)); mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)


def d2_phys(eng, x):
    sd = mp.mpf(SD); x = [mp.mpf(v) for v in x]
    if eng in ('D2Tp', 'D2Tm'):
        n, E, P, K = x; sg = 1 if eng == 'D2Tp' else -1
        e, p, k, tp = sd * E * n, sd * P * n, sd * K * n, sg * sd * n
    else:
        n, r1, r2, x3 = x; T = -mp.mpf(1.05) + 2 * mp.mpf(1.05) * x3
        if eng == 'D2E':
            e, p, k = sd * n, sd * r1 * n, sd * r2 * n
        elif eng == 'D2P':
            e, p, k = sd * r1 * n, sd * n, sd * r2 * n
        else:
            e, p, k = sd * r1 * n, sd * r2 * n, sd * n
        tp = sd * T * n
    tau = mp.mpf(1) / 2 + tp
    return e, p, 1 - k, e * e * tau / (1 - tau)


def roots(e, p, q, t):
    b = mp.sqrt(1 - e * e); N = 32; vals = []
    for n in range(N):
        u = 2 * mp.pi * n / N
        r = 1 - e * mp.cos(u); Z = p * (mp.cos(u) - e) + q * b * mp.sin(u); Zp = -p * mp.sin(u) + q * b * mp.cos(u)
        A = 2 * e * mp.sin(u) * r; rho2 = r * r - Z * Z; d = A - 2 * Z * Zp
        vals.append((1 - t) * A * A * rho2 - t * d * d)
    c = [mp.fsum(vals[n] * mp.exp(-1j * k * 2 * mp.pi * n / N) for n in range(N)) / N for k in range(-6, 7)]
    sc = max(abs(v) for v in c)
    lo = 0
    while abs(c[lo]) < mp.mpf(10) ** -70 * sc: lo += 1
    hi = 12
    while abs(c[hi]) < mp.mpf(10) ** -70 * sc: hi -= 1
    return list(mp.polyroots(c[lo:hi + 1][::-1], maxsteps=500, extraprec=400)) + [mp.mpf(0)] * lo


for fn in files:
    if not os.path.exists(fn):
        print(fn, 'not present'); continue
    D = json.load(open(fn)); H = D['header']; eng = H['engine']; S = H['scale']
    m = importlib.import_module('nbeng_' + eng)
    affine = eng in ('D2E', 'D2P', 'D2K')
    leaves = []; prunes = []; bad = 0
    for r in D['results']:
        stack = [tuple(r['box'])]; wi = 0
        for c in r['tree']:
            bx = stack.pop()
            if c == '+':
                stack.extend(split(bx, S))
            elif c == 'u':
                prunes.append(bx)
            elif '1' <= c <= '6':
                leaves.append((bx, 12 - 2 * int(c), 'w', r['wL'][wi], int(c))); wi += 1
            elif 'a' <= c <= 'k':
                leaves.append((bx, ord(c) - ord('a'), 'r', None, None))
            elif 'A' <= c <= 'K':
                leaves.append((bx, ord(c) - ord('A'), 'd', None, None))
            elif c == '!':
                pass
            else:
                bad += 1
        if wi != len(r['wL']) or (r['done'] and stack):
            bad += 1
    rc = dict(leaves=len(leaves), prunes=len(prunes), wind=sum(1 for l in leaves if l[2] == 'w'), reg_mismatch=0,
              wind_mismatch=0, deep_mismatch=0, prune_viol=0, spot_viol=0, replay_bad=bad)
    for bx, bound, kind, L, W in leaves:
        if kind == 'r':
            ok, tot = m.certify_box(*bx, H['J'], H['M0'], H['rounds'])
            if not ok or int(tot) != bound: rc['reg_mismatch'] += 1
        elif kind == 'd':
            ok, tot = m.certify_box_deep(*bx, H['J'], H['M0'], H['deep_rounds'], H['deep_cap'])
            if not ok or int(tot) > bound: rc['deep_mismatch'] += 1
        else:
            cm = m.box_coeffs(*bx)
            ok, W2 = wind_cert.wind_box(*cm, bx[6], bx[7], L, H.get('wind_m0', 64), H.get('wind_rounds', 3), H.get('wind_cap', 1024))
            if not ok or W2 != W: rc['wind_mismatch'] += 1
    for bx in prunes:                          # every sampled point must have l > 1
        for x in list(itertools.product(*[bx[2 * i:2 * i + 2] for i in range(4)]))[:16] + \
                [tuple(random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(4)) for _ in range(8)]:
            e, p, q, t = d2_phys(eng, x)
            if p * p + q * q <= 1: rc['prune_viol'] += 1
    for bx, bound, kind, L, W in random.sample(leaves, min(NSPOT, len(leaves))):
        x = [random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(4)]
        x[0] = max(x[0], 1e-9)
        if eng not in ('D2E',):
            x[1] = max(x[1], 1e-9)               # keep e > 0 (E ratio) for the physical evaluation
        e, p, q, t = d2_phys(eng, x)
        z = roots(e, p, q, t)
        n_on = sum(1 for r in z if r != 0 and abs(mp.log(abs(r))) < mp.mpf(10) ** -25)
        if n_on > bound: rc['spot_viol'] += 1; print('   SPOT VIOLATION', eng, bx, x, n_on, bound)
        if kind == 'w':
            inside = sum(1 for r in z if abs(r) < mp.e ** L)
            if inside != W: rc['spot_viol'] += 1; print('   WIND SPOT VIOLATION', eng, bx, x, inside, W)
    print(os.path.basename(fn), rc, flush=True)

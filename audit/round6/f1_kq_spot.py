"""AUDIT-6 f1: own replay + 100-digit spot checks for e -> 1 runs in f-chart variables (e, l, w, t):
bb_run 'f' JSONs (F0) and kq_run outputs (final JSON or the partial .records.jsonl of completed initial boxes).
At sampled points of sampled leaves: N_mult(P) = number of real zeros, with multiplicity, of the real trigonometric
polynomial P(f) = mu^2 k^2 sin^2 f sigma^2 - D^2 B^2 (R-normalised, mu = sqrt((1-t)/t)), computed from MY formula
at 100 digits (24-point DFT -> degree-12 companion polynomial -> roots on |z| = 1).  Checks N_mult <= bound (11 -> 10
by parity) and that N_mult is even.  (P is a positive multiple of the f-chart polynomial, so this is the stationary-
point polynomial of the configuration.)
usage: python f1_kq_spot.py NLEAF NPT SEED FILE [...]"""
import sys, json, math, random
import numpy as np
import mpmath as mp
mp.mp.dps = 100
NLEAF = int(sys.argv[1]); NPT = int(sys.argv[2]); random.seed(int(sys.argv[3])); files = sys.argv[4:]


def split(bx, S):
    w = [(bx[2 * i + 1] - bx[2 * i]) / S[i] for i in range(4)]
    i = int(np.argmax(w)); mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)


def bound_of(ch):
    if 'a' <= ch <= 'l': return ord(ch) - ord('a')
    if 'A' <= ch <= 'L': return ord(ch) - ord('A')
    if '1' <= ch <= '6': return 12 - 2 * int(ch)
    return None


def n_real_P(e, l, w, t):
    e, l, w, t = [mp.mpf(v) for v in (e, l, w, t)]
    pz, qz = l * mp.sin(w), l * mp.cos(w)
    mu = mp.sqrt((1 - t) / t)
    N = 24; vals = []
    for n in range(N):
        f = 2 * mp.pi * n / N
        zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
        B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta); D = 1 + e * mp.cos(f); k = 2 * e * (1 + e)
        vals.append(mu ** 2 * k ** 2 * mp.sin(f) ** 2 * (1 - zeta ** 2) - D ** 2 * B ** 2)
    c = [mp.fsum(vals[n] * mp.exp(-1j * kk * 2 * mp.pi * n / N) for n in range(N)) / N for kk in range(-6, 7)]
    sc = max(abs(v) for v in c)
    lo = 0
    while abs(c[lo]) < mp.mpf(10) ** -80 * sc: lo += 1
    hi = 12
    while abs(c[hi]) < mp.mpf(10) ** -80 * sc: hi -= 1
    z = mp.polyroots(c[lo:hi + 1][::-1], maxsteps=600, extraprec=500)
    return sum(1 for r in z if abs(mp.log(abs(r))) < mp.mpf(10) ** -30)


for fn in files:
    if fn.endswith('.jsonl'):
        recs = []
        for line in open(fn):
            try:
                r = json.loads(line)
            except Exception:
                continue
            if r.get('done'):
                recs.append(r)
        S = (1.0, 1.0, math.pi / 2, 1.0); eng = 'KQ (partial records)'
    else:
        D = json.load(open(fn)); recs = D['results']; S = D['header']['scale']; eng = D['header'].get('engine')
    leaves = []; bad = 0
    for r in recs:
        stack = [tuple(r['box'])]
        for ch in r['tree']:
            bx = stack.pop()
            if ch == '+':
                stack.extend(split(bx, S))
            elif ch == '!':
                continue
            else:
                b = bound_of(ch)
                if b is None:
                    bad += 1
                else:
                    leaves.append((bx, b))
        if r.get('done') and stack:
            bad += 1
    hist = {}
    for _, b in leaves:
        hist[b] = hist.get(b, 0) + 1
    sample = random.sample(leaves, min(NLEAF, len(leaves)))
    viol = 0; odd = 0; n = 0; mx = 0
    for bx, b in sample:
        for k in range(NPT):
            x = [random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(4)]
            x[0] = min(x[0], 1 - 1e-12); x[3] = min(max(x[3], 1e-9), 1 - 1e-15)
            N = n_real_P(*x)
            n += 1; mx = max(mx, N)
            if N % 2: odd += 1
            if N > min(b, 10):
                viol += 1; print('   VIOLATION', fn, bx, x, N, b)
    print(f"{fn.split('/')[-1]} [{eng}]: records {len(recs)}, leaves {len(leaves)}, bad codes/stack {bad}, bound histogram "
          f"{dict(sorted(hist.items()))}; spot points {n}: violations {viol}, odd N_mult {odd}, max N_mult {mx}", flush=True)

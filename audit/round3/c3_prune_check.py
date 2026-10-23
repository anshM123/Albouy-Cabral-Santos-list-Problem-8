"""AUDIT-3 c3: are the pruned C3 leaves ('u' = l > 1 everywhere, 'o' = outside the cell everywhere) really so?
Samples corners + random points of each pruned leaf AND of the root variable where it enters (Phi: m; L2*Y: mp),
maps them with my own chart maps to (Qz, kappa, eps, mu) and checks the predicate at every sample (nu > 0)."""
import sys, json, math, random, itertools
import numpy as np
def split(bx, S):
    w = [(bx[2 * i + 1] - bx[2 * i]) / S[i] for i in range(4)]
    i = int(np.argmax(w))
    mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)
QZMAX = math.cos(9 * math.pi / 20); KAPMAX = 1 - 0.95 * math.sin(9 * math.pi / 20); EPSMAX = 0.01; MUMAX = 0.0505
Y2, E2, M2 = 0.25, 0.01, 0.01
random.seed(5)
def phys(ch, x, r):
    nu = x[0]
    if ch == 'Ca': return nu, nu**2 * x[1], nu**2 * x[2], nu**3 * x[3]
    if ch == 'CK': return nu * x[1], nu**2, nu**2 * x[2], nu**3 * x[3]
    if ch == 'CE': return nu * x[1], nu**2 * x[2], nu**2, nu**3 * x[3]
    if ch == 'Cm': return nu * x[1], nu**2 * x[2], nu**2 * x[3], nu**3
    if ch in ('PhiP', 'PhiM'): return nu * x[1], nu**2 * x[2], nu**2 * x[3], nu**3 * r
    rv, eta, w4 = x[1], x[2], x[3]
    qz, kap = (nu, nu**2 * rv) if ch.startswith('L2a') else (nu * rv, nu**2)
    k = ch[3:]
    if k == 'E': E, m = E2 * eta**2, M2 * eta**2 * w4
    elif k == 'm': E, m = E2 * eta**2 * w4, M2 * eta**2
    else: E, m = E2 * eta**2 * w4, M2 * eta**2 * r
    return qz, kap, nu**2 * E, nu**3 * m
RROOT = {'PhiP': (0, 1.001**3), 'PhiM': (0, 1.001**3), 'L2aYp': (0, 1.001**2), 'L2aYm': (0, 1.001**2),
         'L2KYp': (0, 1.001**2), 'L2KYm': (0, 1.001**2)}
for fn in sys.argv[1:]:
    D = json.load(open(fn)); H = D['header']; ch = H['chart']; S = H['scale']
    pruned = []
    for r in D['results']:
        stack = [tuple(r['box'])]
        for c in r['tree']:
            bx = stack.pop()
            if c == '+': stack.extend(split(bx, S))
            elif c in 'uo': pruned.append((bx, c))
    bad = {'u': 0, 'o': 0}; n = {'u': 0, 'o': 0}
    for bx, c in pruned:
        pts = [tuple(bx[2*i + b[i]] for i in range(4)) for b in itertools.product((0, 1), repeat=4)]
        pts += [tuple(random.uniform(bx[2*i], bx[2*i+1]) for i in range(4)) for _ in range(40)]
        rlo, rhi = RROOT.get(ch, (0, 0))
        for x in pts:
            if x[0] <= 0: continue
            for r in (rlo, rhi, random.uniform(rlo, rhi)):
                qz, kap, ep, mu = phys(ch, x, r)
                n[c] += 1
                if c == 'u' and not (qz * qz > 2 * kap - kap * kap):          # l <= 1 somewhere -> physical point dropped
                    bad['u'] += 1
                if c == 'o' and not (qz > QZMAX or kap > KAPMAX or ep > EPSMAX or mu > MUMAX):
                    bad['o'] += 1
    print(f"{fn.split('/')[-1]}: pruned leaves {len(pruned)} ('u' {sum(1 for p in pruned if p[1]=='u')}, 'o' {sum(1 for p in pruned if p[1]=='o')}); "
          f"sample checks u {n['u']} (violations {bad['u']}), o {n['o']} (violations {bad['o']})")

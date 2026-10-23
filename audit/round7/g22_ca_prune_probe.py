"""AUDIT-7 g22: locate and examine pruned ('u'/'o') leaves of the C3 Ca base file for which a sampled point is not
unphysical / outside the C3 cell (g2 prune check), using exact (mpmath) arithmetic for the predicate."""
import sys, os, json, math, itertools, random
import mpmath as mp
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import g2_replay_all as g2
mp.mp.dps = 40
D = json.load(open(os.path.join(g2.ROOT, 'c3/runs/Ca.json'))); H = D['header']; S = H['scale']
QZ = mp.cos(mp.mpf('1.413716694115407')); KAP = 1 - mp.mpf('0.95') * mp.sin(mp.mpf('1.413716694115407'))
EPS, MU = mp.mpf('0.01'), mp.mpf('0.0505')
print('header scale', S, 'engine_module', H.get('engine_module'))
bad = []
npr = 0
for r in D['results']:
    stack = [tuple(r['box'])]
    for ch in r['tree']:
        bx = stack.pop()
        if ch == '+':
            w = [(bx[2 * i + 1] - bx[2 * i]) / S[i] for i in range(4)]
            i = max(range(4), key=lambda j: w[j]); mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
            a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
            stack += [tuple(a), tuple(c)]
        elif ch in 'uo':
            npr += 1
            pts = list(itertools.product(*[bx[2 * i:2 * i + 2] for i in range(4)])) + \
                [tuple(random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(4)) for _ in range(30)]
            for x in pts:
                if x[0] <= 0:
                    continue
                nu = mp.mpf(x[0]); qz, kap, eps, mu = nu, nu * nu * mp.mpf(x[1]), nu * nu * mp.mpf(x[2]), nu ** 3 * mp.mpf(x[3])
                l2 = (1 - kap) ** 2 + qz * qz
                ok = (l2 > 1) if ch == 'u' else (qz > QZ or kap > KAP or eps > EPS or mu > MU)
                if not ok:
                    bad.append((ch, bx, x, float(l2 - 1), float(qz - QZ), float(kap - KAP), float(eps - EPS), float(mu - MU)))
                    break
print('pruned leaves', npr, '; leaves with a sampled point violating the predicate', len(bad))
for b in bad[:5]:
    print('  code', b[0], 'box', b[1]); print('     point', b[2], ' l^2-1 %.3e  qz-QZ %.3e  kap-KAP %.3e  eps-EPS %.3e  mu-MU %.3e' % b[3:])

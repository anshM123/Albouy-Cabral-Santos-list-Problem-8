"""AUDIT-8 h5: (1) l*(w) on the whole range l in [0, 1]: the smallest l for which B (e = 1) has 4 real zeros, i.e. the
corner limit e = 1, mu = 0 has N_mult = 12; (2) where the currently open boxes of the KQ chains (KQK1, KQK2, KQB) lie:
how many touch the corner face e = 1 and t = 1, and whether their (l, w) range meets {l >= l*(w)}."""
import json, os, sys, math
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h4_corner_limit import nreal
CODE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code'))


def lstar(w):
    ls = np.linspace(0.0, 1.0, 4001)
    four = [l for l in ls if nreal(l, w) >= 4]
    return min(four) if four else None


W = np.linspace(0.0, math.pi / 2, 41)
LS = [lstar(w) for w in W]
print('w      l*(w)   (corner e = 1, mu = 0 has 12 real zeros for l >= l*)')
for w, l in zip(W, LS):
    print(f'{w:6.3f}  {l if l is None else round(l, 4)}')


def lstar_min(w0, w1):
    v = [l for w, l in zip(W, LS) if w0 - 0.04 <= w <= w1 + 0.04 and l is not None]
    return min(v) if v else 2.0


for name in ('KQK1', 'KQK2', 'KQB'):
    files = [f'{name}.json'] + [f'{name}.r{k}.json' for k in range(1, 10)]
    owed = set(); 
    for k, fn in enumerate(files):
        p = os.path.join(CODE, 'runs', fn)
        if not os.path.exists(p):
            break
        D = json.load(open(p))
        for r in D['results']:
            root = tuple(r['box'])
            owed.discard(root)
            for b in r['left']:
                owed.add(tuple(b))
            for f in r['fails']:
                owed.add(tuple(f[0]))
    corner = [b for b in owed if b[1] >= 1.0 and b[7] >= 1.0]
    deg = [b for b in corner if b[3] >= lstar_min(b[4], b[5])]
    print(f'{name}: open boxes {len(owed)}; touching e = 1 and t = 1: {len(corner)}; of these with l_hi >= l*(w) '
          f'(corner limit has 12 real zeros somewhere in the box): {len(deg)}')
    if owed:
        for i, nm in enumerate('elwt'):
            print(f'   {nm} range of open boxes: [{min(b[2*i] for b in owed):.6g}, {max(b[2*i+1] for b in owed):.6g}]')

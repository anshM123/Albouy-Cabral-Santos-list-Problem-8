"""AUDIT-2 b2: re-certification of recorded leaves with the CORRECT engine (bb_replay.py --recert falls through to
nb_corner for B*/L2* engines), plus brute-force root counts inside each re-certified leaf.

For each run JSON: replay the tree (own replay, same split rule), collect the leaves, then for a sample of leaves
(or all if N < 0):
  * regular leaf ('a'..'k'): nbeng_<chart>.certify_box(box, J, M0, rounds) must return (True, recorded bound) exactly;
    certify_box_deep(box, J, M0, deep_rounds, deep_cap) must return True with bound <= recorded;
  * deep leaf ('A'..'K'): certify_box_deep must reproduce (True, bound);
  * brute force: 16 points of the leaf (corners mixed with random) -> real-root count of the chart function
    (coef_val at the atoms; near-circle roots |log|z|| < 1e-3 counted as real) must be <= recorded bound.
usage: python b2_recert.py N SEED RUN1.json [RUN2.json ...]   (run with the .venv python)
"""
import os
import sys
import json
import math
import random
import importlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.join(HERE, '..', '..', 'code')
sys.path.insert(0, CODE)

N = int(sys.argv[1]); random.seed(int(sys.argv[2])); files = sys.argv[3:]
S0 = 0.01


def split(bx, SCALE):
    w = [(bx[2 * i + 1] - bx[2 * i]) / SCALE[i] for i in range(4)]
    i = int(np.argmax(w))
    mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx)
    a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)


def atoms(ch, names, x):
    """atom values at a box point (own implementation)"""
    v = {}
    if ch == 'B1':
        v.update(e=x[0], ph=x[1], qh=x[2]); e = x[0]
    elif ch == 'B2':
        v.update(th=x[0], eh=x[1], sw=math.sin(x[2]), cw=math.cos(x[2])); e = x[0] * x[1]
    elif ch in ('B3p', 'B3m'):
        v.update(eta=x[0], eh=x[1], thh=x[2], sw=math.sin(x[3]), cw=math.cos(x[3])); e = x[0] * x[1]
    elif ch == 'L2E':
        v.update(nu=x[0], P=x[1], A=x[2]); e = S0 * x[0]
    elif ch == 'L2P':
        v.update(nu=x[0], E=x[1], A=x[2]); e = S0 * x[1] * x[0]
    elif ch in ('L2Ap', 'L2Am'):
        v.update(nu=x[0], E=x[1], P=x[2]); e = S0 * x[1] * x[0]
    else:
        v.update(nu=x[0], E=x[1], P=x[2], A=x[3]); e = S0 * x[1] * x[0]
    b = math.sqrt(1 - e * e); v.update(b=b, ib=1 / b)
    return np.array([v[n] for n in names])


def count_real(a, s):
    c = np.zeros(13, complex); c[6] = a[0]
    for k in range(1, 7):
        c[6 + k] = (a[k] - 1j * s[k]) / 2; c[6 - k] = (a[k] + 1j * s[k]) / 2
    mx = np.abs(c).max()
    if mx == 0:
        return 99
    top = 12
    while top > 0 and abs(c[top]) < 1e-15 * mx:
        top -= 1
    low = 0
    while low < top and abs(c[low]) < 1e-15 * mx:
        low += 1
    z = np.roots(c[low:top + 1][::-1])
    return int(np.sum(np.abs(np.log(np.abs(z))) < 1e-3))


for fn in files:
    D = json.load(open(fn)); H = D['header']; ch = H['engine']; SCALE = H['scale']
    eng = importlib.import_module('nbeng_' + ch)
    cm = importlib.import_module('nbcoef_' + ch)
    names = cm.ATOMS; nvar = cm.NVAR; affine = (cm.NF == 2)
    leaves = []
    for r in D['results']:
        stack = [tuple(r['box'])]
        for c in r['tree']:
            bx = stack.pop()
            if c == '+':
                stack.extend(split(bx, SCALE))
            elif c != '!':
                deep = c.isupper()
                leaves.append((bx, ord(c) - (ord('A') if deep else ord('a')), deep))
    sample = leaves if N < 0 else random.sample(leaves, min(N, len(leaves)))
    st = dict(leaves=len(leaves), checked=0, reg_mismatch=0, deep_worse=0, brute_viol=0, deep_leaves=sum(d for _, _, d in leaves))
    ex = None
    for bx, bound, deep in sample:
        st['checked'] += 1
        okd, td = eng.certify_box_deep(*bx, H['J'], H['M0'], H['deep_rounds'], H['deep_cap'])
        if not okd or int(td) > bound:
            st['deep_worse'] += 1; ex = ('deep', bx, bound, okd, td)
        if not deep:
            ok, tot = eng.certify_box(*bx, H['J'], H['M0'], H['rounds'])
            if (not ok) or int(tot) != bound:
                st['reg_mismatch'] += 1; ex = ('reg', bx, bound, ok, tot)
        # brute force
        for k in range(16):
            x = [bx[2 * i] + (bx[2 * i + 1] - bx[2 * i]) * (random.random() if k >= 8 else ((k >> i) & 1 if i < 3 else random.random()))
                 for i in range(4)]
            out = np.zeros((nvar + 1, 2, 2, 7))
            cm.coef_val(atoms(ch, names, x), out)
            if affine:
                a = (1 - x[3]) * out[0, 0, 0] - x[3] * out[0, 1, 0]; s = (1 - x[3]) * out[0, 0, 1] - x[3] * out[0, 1, 1]
            else:
                a = out[0, 0, 0]; s = out[0, 0, 1]
            n = count_real(a, s)
            if n > bound and n != 99:
                st['brute_viol'] += 1; ex = ('brute', bx, x, n, bound)
    print(os.path.basename(fn), st, ('EXAMPLE ' + repr(ex)) if ex else '', flush=True)

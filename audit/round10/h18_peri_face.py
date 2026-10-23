"""AUDIT-10 h18: face mode of --peri (certify_outer -> count_peri(..., PHI_0, PHI_P, t_lo = T_X)): on random boxes of the
KQF domain (e in [0.9996, 1), t in [0.99999983, 1), all l, w outside the C3 cell) the certified outer total must bound
the zeros of P with |th| >= T_X = 0.06 at 100 digits (MY OWN P), at 4 points per box.
usage: python h18_peri_face.py NBOX SEED"""
import os, sys, math, random
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.normpath(os.path.join(HERE, '..', '..', 'code'))
sys.path.insert(0, CODE); sys.path.insert(0, os.path.join(CODE, 'c3')); sys.path.insert(0, HERE); os.chdir(CODE)
import kq_run
from h17_peri_parts import zeros_f
NBOX, SEED = int(sys.argv[1]), int(sys.argv[2]); random.seed(SEED)
kq_run.init(face=True, peri=kq_run.PHI_P)
mp.mp.dps = 100
W9 = 1.413716694115407
nb = cert = viol = 0; mx = 0
while nb < NBOX:
    wd = [10 ** random.uniform(-6, -4.5), 10 ** random.uniform(-3, -1.5), 10 ** random.uniform(-3, -1.5), 10 ** random.uniform(-9, -7.5)]
    e0 = random.uniform(0.9996, 1 - wd[0]); l0 = random.uniform(0, 1 - wd[1]); w0 = random.uniform(0, math.pi / 2 - wd[2])
    if l0 >= 0.95 and w0 >= W9:
        continue
    t0 = random.uniform(0.99999983, 1 - wd[3])
    bx = (e0, e0 + wd[0], l0, l0 + wd[1], w0, w0 + wd[2], t0, t0 + wd[3])
    ok, tot = kq_run.count_peri(bx, 16, 8192, kq_run.PHI_0, kq_run.PHI_P, kq_run.T_X)
    nb += 1
    if not ok:
        continue
    cert += 1
    for k in range(4):
        x = [0.5 * (bx[2 * i] + bx[2 * i + 1]) if k == 0 else random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(4)]
        fz = zeros_f(*x)
        outer = sum(1 for f in fz if abs(mp.tan((f - mp.pi) / 2)) >= mp.mpf(kq_run.T_X))
        mx = max(mx, outer)
        if outer > tot:
            viol += 1; print('   VIOLATION', bx, x, outer, tot, flush=True)
print(f'face mode (--peri): {nb} random KQF boxes, {cert} certified (total <= 6); 100-digit outer counts at 4 points each: '
      f'violations {viol}, max outer zeros {mx}')

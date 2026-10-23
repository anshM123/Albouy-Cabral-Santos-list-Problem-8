"""AUDIT-8 h8: kq_run.certify (current code, re-certification mode) on boxes that contain the corner point
e = 1, t = 1 at (l, w) with l >= l*(w), for shrinking widths, compared with boxes of the same size that stop short of
e = 1 or of t = 1.  Expected: corner boxes never certify (total >= 12 for every sound count, since P has 12 real zeros
at the corner), the others do once small enough."""
import os, sys, math
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.normpath(os.path.join(HERE, '..', '..', 'code'))
sys.path.insert(0, CODE); sys.path.insert(0, os.path.join(CODE, 'c3')); os.chdir(CODE)
import kq_run
kq_run.init(); kq_run.DEEP_ALWAYS = True
for (l, w) in ((1.0, 0.2), (0.999, 0.3), (0.97, 0.6), (0.9, 0.1)):
    for h in (1e-3, 1e-5, 1e-7):
        corner = (1 - h, 1.0, l - h, l, w, w + h, 1 - h * 1e-3, 1.0)
        off_e = (1 - 2 * h, 1 - h, l - h, l, w, w + h, 1 - h * 1e-3, 1.0)
        off_t = (1 - h, 1.0, l - h, l, w, w + h, 1 - 2 * h * 1e-3, 1 - h * 1e-3)
        r = [kq_run.certify(b) for b in (corner, off_e, off_t)]
        print(f'l = {l}, w = {w}, width {h:g}: corner box {r[0]}; box below e = 1: {r[1]}; box below t = 1: {r[2]}', flush=True)

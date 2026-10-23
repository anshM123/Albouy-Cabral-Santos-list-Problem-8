"""AUDIT-3 c9: R2-2 re-check with the current bb_run fallback rule
   deep if not good and (tot < 999 or maxw < 1e-3) and maxw < DEEP_MAXW (1e-2)."""
import os, sys, math
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(HERE, '..', '..', 'code'))
import nb_cert
SCALE = (math.pi / 2, 1.0, math.pi / 2, 1.0)
for e in (0.9, 0.93, 0.95):
    c = (math.asin(e), math.sqrt(e), math.pi / 2, 1.0)
    for h in (1e-3, 1e-4, 1e-5, 2e-6):
        bx = (c[0] - h, c[0] + h, c[1] - h, c[1] + h, c[2] - h, c[2], c[3] - h, c[3])
        ok, tot = nb_cert.certify_box(0, *bx, 6, 64, 8)
        maxw = max((bx[2*i+1] - bx[2*i]) / SCALE[i] for i in range(4))
        calls = (not ok) and (tot < 999 or maxw < 1e-3) and maxw < 1e-2
        res = nb_cert.certify_box_deep(0, *bx, 6, 64, 16, 8192) if calls else None
        print(f"e={e} h={h:.0e} maxw={maxw:.1e}: regular ({ok}, {tot}); deep called: {calls}; deep result: {res}")

"""AUDIT a6c: can certify_box (J=6) ever certify boxes touching the 6-fold curve (w = pi/2, e = l^2, t = 1)
or the near-triple cluster?  Shrink boxes around a point and report the certified total vs half-width."""
import os, sys, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, 'numba_cache'))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(HERE, '..', 'code')); sys.path.insert(0, HERE)
import nb_cert
from a6c_helpers import port_trace

def boxes_at(c, h, faces):
    x0, l, w, t = c
    lo = [x0 - h, l - h, w - h, t - h]; hi = [x0 + h, l + h, w + h, t + h]
    for i, side in faces:
        if side == 'hi': hi[i] = c[i]
        else: lo[i] = c[i]
    return (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], lo[3], hi[3])

for e in (0.3, 0.6, 0.9):
    c = (math.asin(e), math.sqrt(e), math.pi / 2, 1.0)
    row = []
    for h in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7):
        bx = boxes_at(c, h, [(2, 'hi'), (3, 'hi')])
        ok, tot = nb_cert.certify_box(0, *bx, 6, 64, 8)
        row.append(f"h={h:.0e}:{'OK' if ok else 'fail'}({tot})")
    print(f"curve6 e={e}: " + "  ".join(row))
    # also a box away from t = 1 but on the curve in (e, l, w)
    row = []
    for tt in (0.999, 0.99, 0.9):
        bx = boxes_at((c[0], c[1], c[2], tt), 1e-4, [(2, 'hi')])
        ok, tot = nb_cert.certify_box(0, *bx, 6, 64, 8)
        row.append(f"t={tt}:{'OK' if ok else 'fail'}({tot})")
    print(f"      same (e,l,w), t below 1, h=1e-4: " + "  ".join(row))
c2 = (0.451082, 0.989772, 0.261529, 0.196923)
row = []
for h in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7):
    ok, tot = nb_cert.certify_box(0, *boxes_at(c2, h, []), 6, 64, 8)
    row.append(f"h={h:.0e}:{'OK' if ok else 'fail'}({tot})")
print("cluster2 point: " + "  ".join(row))
port_trace(0, boxes_at((math.asin(0.6), math.sqrt(0.6), math.pi/2, 1.0), 1e-5, [(2, 'hi'), (3, 'hi')]))

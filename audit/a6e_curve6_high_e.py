"""AUDIT a6e: required box size near the 6-fold curve (w = pi/2, e = l^2, t = 1) at high e, compared with the
nb_run stopping width minw = 2e-5 (relative to scale pi/2 in eta and w, 1 in l and t).  A box of half-width h
in every coordinate has min relative width 2h/(pi/2) in eta/w; nb_run stops splitting below 2e-5."""
import os, sys, math
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, 'numba_cache'))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(HERE, '..', 'code'))
import nb_cert
for e in (0.8, 0.9, 0.93, 0.95):
    c = (math.asin(e), math.sqrt(e), math.pi / 2, 1.0)
    row = []
    for h in (1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 1e-7, 1e-8):
        bx = (c[0] - h, min(c[0] + h, math.asin(0.95) if e == 0.95 else c[0] + h), c[1] - h, c[1] + h, c[2] - h, c[2], c[3] - h, c[3])
        ok, tot = nb_cert.certify_box(0, *bx, 6, 64, 8)
        row.append(f"{h:.0e}:{'OK' if ok else 'x'}({tot})")
    print(f"e={e}: " + " ".join(row))
print("nb_run: a box is recorded as FAILED once 2h/(pi/2) < 2e-5, i.e. h < %.1e (eta, w) or 2h < 2e-5, h < 1e-5 (l, t)" % (1e-5 * math.pi / 2))

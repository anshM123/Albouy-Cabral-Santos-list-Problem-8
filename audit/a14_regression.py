"""AUDIT a14: effect of the 23:36 change to nb_cert.certify_box (stall rule + 4096-arc cap) on certification,
compared with the old refinement loop (port).  Soundness is unaffected (early exits only return False)."""
import os, sys, math, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, 'numba_cache'))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(HERE, '..', 'code')); sys.path.insert(0, HERE)
import nb_cert
from oldcert_port import port_certify
rng = np.random.default_rng(14)
ETA = (math.asin(0.2), math.asin(0.95))
cnt = dict(both=0, old_only=0, new_only=0, neither=0)
ex = []
t0 = time.time()
for it in range(int(sys.argv[1]) if len(sys.argv) > 1 else 400):
    c = [rng.uniform(*ETA), rng.uniform(0, 1), rng.uniform(0, np.pi / 2), rng.uniform(0, 1)]
    k = it % 4
    if k == 1: c[2] = np.pi / 2
    if k == 2: c[1] = 1.0
    if k == 3: e = rng.uniform(0.2, 0.9); c = [math.asin(e), math.sqrt(e), np.pi / 2, 1.0]
    r = np.array([0.01, 0.01, 0.015, 0.01]) * np.exp(rng.uniform(-5, 0, 4))
    bx = (max(c[0] - r[0], ETA[0]), min(c[0] + r[0], ETA[1]), max(c[1] - r[1], 0), min(c[1] + r[1], 1),
          max(c[2] - r[2], 0), min(c[2] + r[2], np.pi / 2), max(c[3] - r[3], 0), min(c[3] + r[3], 1))
    new_ok, new_t = nb_cert.certify_box(0, *bx, 6, 64, 8)
    old_ok, old_t, *_ = port_certify(0, bx)
    key = 'both' if new_ok and old_ok else 'old_only' if old_ok else 'new_only' if new_ok else 'neither'
    cnt[key] += 1
    if key == 'old_only' and len(ex) < 5:
        ex.append((['%.6g' % v for v in bx], int(old_t), int(new_t)))
print(cnt, f"{time.time()-t0:.0f}s")
for x in ex:
    print("  certified by OLD loop only:", x)

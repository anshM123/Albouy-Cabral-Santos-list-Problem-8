"""AUDIT a6d: empirical soundness of the python v5 engine certify5_ce.certify_box (used for costmap.json):
random a-chart boxes; when certified with bound B, sampled true counts (upward-biased) must be <= B."""
import os, sys, itertools, time
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'code'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import certify5_ce
from a6b_faces import count_real
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 6)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 120
st = dict(boxes=0, cert=0, viol=0); hist = {}; t0 = time.time()
for it in range(N):
    c = (rng.uniform(0.0502, 1.2532), rng.uniform(0, 1), rng.uniform(0, np.pi/2), rng.uniform(0, 1))
    r = np.array([0.01, 0.006, 0.01, 0.006]) * np.exp(rng.uniform(-2.5, 0, 4))
    if it % 4 == 0: c = (c[0], c[1], np.pi/2 - r[2], c[3])      # touch w = pi/2
    if it % 4 == 1: c = (c[0], 1 - r[1], c[2], c[3])            # touch l = 1
    bx = (c[0]-r[0], c[0]+r[0], max(c[1]-r[1], 0), min(c[1]+r[1], 1), max(c[2]-r[2], 0), min(c[2]+r[2], np.pi/2), max(c[3]-r[3], 0), min(c[3]+r[3], 1))
    ok, B = certify5_ce.certify_box(bx)
    st['boxes'] += 1
    if not ok: continue
    st['cert'] += 1; hist[B] = hist.get(B, 0) + 1
    pts = list(itertools.product(*[bx[2*i:2*i+2] for i in range(4)])) + [tuple(rng.uniform(bx[2*i], bx[2*i+1]) for i in range(4)) for _ in range(60)]
    for (x0, l, w, t) in pts:
        n = count_real(0, x0, l, w, t)
        if n > B and n != 99:
            st['viol'] += 1; print("VIOLATION", bx, (x0, l, w, t), n, B)
print(st, "bounds", dict(sorted(hist.items())), f"{time.time()-t0:.0f}s")

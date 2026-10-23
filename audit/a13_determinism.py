"""AUDIT a13: is nb_cert.certify_box deterministic?  Same boxes, repeated calls in one process; prints a digest
so that several processes can be compared."""
import os, sys, math, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, 'numba_cache'))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(HERE, '..', 'code'))
import nb_cert
boxes = []
for e in (0.3, 0.6, 0.9):
    c = (math.asin(e), math.sqrt(e), math.pi / 2, 1.0)
    for h in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
        boxes.append((c[0] - h, c[0] + h, c[1] - h, c[1] + h, c[2] - h, c[2], c[3] - h, c[3]))
c2 = (0.451082, 0.989772, 0.261529, 0.196923)
for h in (1e-3, 1e-4, 1e-5):
    boxes.append(tuple(v for i in range(4) for v in (c2[i] - h, c2[i] + h)))
res = []
for rep in range(int(sys.argv[1]) if len(sys.argv) > 1 else 5):
    res.append(tuple(nb_cert.certify_box(0, *bx, 6, 64, 8) for bx in boxes))
print("within-process repeats identical:", all(r == res[0] for r in res))
print("results:", [(int(ok), int(t)) for ok, t in res[0]])
# lower-level: signs of the first round for one box, repeated
bx = boxes[8]
outs = []
for rep in range(5):
    cmid, crad, dcen, dwid, cmax, rr = nb_cert.box_coeffs(0, *bx[:6])
    um = np.arange(64) * (2 * math.pi / 64) + math.pi / 64
    sg, hop = nb_cert.arc_signs(cmid, crad, dcen, dwid, cmax, rr, bx[6], bx[7], um, np.full(64, math.pi / 64), 6)
    outs.append((hashlib.md5(np.ascontiguousarray(cmid).tobytes() + np.ascontiguousarray(dwid).tobytes()).hexdigest()[:10],
                 hashlib.md5(sg.tobytes()).hexdigest()[:10], bool(hop)))
print("box_coeffs/arc_signs digests:", outs)
print("digest of all results:", hashlib.md5(repr(res[0]).encode()).hexdigest()[:12])

"""AUDIT a6b: certify_box on boxes touching the degenerate faces (exact multiple roots) and on the two failure
clusters of the v5 cost map.  Sampling INCLUDES points exactly on the face; the root count counts near-circle
roots (|log|z|| < 1e-3) as real, so multiple roots are counted with multiplicity (upward bias).
Faces: w = pi/2, w = 0 (P_z Q_z = 0: double roots at u = 0, pi), l = 0 (planar), l = 1 (rho = 0 twice),
t = 0, t = 1 (perfect squares), curve w = pi/2, e = l^2, t -> 1 (6-fold root), cluster-2 (near-triple roots).
usage: python a6b_faces.py NPERFACE SEED
"""
import os
import sys
import math
import time
import itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, 'numba_cache'))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(HERE, '..', 'code'))
sys.path.insert(0, HERE)
import nb_cert
from nbport import phipsi_coeffs

NP = int(sys.argv[1]) if len(sys.argv) > 1 else 40
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 2)
J = 6
ETA = (math.asin(0.2), math.asin(0.95))


def count_real(chart, x0, l, w, t):
    (aP, sP), (aS, sS) = phipsi_coeffs(chart, x0, l, w)
    a = (1 - t) * aP - t * aS; s = (1 - t) * sP - t * sS
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


def clip_box(c, r):
    x0, l, w, t = c
    return (max(x0 - r[0], ETA[0]), min(x0 + r[0], ETA[1]), max(l - r[1], 0.0), min(l + r[1], 1.0),
            max(w - r[2], 0.0), min(w + r[2], np.pi / 2), max(t - r[3], 0.0), min(t + r[3], 1.0))


def face_centre(face):
    c = [rng.uniform(*ETA), rng.uniform(0, 1), rng.uniform(0, np.pi / 2), rng.uniform(0, 1)]
    if face == 'w=pi/2': c[2] = np.pi / 2
    if face == 'w=0': c[2] = 0.0
    if face == 'l=0': c[1] = 0.0
    if face == 'l=1': c[1] = 1.0
    if face == 't=0': c[3] = 0.0
    if face == 't=1': c[3] = 1.0
    if face == 'curve6':
        e = rng.uniform(0.2, 0.9); c = [math.asin(e), math.sqrt(e), np.pi / 2, 1.0]
    if face == 'cluster2':
        c = [rng.uniform(0.4501, 0.4511), rng.uniform(0.9891, 0.9900), rng.uniform(0.2598, 0.2618), rng.uniform(0.19688, 0.19785)]
    return c


faces = ['w=pi/2', 'w=0', 'l=0', 'l=1', 't=0', 't=1', 'curve6', 'cluster2']
t0 = time.time()
for face in faces:
    st = dict(boxes=0, cert=0, viol=0, hist={})
    ex = None
    for it in range(NP):
        c = face_centre(face)
        r = np.array([0.01, 0.01, 0.015, 0.01]) * np.exp(rng.uniform(-4, 0, 4))
        box = clip_box(c, r)
        ok, B = nb_cert.certify_box(0, *box, J, 64, 8)
        st['boxes'] += 1
        if not ok:
            continue
        st['cert'] += 1
        st['hist'][int(B)] = st['hist'].get(int(B), 0) + 1
        pts = list(itertools.product(box[0:2], box[2:4], box[4:6], box[6:8]))
        pts += [tuple(c)] if all(box[2*i] <= c[i] <= box[2*i+1] for i in range(4)) else []
        for _ in range(80):
            pts.append(tuple(rng.uniform(box[2*i], box[2*i+1]) for i in range(4)))
        for (x0, l, w, t) in pts:
            n = count_real(0, x0, l, w, t)
            if n > B and n != 99:
                st['viol'] += 1
                ex = (box, (x0, l, w, t), n, int(B))
    print(f"[{face:9s}] boxes {st['boxes']}, certified {st['cert']}, violations {st['viol']}, bounds {dict(sorted(st['hist'].items()))}"
          + (f"  EXAMPLE {ex}" if ex else ""), flush=True)
print(f"{time.time()-t0:.0f}s")

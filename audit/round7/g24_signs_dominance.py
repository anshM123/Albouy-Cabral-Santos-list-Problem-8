"""AUDIT-7 g24: after the 09:41 patch of nb_cert2 (radius = min of two complete radii), every sign certified by the
first-order nb_cert.arc_signs must be certified, with the same sign, by nb_cert2.arc_signs2 (charts 0 and 1), and every
half-plane label of wind_cert.labels must be reproduced by nb_cert2.labels2 (same arcs, random boxes, t-ranges, L).
usage: python g24_signs_dominance.py NBOX SEED"""
import os, sys, math, random
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.normpath(os.path.join(HERE, '..', '..', 'code')); sys.path.insert(0, CODE)
import nb_cert, nb_cert2, wind_cert
NBOX, SEED = int(sys.argv[1]), int(sys.argv[2]); random.seed(SEED)
ns = nlost = ngain = nl = nllost = 0
for it in range(NBOX):
    ch = it % 2
    wd = [10 ** random.uniform(-5, -0.8) for _ in range(3)]
    if ch == 0:
        x0 = random.uniform(0.2, 1.25 - wd[0])
    else:
        x0 = random.uniform(0.9, 1.0 - wd[0])
    l0 = random.uniform(0, 1 - wd[1]); w0 = random.uniform(0, math.pi / 2 - wd[2])
    bx = (x0, x0 + wd[0], l0, l0 + wd[1], w0, w0 + wd[2])
    t0 = random.uniform(0, 1); t1 = min(1.0, t0 + 10 ** random.uniform(-6, -0.5))
    M = 64
    um = (np.arange(M) + 0.5) * (2 * math.pi / M); hh = np.full(M, math.pi / M + nb_cert.ARC_PAD)
    d1 = nb_cert.box_coeffs(ch, *bx); d2 = nb_cert2.box_data2(ch, *bx)
    s1, _ = nb_cert.arc_signs(*d1, t0, t1, um, hh, 6)
    s2, _ = nb_cert2.arc_signs2(*d2, t0, t1, um, hh, 6)
    m1 = s1 != 0
    ns += int(m1.sum()); nlost += int((s2[m1] != s1[m1]).sum()); ngain += int(((s1 == 0) & (s2 != 0)).sum())
    L = -random.uniform(0.05, 3.0)
    l1 = wind_cert.labels(*d1, t0, t1, L, um, hh)
    l2 = nb_cert2.labels2(*d2, t0, t1, L, um, hh)
    lab1 = l1 != 0
    nl += int(lab1.sum())
    # a label certified first-order must still be certified (labels2 may pick an earlier half-plane in its own order,
    # which is equally valid; count only arcs that lose every label)
    nllost += int((l2[lab1] == 0).sum())
print(f'arc_signs: {ns} certified (arc, j) signs of nb_cert over {NBOX} boxes; lost or changed in nb_cert2 {nlost}; '
      f'additionally certified by nb_cert2 {ngain}')
print(f'labels: {nl} labelled arcs of wind_cert; unlabelled by nb_cert2.labels2 {nllost}')

"""AUDIT-11 i4: why the KQF2 re-certification must use the header T_X.  A small KQ box at the authors' example
point (sigma = 0.01877, E = 0.25, Pz = 0.9863, Qz = 0.1646, m = 0.9975) has outer count 6 from |tau| >= 0.06 but 8 from
|tau| >= 0.030015 (zeros at tau = -0.041 and +0.056).  kq_run.certify_outer with the default T_X = 0.06 (what
bb_replay.py uses for every KQF-engine file) certifies it with bound <= 6; with the KQF2 header T_X it must not."""
import os, sys, math
CODE = r'<repo>\code'
sys.path.insert(0, CODE); sys.path.insert(0, os.path.join(CODE, 'c3')); os.chdir(CODE)
import kq_run
sg, E, pz, qz, m = 0.01877, 0.25, 0.9863, 0.1646, 0.9975
e = 1 - sg * sg * E; mu = sg * sg * m; t = 1 / (1 + mu * mu); l = math.hypot(pz, qz); w = math.atan2(pz, qz)
bx = (e - 1e-9, e + 1e-9, l - 1e-7, l + 1e-7, w - 1e-7, w + 1e-7, t - 1e-12, t + 1e-12)
for tx in (0.06, 0.030015):
    kq_run.init([0.003, 1.0, math.pi / 2, 3e-7], True, faces=tx / 3, tx=tx)
    kq_run.DEEP_ALWAYS = True
    print('T_X =', tx, ': certify_outer ->', kq_run.certify_outer(bx), flush=True)

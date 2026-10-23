"""AUDIT-2 b8: direct soundness test of the 4-variable `signs` (MV4) enclosures and of the affine-chart enclosures in
the blow-up engines: for random boxes, every arc/derivative sign claimed by nbeng_<chart>.signs is checked at sampled
(x in box, u in padded arc) points, with the chart function evaluated from coef_val (exact up to rounding).
usage: python b8_mv4_signs.py NBOX SEED   (.venv python)"""
import os, sys, math, importlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(HERE, '..', '..', 'code')); sys.path.insert(0, HERE)
from b2_recert import atoms   # own atom map (b2_recert has no side effects at import except argv parsing)
NB = int(sys.argv[1]) if len(sys.argv) > 1 else 40
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 8)
J = 6
DOMS = {'B3p': [(0, 0.13), (0, 0.5654), (0, 1.6154), (0, math.pi / 2)],
        'L2Tp': [(0, 1.01), (0, 1.05), (0, 1.05), (-1.05, 1.05)],
        'B2': [(0, 0.2), (0, 0.3675), (0, math.pi / 2), (0, 1)],
        'L2P': [(0, 1.01), (0, 1.05), (-1.05, 1.05), (0, 1)]}


def derivs(a, s, u, jmax):
    ks = np.arange(7.0)
    out = np.zeros((jmax + 1, u.size))
    for j in range(jmax + 1):
        ang = np.outer(u, ks) + j * np.pi / 2
        out[j] = (np.cos(ang) * ks ** j) @ a + (np.sin(ang) * ks ** j) @ s
    return out, np.sum(np.abs(a)) + np.sum(np.abs(s))


for ch, dom in DOMS.items():
    eng = importlib.import_module('nbeng_' + ch); cm = importlib.import_module('nbcoef_' + ch)
    names = cm.ATOMS; nvar = cm.NVAR; affine = cm.NF == 2
    nsign = 0; nchk = 0; viol = 0
    for ib in range(NB):
        c = [rng.uniform(*d) for d in dom]
        r = np.array([d[1] - d[0] for d in dom]) * 0.03 * np.exp(rng.uniform(-3, 0, 4))
        bx = []
        for i in range(4):
            bx += [max(c[i] - r[i], dom[i][0]), min(c[i] + r[i], dom[i][1])]
        cmid, crad, dcen, dwid, cmax, rr = eng.box_coeffs(*bx)
        M = 64
        left = np.arange(M) * (2 * math.pi / M); width = np.full(M, 2 * math.pi / M)
        um = left + width / 2; hh = width / 2 + eng.ARC_PAD
        sg, _ = eng.signs(cmid, crad, dcen, dwid, cmax, rr, bx[6], bx[7], um, hh, J)
        nsign += int(np.sum(sg != 0))
        for k in range(30):
            x = [bx[2 * i] + (bx[2 * i + 1] - bx[2 * i]) * (rng.uniform() if k >= 16 else ((k >> i) & 1)) for i in range(4)]
            out = np.zeros((nvar + 1, 2, 2, 7))
            cm.coef_val(atoms(ch, names, x), out)
            if affine:
                a = (1 - x[3]) * out[0, 0, 0] - x[3] * out[0, 1, 0]; s = (1 - x[3]) * out[0, 0, 1] - x[3] * out[0, 1, 1]
            else:
                a = out[0, 0, 0]; s = out[0, 0, 1]
            uu = np.concatenate([um - hh, um + hh, um + rng.uniform(-1, 1, M) * hh]); idx = np.concatenate([np.arange(M)] * 3)
            vals, mag = derivs(a, s, uu, J)
            noise = 1e-12 * (mag * 6.0 ** np.arange(J + 1))[:, None]
            S = sg[:, idx]
            bad = ((S == 1) & (vals < -noise)) | ((S == -1) & (vals > noise))
            nchk += int(np.sum(S != 0)); viol += int(bad.sum())
    print(f"{ch:5s}: {NB} boxes, {nsign} signed (arc, j) pairs, {nchk} point checks, violations {viol}", flush=True)

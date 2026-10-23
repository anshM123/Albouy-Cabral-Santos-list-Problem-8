"""AUDIT a5: soundness of the nb_cert enclosures, tested against an independent evaluator.

Run with the .venv python (numba).  Set NUMBA_CACHE_DIR to audit/numba_cache so nothing is written to code/.
  A. the numpy port (nbport.py) reproduces the numba signs exactly on random boxes (both charts);
  B. pointwise parameter enclosure:  f_t^{(j)}(u_m; x) in [vc - rad, vc + rad]  for x in the box, t = t0, t1;
     reports max |f - vc| / rad over box corners, edge points and random points (must be <= 1);
  C. full arc x box enclosure: f_t^{(j)}(u; x) in [lo, hi] for u in the arc, x in the box, t in [t0, t1];
  D. sign claims of the numba arc_signs: every arc with sg[j,i] = +-1 is checked at sampled (u, x, t).
usage: python a5_enclosures.py NBOX SEED
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
from nbcoef_a import coef_val as cval_a, coef_abs as cabs_a
from nbcoef_f import coef_val as cval_f, coef_abs as cabs_f
from nbport import box_coeffs_port, arc_bounds_port, signs_port, ft_derivs

NBOX = int(sys.argv[1]) if len(sys.argv) > 1 else 60
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
J = 6


def random_box(chart, scale):
    if chart == 0:
        x0c = rng.uniform(0.2014, 1.2532); lc = rng.uniform(0, 1); wc = rng.uniform(0, np.pi / 2); tc = rng.uniform(0, 1)
    else:
        x0c = rng.uniform(0.9, 1.0); lc = rng.uniform(0, 1); wc = rng.uniform(0, np.pi / 2); tc = rng.uniform(0, 0.9)
    r = scale * np.exp(rng.uniform(-1.5, 0, 4))
    box = (x0c - r[0], x0c + r[0], max(lc - r[1], 0.0), min(lc + r[1], 1.0), max(wc - r[2], 0.0),
           min(wc + r[2], np.pi / 2), max(tc - r[3], 0.0), min(tc + r[3], 1.0 if chart == 0 else 0.9))
    if chart == 1:
        box = (max(box[0], 0.9), min(box[1], 1.0)) + box[2:]
    return box


def sample_params(box, nrand=40):
    x0lo, x0hi, llo, lhi, wlo, whi = box[:6]
    pts = [c for c in itertools.product((x0lo, x0hi), (llo, lhi), (wlo, whi))]
    for i in range(3):            # edge midpoints / face centres
        for side in (0, 1):
            m = [0.5 * (x0lo + x0hi), 0.5 * (llo + lhi), 0.5 * (wlo + whi)]
            m[i] = box[2 * i + side]
            pts.append(tuple(m))
    for _ in range(nrand):
        pts.append((rng.uniform(x0lo, x0hi), rng.uniform(llo, lhi), rng.uniform(wlo, whi)))
    return pts


stats = dict(A_boxes=0, A_mismatch=0, B_checks=0, B_maxratio=0.0, B_viol=0, C_checks=0, C_viol=0, C_minmargin=np.inf,
             D_signed=0, D_checks=0, D_viol=0)
worstB = None; worstC = None; worstD = None
t_start = time.time()
for ib in range(NBOX):
    chart = ib % 2
    scale = [0.002, 0.01, 0.04, 0.15][ib % 4] if chart == 0 else [0.002, 0.01, 0.03, 0.1][ib % 4]
    box = random_box(chart, scale)
    cval, cabs = (cval_a, cabs_a) if chart == 0 else (cval_f, cabs_f)
    M = 64
    left = np.arange(M) * (2 * math.pi / M); width = np.full(M, 2 * math.pi / M)
    um = left + width / 2; hh = width / 2
    # ---- A: numba signs vs port signs
    cmid, crad, dcen, dwid, cmax, rr = nb_cert.box_coeffs(chart, *box[:6])
    sg_nb, _ = nb_cert.arc_signs(cmid, crad, dcen, dwid, cmax, rr, box[6], box[7], um, hh, J)
    sg_po = signs_port(chart, box, um, hh, J, cval, cabs)
    stats['A_boxes'] += 1
    stats['A_mismatch'] += int(np.sum(sg_nb != sg_po))
    # ---- B, C: enclosures from the port
    pieces = box_coeffs_port(chart, *box[:6], cval, cabs)
    prm = sample_params(box)
    for t in (box[6], box[7]):
        bd = arc_bounds_port(*pieces, t, um, hh, J)
        for (x0, l, w) in prm:
            vals, mag = ft_derivs(chart, x0, l, w, t, um, J + 1)
            dev = np.abs(vals - bd['vc'])                       # [j, i]
            ratio = dev / np.maximum(bd['rad'], 1e-300)
            noise = 1e-13 * (mag * 6.0 ** np.arange(J + 2))[:, None]
            stats['B_checks'] += ratio.size
            rmax = float(ratio.max())
            if rmax > stats['B_maxratio']:
                stats['B_maxratio'] = rmax
            bad = (dev > bd['rad'] + noise)
            if bad.any():
                stats['B_viol'] += int(bad.sum())
                j, i = np.argwhere(bad)[0]
                worstB = (chart, box, t, (x0, l, w), int(j), float(um[i]), float(dev[j, i]), float(bd['rad'][j, i]))
    # C: points inside arcs, t inside [t0, t1] (affine: check both ends + random t)
    for (x0, l, w) in prm[::3]:
        for t in (box[6], box[7], rng.uniform(box[6], box[7])):
            uu = np.concatenate([um - hh, um + hh, um + rng.uniform(-1, 1, M) * hh])
            idx = np.concatenate([np.arange(M)] * 3)
            vals, mag = ft_derivs(chart, x0, l, w, t, uu, J)
            # enclosure claims are per t-end; for interior t use the affine combination of the two end enclosures
            bd0 = arc_bounds_port(*pieces, box[6], um, hh, J); bd1 = arc_bounds_port(*pieces, box[7], um, hh, J)
            lam = 0.0 if box[7] == box[6] else (t - box[6]) / (box[7] - box[6])
            lo = (1 - lam) * bd0['lo'] + lam * bd1['lo']; hi = (1 - lam) * bd0['hi'] + lam * bd1['hi']
            noise = 1e-13 * (mag * 6.0 ** np.arange(J + 1))[:, None]
            v = vals; L = lo[:, idx]; H = hi[:, idx]
            stats['C_checks'] += v.size
            margin = np.minimum(v - L, H - v) / np.maximum(H - L, 1e-300)
            stats['C_minmargin'] = min(stats['C_minmargin'], float(margin.min()))
            bad = (v < L - noise) | (v > H + noise)
            if bad.any():
                stats['C_viol'] += int(bad.sum())
                j, k = np.argwhere(bad)[0]
                worstC = (chart, box, t, (x0, l, w), int(j), float(uu[k]), float(v[j, k]), float(L[j, k]), float(H[j, k]))
    # ---- D: numba sign claims at sampled points
    signed = np.argwhere(sg_nb != 0)
    stats['D_signed'] += len(signed)
    for (x0, l, w) in prm[::2]:
        t = rng.uniform(box[6], box[7]) if rng.uniform() < 0.6 else box[6 + rng.integers(2)]
        uu = np.concatenate([um - hh, um + hh, um + rng.uniform(-1, 1, M) * hh])
        idx = np.concatenate([np.arange(M)] * 3)
        vals, mag = ft_derivs(chart, x0, l, w, t, uu, J)
        noise = 1e-13 * (mag * 6.0 ** np.arange(J + 1))[:, None]
        S = sg_nb[:, idx]
        bad = ((S == 1) & (vals < -noise)) | ((S == -1) & (vals > noise))
        stats['D_checks'] += int(np.sum(S != 0))
        if bad.any():
            stats['D_viol'] += int(bad.sum())
            j, k = np.argwhere(bad)[0]
            worstD = (chart, box, t, (x0, l, w), int(j), float(uu[k]), float(vals[j, k]), int(S[j, k]))
print(f"{NBOX} random boxes (both charts, 4 size scales), {time.time()-t_start:.0f}s")
for k, v in stats.items():
    print(f"  {k:14s} {v}")
print("  worst B violation:", worstB)
print("  worst C violation:", worstC)
print("  worst D violation:", worstD)

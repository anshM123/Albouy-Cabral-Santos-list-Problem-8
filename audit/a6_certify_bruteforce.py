"""AUDIT a6: end-to-end empirical soundness of nb_cert.certify_box.

For random boxes (both charts) and for boxes centred at configurations with 10 real roots, run
certify_box(chart, box, J=6, M0=64, rounds=8) exactly as nb_run does.  When it returns (True, B), sample the box
(corners, face centres, random points; t at both ends and inside) and count the real roots of f_t on the circle
with an independent evaluator (closed forms + FFT + companion polynomial; ambiguous near-circle roots counted as
real, i.e. the count is biased UPWARDS so that a violation cannot be hidden).  Any sampled count > B is a violation.
Also re-runs the refinement loop of certify_box in python (port) and checks every final arc sign at sampled points.
usage: python a6_certify_bruteforce.py NRANDOM NTEN SEED
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
from nbport import box_coeffs_port, arc_bounds_port, phipsi_coeffs, ft_derivs

NR = int(sys.argv[1]) if len(sys.argv) > 1 else 100
NT = int(sys.argv[2]) if len(sys.argv) > 2 else 20
rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
J = 6


def count_real(chart, x0, l, w, t):
    (aP, sP), (aS, sS) = phipsi_coeffs(chart, x0, l, w)
    a = (1 - t) * aP - t * aS; s = (1 - t) * sP - t * sS
    c = np.zeros(13, complex); c[6] = a[0]
    for k in range(1, 7):
        c[6 + k] = (a[k] - 1j * s[k]) / 2; c[6 - k] = (a[k] + 1j * s[k]) / 2
    mx = np.abs(c).max()
    top = 12
    while top > 0 and abs(c[top]) < 1e-15 * mx:
        top -= 1
    low = 0
    while low < top and abs(c[low]) < 1e-15 * mx:
        low += 1
    z = np.roots(c[low:top + 1][::-1])
    lz = np.abs(np.log(np.abs(z)))
    n_on = int(np.sum(lz < 1e-6)); n_amb = int(np.sum((lz >= 1e-6) & (lz < 1e-3)))
    return n_on, n_on + n_amb


def port_certify(chart, box):
    """python replica of certify_box's refinement loop, returns (ok, total, left, width, sg)"""
    cval, cabs = (cval_a, cabs_a) if chart == 0 else (cval_f, cabs_f)
    pieces = box_coeffs_port(chart, *box[:6], cval, cabs)
    M0, rounds = 64, 8
    left = np.arange(M0) * (2 * math.pi / M0); width = np.full(M0, 2 * math.pi / M0)
    total = 999
    for rd in range(rounds + 1):
        um = left + width / 2; hh = width / 2
        pos = np.ones((J + 1, um.size), bool); neg = np.ones((J + 1, um.size), bool)
        hopeless = np.ones(um.size, bool)
        for t in (box[6], box[7]):
            bd = arc_bounds_port(*pieces, t, um, hh, J)
            pos &= bd['lo'] > 0; neg &= bd['hi'] < 0
            hopeless &= np.all(np.abs(bd['vc'][:J + 1]) <= bd['rad'][:J + 1], axis=0)
        sg = np.where(pos, 1, np.where(neg, -1, 0)).astype(np.int64)
        total = nb_cert.dp_count(sg)
        if total <= 10:
            return True, total, left, width, sg
        if rd == rounds or hopeless.any():
            break
        sel = sg[0] == 0
        if not sel.any():
            break
        nl = []; nw = []
        for i in range(left.size):
            if sel[i]:
                nl += [left[i], left[i] + width[i] / 2]; nw += [width[i] / 2, width[i] / 2]
            else:
                nl.append(left[i]); nw.append(width[i])
        left = np.array(nl); width = np.array(nw)
    return False, total, left, width, sg


def box_samples(box, nrand):
    pts = list(itertools.product(box[0:2], box[2:4], box[4:6]))
    for _ in range(nrand):
        pts.append((rng.uniform(*box[0:2]), rng.uniform(*box[2:4]), rng.uniform(*box[4:6])))
    return pts


def make_box(chart, centre, r):
    x0, l, w, t = centre
    lim = (0.2014, 1.2532) if chart == 0 else (0.9, 1.0)
    tmax = 1.0 if chart == 0 else 0.9
    return (max(x0 - r[0], lim[0]), min(x0 + r[0], lim[1]), max(l - r[1], 0.0), min(l + r[1], 1.0),
            max(w - r[2], 0.0), min(w + r[2], np.pi / 2), max(t - r[3], 0.0), min(t + r[3], tmax))


def find_ten(chart):
    """random search for a parameter point with 10 real roots (strict count)"""
    while True:
        if chart == 0:
            x0 = rng.uniform(0.4, 1.25); l = rng.uniform(0.3, 1); w = rng.uniform(0, np.pi / 2)
        else:
            x0 = rng.uniform(0.9, 0.999); l = rng.uniform(0.3, 1); w = rng.uniform(0, np.pi / 2)
        for t in rng.uniform(0.02, 0.98 if chart == 0 else 0.9, 40):
            n, nb = count_real(chart, x0, l, w, t)
            if n == 10 and nb == 10:
                return (x0, l, w, t)


stats = dict(boxes=0, certified=0, viol=0, sign_viol=0, nb_port_total_mismatch=0, hist={}, slack={})
worst = []
t0 = time.time()
jobs = []
for i in range(NR):
    chart = i % 2
    if chart == 0:
        centre = (rng.uniform(0.2014, 1.2532), rng.uniform(0, 1), rng.uniform(0, np.pi / 2), rng.uniform(0, 1))
        r = np.array([0.02, 0.012, 0.02, 0.012]) * np.exp(rng.uniform(-2.5, 0, 4))
    else:
        centre = (rng.uniform(0.9, 1.0), rng.uniform(0, 1), rng.uniform(0, np.pi / 2), rng.uniform(0, 0.9))
        r = np.array([0.012, 0.012, 0.02, 0.012]) * np.exp(rng.uniform(-2.5, 0, 4))
    jobs.append((chart, make_box(chart, centre, r), 'random'))
for i in range(NT):
    chart = i % 2
    centre = find_ten(chart)
    r = np.array([0.004, 0.004, 0.006, 0.004]) * np.exp(rng.uniform(-3, 0, 4))
    jobs.append((chart, make_box(chart, centre, r), 'ten'))

for chart, box, kind in jobs:
    ok, B = nb_cert.certify_box(chart, *box, J, 64, 8)
    okp, Bp, left, width, sg = port_certify(chart, box)
    stats['boxes'] += 1
    if (ok, int(B)) != (okp, int(Bp)):
        stats['nb_port_total_mismatch'] += 1
    if not ok:
        continue
    stats['certified'] += 1
    stats['hist'][(kind, int(B))] = stats['hist'].get((kind, int(B)), 0) + 1
    mx = 0
    for (x0, l, w) in box_samples(box, 60):
        for t in (box[6], box[7], rng.uniform(box[6], box[7])):
            n, nb_ = count_real(chart, x0, l, w, t)
            mx = max(mx, nb_)
            if nb_ > B:
                stats['viol'] += 1
                worst.append((chart, box, (x0, l, w, t), n, nb_, int(B)))
    stats['slack'][int(B) - mx] = stats['slack'].get(int(B) - mx, 0) + 1
    # sign claims on the final refined arcs (port arcs; port == numba checked above)
    um = left + width / 2; hh = width / 2
    for (x0, l, w) in box_samples(box, 6)[::2]:
        t = rng.uniform(box[6], box[7])
        uu = np.concatenate([um - hh, um + hh, um + rng.uniform(-1, 1, um.size) * hh])
        idx = np.concatenate([np.arange(um.size)] * 3)
        vals, mag = ft_derivs(chart, x0, l, w, t, uu, J)
        noise = 1e-13 * (mag * 6.0 ** np.arange(J + 1))[:, None]
        S_ = sg[:, idx]
        bad = ((S_ == 1) & (vals < -noise)) | ((S_ == -1) & (vals > noise))
        stats['sign_viol'] += int(bad.sum())
    # coverage of the circle by the final arcs (bookkeeping)
    gaps = np.max(np.abs((left[1:] - (left[:-1] + width[:-1]))))
    end_gap = 2 * math.pi - (left[-1] + width[-1])
    stats['max_internal_gap'] = max(stats.get('max_internal_gap', 0.0), float(gaps))
    stats['max_end_gap'] = max(stats.get('max_end_gap', 0.0), float(end_gap))

print(f"{len(jobs)} boxes in {time.time()-t0:.0f}s")
for k, v in stats.items():
    if isinstance(v, dict):
        v = dict(sorted(v.items(), key=lambda kv: str(kv[0])))
    print(f"  {k}: {v}")
for wv in worst[:10]:
    print("  VIOLATION", wv)

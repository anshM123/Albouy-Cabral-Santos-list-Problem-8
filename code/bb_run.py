"""Unified parallel branch-and-bound driver with full split-tree recording (replaces nb_run.py, nb_run_chart.py and
nb_run_corner.py for the final runs).

usage:
  python bb_run.py --engine E --dom x0lo x0hi x1lo x1hi x2lo x2hi x3lo x3hi --n n0 n1 n2 n3 --tmax SEC --workers W
                   --out OUT.json [--minw 2e-5] [--exclude 8 numbers]... [--boxes BOXES.json]
  E in {a, f, C1, C1p, C1pp, K}.  K = e -> 1 corner engine (nb_corner, PHI0 = 0.9, PHI_IN = 0.12).

Split rule (deterministic): split the box at the float midpoint of the coordinate with the largest normalised width
(width / SCALE[i]; ties -> lowest index); the children are pushed (lower, upper), so the upper child is processed first.

Tree code, one character per popped box in processing order (a pre-order traversal):
  '+'  split            '!'  failed (normalised width below minw)
  'a'..'l'  certified leaf, bound 0..11, regular certifier (certify_box, M0 = 64, rounds = 8); a bound of 11
            means at most 10 by the parity lemma (the real-zero count with multiplicity is even)
  'A'..'L'  certified leaf, bound 0..11, deep certifier (certify_box_deep, M0 = 64, rounds = DEEP_ROUNDS, cap = DEEP_CAP)
  '1'..'6'  certified leaf by the winding certificate (wind_cert.wind_box): W = digit roots of Q in |z| < e^L, so at
            most 12 - 2W real roots; the test-circle log-radii L are stored in order in the record's 'wL' list
Replaying the code from the initial box with the same split rule reproduces every leaf exactly.  An initial box that
hits the time cap is stored with done = False and its complete remaining stack ('left'); resume jobs start from it.
Every leaf certificate (regular or deep) is reproduced by certify_box_deep with the recorded parameters, because both
use the same arc-refinement sequence and the regular one only stops earlier.
"""
import argparse
import json
import math
import sys
import time
import numpy as np
from multiprocessing import Pool

JJ = 6
M0 = 64
ROUNDS = 8
DEEP_ROUNDS = 16
DEEP_CAP = 8192
DEEP_MAXW = 1e-2     # deep fallback only for boxes whose largest normalised width is below this
PHI0 = 0.9
PHI_IN = 0.12
PHI_A = 2 * math.atan(3.0 * 0.37) - 1e-9      # away arc for the C3 box (c3/c3_run.py: Y0 = 3, NU1 = 0.37)

ENG = None
SCALE = None


WIND_M0, WIND_ROUNDS, WIND_CAP = 64, 8, 16384


def init(engine, scale, deep_cap=None, deep_rounds=None, order2=False):
    global ENG, SCALE, reg, deep, wind, DEEP_CAP, DEEP_ROUNDS
    wind = None
    ENG = engine; SCALE = scale
    if deep_cap:
        DEEP_CAP = deep_cap
    if deep_rounds:
        DEEP_ROUNDS = deep_rounds
    if order2 and engine in ('a', 'f'):
        # second-order parameter enclosure (nb_cert2): min(mean-value bound, second-order bound) everywhere
        import nb_cert2
        ch = {'a': 0, 'f': 1}[engine]
        reg = lambda bx: nb_cert2.certify_box2(ch, *bx, JJ, M0, ROUNDS)
        deep = lambda bx: nb_cert2.certify_box_deep2(ch, *bx, JJ, M0, DEEP_ROUNDS, DEEP_CAP)
        wind = make_wind(lambda bx: nb_cert2.box_data2(ch, *bx[:6]), second=True)
    elif engine in ('a', 'f', 'fm5', 'fm8'):
        import nb_cert
        ch = {'a': 0, 'f': 1, 'fm5': 3, 'fm8': 4}[engine]
        reg = lambda bx: nb_cert.certify_box(ch, *bx, JJ, M0, ROUNDS)
        deep = lambda bx: nb_cert.certify_box_deep(ch, *bx, JJ, M0, DEEP_ROUNDS, DEEP_CAP)
        wind = make_wind(lambda bx: nb_cert.box_coeffs(ch, *bx[:6]))
    elif engine in ('C1', 'C1p', 'C1pp', 'B1', 'B2', 'B3p', 'B3m', 'L2E', 'L2P', 'L2Ap', 'L2Am', 'L2Tp', 'L2Tm',
                    'D2E', 'D2E2', 'D2P', 'D2K', 'D2Tp', 'D2Tm'):
        import importlib
        m = importlib.import_module('nbeng_' + engine)
        reg = lambda bx: m.certify_box(*bx, JJ, M0, ROUNDS)
        deep = lambda bx: m.certify_box_deep(*bx, JJ, M0, DEEP_ROUNDS, DEEP_CAP)
        if engine in ('C1', 'C1p', 'C1pp'):
            wind = make_wind(lambda bx: m.box_coeffs(*bx[:6]))
        elif engine in ('B1', 'B2', 'L2E', 'L2P', 'L2Ap', 'L2Am', 'D2E', 'D2E2', 'D2P', 'D2K'):
            wind = make_wind(lambda bx: m.box_coeffs(*bx))
    elif engine == 'K':
        import nb_corner

        def reg(bx):
            r = nb_corner.certify_corner_box(*bx, PHI0, PHI_IN)
            return r[0], r[1]
        deep = None
    elif engine == 'KA':
        # C3 box, away arc |f - pi| >= PHI_A (= 2 atan(Y0 * NU1) of c3/c3_run.py, minus 1e-9): at most 2 zeros
        import c3_away

        def reg(bx):
            t0 = 1.0 / (1.0 + bx[7] * bx[7]); t1 = 1.0 / (1.0 + bx[6] * bx[6])
            ok, c = c3_away.away_count_t(bx[0], bx[1], bx[2], bx[3], bx[4], bx[5], t0, t1, PHI_A, 64, 16, JJ, 2.0, 8192)
            return ok, c
        deep = None
    else:
        raise SystemExit('unknown engine ' + engine)
    # compile once per worker on a box inside the respective domain
    probe = {'a': (0.5, 0.51, 0.5, 0.51, 0.5, 0.51, 0.3, 0.31), 'f': (0.96, 0.97, 0.5, 0.51, 0.5, 0.51, 0.3, 0.31),
             'fm5': (0.96, 0.97, 0.5, 0.51, 0.5, 0.51, 0.3, 0.31), 'fm8': (0.96, 0.97, 0.5, 0.51, 0.5, 0.51, 0.3, 0.31),
             'C1': (0.05, 0.06, 0.3, 0.31, 0.5, 0.51, 0.3, 0.31), 'C1p': (0.05, 0.06, 0.5, 0.51, 0.5, 0.51, 0.3, 0.31),
             'C1pp': (0.2, 0.21, 0.3, 0.31, 0.5, 0.51, 0.3, 0.31), 'K': (0.99, 0.995, 0.5, 0.55, 0.2, 0.25, 0.0, 0.02),
             'KA': (0.995, 0.996, 0.96, 0.97, 1.45, 1.46, 0.0, 0.01)}
    probe.update({c: (0.5, 0.51, 0.3, 0.31, 0.2, 0.21, 0.3, 0.31) for c in ('B1', 'L2E', 'L2P', 'L2Ap', 'L2Am', 'D2E', 'D2E2', 'D2P',
                                                                         'D2K', 'D2Tp', 'D2Tm')})
    probe.update({'B2': (0.05, 0.06, 0.2, 0.21, 0.5, 0.51, 0.3, 0.31), 'B3p': (0.05, 0.06, 0.2, 0.21, 0.5, 0.51, 0.5, 0.51),
                  'B3m': (0.05, 0.06, 0.2, 0.21, 0.5, 0.51, 0.5, 0.51), 'L2Tp': (0.5, 0.51, 0.3, 0.31, 0.2, 0.21, 0.3, 0.31),
                  'L2Tm': (0.5, 0.51, 0.3, 0.31, 0.2, 0.21, 0.3, 0.31)})
    reg(probe[engine])
    if deep is not None:
        deep(probe[engine])


def make_wind(bc, second=False):
    import wind_cert
    if second:
        import nb_cert2
        wbox = nb_cert2.wind_box2
    else:
        wbox = wind_cert.wind_box

    def w(bx):
        cm = bc(bx)
        L = wind_cert.choose_L(cm[0], 0.5 * (bx[6] + bx[7]))
        if L is None:
            return False, -1, None
        ok, W = wbox(*cm, bx[6], bx[7], L, WIND_M0, WIND_ROUNDS, WIND_CAP)
        return ok and W >= 1, W, L
    probe_box = (0.5, 0.51, 0.5, 0.51, 0.5, 0.51, 0.3, 0.31)
    try:
        wbox(*bc(probe_box), 0.3, 0.31, -1.0, 64, 1, 100)
    except Exception:
        pass
    return w


SD2 = 5e-3


def prune(bx):
    """D2b charts only: True if every point of the box has l > 1 (p^2 + (1 - kappa)^2 > 1), i.e. is unphysical."""
    if not ENG.startswith('D2'):
        return False
    n0, n1 = bx[0], bx[1]
    if ENG in ('D2E', 'D2E2'):
        p_lo = SD2 * bx[2] * n0; k_hi = SD2 * bx[5] * n1
    elif ENG == 'D2P':
        p_lo = SD2 * n0; k_hi = SD2 * bx[5] * n1
    elif ENG == 'D2K':
        p_lo = SD2 * bx[4] * n0; k_hi = SD2 * n1
    else:
        p_lo = SD2 * bx[4] * n0; k_hi = SD2 * bx[7] * n1
    return p_lo * p_lo + (1 - k_hi) ** 2 > 1 + 1e-12


def split(bx):
    w = [(bx[2 * i + 1] - bx[2 * i]) / SCALE[i] for i in range(4)]
    i = int(np.argmax(w))
    mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx)
    a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)


def work(args):
    bx0, tmax, minw = args
    t0 = time.time()
    stack = [tuple(bx0)]; code = []; nleaf = 0; ndeep = 0; nwind = 0; wL = []; vol = 0.0; hist = {}; fails = []
    while stack:
        if time.time() - t0 > tmax:
            return dict(box=list(bx0), done=False, tree=''.join(code), nleaf=nleaf, ndeep=ndeep, nwind=nwind, wL=wL,
                        vol=vol, hist=hist, fails=fails, left=[list(b) for b in stack], sec=time.time() - t0)
        bx = stack.pop()
        if prune(bx):
            code.append('u'); nleaf += 1; hist['u'] = hist.get('u', 0) + 1
            continue
        good, tot = reg(bx)
        tot = int(tot)
        method = 0
        if not good and wind is not None:
            okw, W, L = wind(bx)
            if okw:
                code.append(chr(ord('0') + W)); wL.append(L)
                nleaf += 1; nwind += 1
                vol += float(np.prod([bx[2 * i + 1] - bx[2 * i] for i in range(4)]))
                hist[12 - 2 * W] = hist.get(12 - 2 * W, 0) + 1
                continue
        if not good and deep is not None and (tot < 999 or max((bx[2 * i + 1] - bx[2 * i]) / SCALE[i] for i in range(4)) < 1e-3) and                 max((bx[2 * i + 1] - bx[2 * i]) / SCALE[i] for i in range(4)) < DEEP_MAXW:
            good, tot = deep(bx)
            tot = int(tot)
            method = 1
        if good:
            code.append(chr((ord('A') if method else ord('a')) + tot))
            nleaf += 1; ndeep += method
            vol += float(np.prod([bx[2 * i + 1] - bx[2 * i] for i in range(4)]))
            hist[tot] = hist.get(tot, 0) + 1
            continue
        if min((bx[2 * i + 1] - bx[2 * i]) / SCALE[i] for i in range(4)) < minw:
            code.append('!')
            fails.append((list(bx), tot))
            continue
        code.append('+')
        stack.extend(split(bx))
    return dict(box=list(bx0), done=True, tree=''.join(code), nleaf=nleaf, ndeep=ndeep, nwind=nwind, wL=wL, vol=vol,
                hist=hist, fails=fails, left=[], sec=time.time() - t0)


def inside(bx, ex):
    return all(ex[2 * i] <= bx[2 * i] and bx[2 * i + 1] <= ex[2 * i + 1] for i in range(4))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('--engine', required=True)
    ap.add_argument('--dom', type=float, nargs=8)
    ap.add_argument('--n', type=int, nargs=4)
    ap.add_argument('--boxes', default=None)
    ap.add_argument('--exclude', type=float, nargs=8, action='append', default=[])
    ap.add_argument('--tmax', type=float, default=7200.0)
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--minw', type=float, default=2e-5)
    ap.add_argument('--scale', type=float, nargs=4, default=None)
    ap.add_argument('--out', required=True)
    ap.add_argument('--resume', action='store_true', help='skip initial boxes completed in OUT.records.jsonl')
    ap.add_argument('--deep-cap', type=int, default=None)
    ap.add_argument('--deep-rounds', type=int, default=None)
    ap.add_argument('--order2', action='store_true', help='second-order enclosure (nb_cert2) for engines a and f')
    A = ap.parse_args()
    if A.boxes:
        boxes = [tuple(b) for b in json.load(open(A.boxes))['boxes']]
        dom = None
    else:
        dom = tuple(A.dom)
        edges = [np.linspace(dom[2 * i], dom[2 * i + 1], A.n[i] + 1) for i in range(4)]
        boxes = [(edges[0][a], edges[0][a + 1], edges[1][b], edges[1][b + 1], edges[2][c], edges[2][c + 1],
                  edges[3][d], edges[3][d + 1])
                 for a in range(A.n[0]) for b in range(A.n[1]) for c in range(A.n[2]) for d in range(A.n[3])]
    boxes = [tuple(float(x) for x in b) for b in boxes]
    skipped = [b for b in boxes if any(inside(b, ex) for ex in A.exclude)]
    boxes = [b for b in boxes if not any(inside(b, ex) for ex in A.exclude)]
    if A.scale:
        scale = tuple(A.scale)
    elif A.engine == 'a':
        scale = (math.pi / 2, 1.0, math.pi / 2, 1.0)
    elif A.engine in ('f', 'fm5', 'fm8'):
        scale = (1.0, 1.0, math.pi / 2, 1.0)
    else:
        if dom is None:
            raise SystemExit('--scale is required with --boxes for this engine')
        scale = tuple(dom[2 * i + 1] - dom[2 * i] for i in range(4))
    if A.deep_cap:
        DEEP_CAP = A.deep_cap
    if A.deep_rounds:
        DEEP_ROUNDS = A.deep_rounds
    header = dict(engine=A.engine, domain=dom, n0=A.n, boxes_file=A.boxes, exclude=A.exclude, skipped=skipped,
                  scale=scale, minw=A.minw, J=JJ, M0=M0, rounds=ROUNDS, deep_rounds=DEEP_ROUNDS, deep_cap=DEEP_CAP,
                  deep_maxw=DEEP_MAXW, wind_m0=WIND_M0, wind_rounds=WIND_ROUNDS, wind_cap=WIND_CAP,
                  phi0=PHI0, phi_in=PHI_IN, phi_a=PHI_A, order2=A.order2, argv=sys.argv)
    res = []
    if A.resume:
        import os
        todo = set(tuple(b) for b in boxes); seen = set()
        if os.path.exists(A.out + '.records.jsonl'):
            for line in open(A.out + '.records.jsonl'):
                try:
                    rec = json.loads(line)
                except Exception:
                    continue
                # keep one completed record per initial box that is still in the run (not excluded)
                if rec['done'] and tuple(rec['box']) in todo and tuple(rec['box']) not in seen:
                    res.append(rec); seen.add(tuple(rec['box']))
        done_roots = set(tuple(r['box']) for r in res)
        boxes = [b for b in boxes if tuple(b) not in done_roots]
        print('resume: kept', len(res), 'completed records', flush=True)
    print('initial boxes', len(boxes), 'skipped (excluded)', len(skipped), flush=True)
    t0 = time.time()
    frec = open(A.out + '.records.jsonl', 'a' if A.resume else 'w')
    flog = open(A.out + '.fails.jsonl', 'a' if A.resume else 'w')
    with Pool(A.workers, initializer=init, initargs=(A.engine, scale, A.deep_cap, A.deep_rounds, A.order2)) as pool:
        for k, r in enumerate(pool.imap_unordered(work, [(b, A.tmax, A.minw) for b in boxes])):
            res.append(r)
            frec.write(json.dumps(r) + chr(10)); frec.flush()
            if r['fails'] or not r['done']:
                flog.write(json.dumps(dict(box=r['box'], done=r['done'], fails=r['fails'], left=r['left'])) + chr(10))
                flog.flush()
            if (k + 1) % 25 == 0 or k + 1 == len(boxes):
                print(f"[{time.time()-t0:7.0f}s] {k+1}/{len(boxes)} returned, {sum(x['done'] for x in res)} complete, "
                      f"{sum(x['nleaf'] for x in res)} leaves ({sum(x['ndeep'] for x in res)} deep, {sum(x['nwind'] for x in res)} wind), "
                      f"{sum(len(x['fails']) for x in res)} failed", flush=True)
    flog.close(); frec.close()
    json.dump(dict(header=header, results=res), open(A.out, 'w'))
    vol = sum(x['vol'] for x in res)
    totvol = sum(float(np.prod([b[2 * i + 1] - b[2 * i] for i in range(4)])) for b in boxes)
    print("certified volume fraction (of processed boxes):", vol / totvol if totvol > 0 else 1.0, "; incomplete initial boxes:",
          sum(not x['done'] for x in res), "; failed leaves:", sum(len(x['fails']) for x in res),
          "; wall", round(time.time() - t0), "s", flush=True)

"""Branch and bound for the C3 core charts (see gen_c3.py, C3_CORE.md), with split-tree recording as in bb_run.py.

usage: python c3_run.py CHART n0 n1 n2 n3 TMAX WORKERS OUT.json

Per chart: box ranges DOM[chart], root interval, target.  Leaf codes: 'a'+bound (certified), 'u' (box entirely
unphysical: l > 1), 'o' (box entirely outside the C3 parameter box), '+' split, '!' failed.
Level-1 charts Ca, CK: a box with E <= E2 and m <= M2 everywhere ("level-2 zone") is certified by counting the roots
of G on |y| in [Y2, Y0] only (target ZONE_OUTER = 5); the roots with |y| <= Y2 are bounded by 4 by the level-2
charts.  All other level-1 boxes: roots on |y| <= Y0, target TARGET[chart] (Ca 9, CK/CE/Cm 8).  Intermediate charts PhiP/PhiM and level-2 charts *Yp/*Ym: no zero (0).
Level-2 charts *E/*m: roots on |Y| <= R2, target 4.
C3 parameter box (the region excluded from the K run): Qz <= QZMAX, kappa <= KAPMAX, eps <= EPSMAX, mu <= MUMAX.
"""
import json
import math
import sys
import time
import numpy as np
from multiprocessing import Pool
import os as _os
sys.path.insert(1, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))   # code/: tau_wind.py

Y0 = 3.0
Y2, E2, M2 = 0.25, 0.01, 0.01
R1 = 1.001          # level-1 / intermediate cone ratio bound (exact cones need 1)
R2 = 1.001          # level-2 cone ratio bound
NU1 = 0.37          # radial range of the level-1 charts (C3 box has weighted radius <= 0.36963)
ETA1 = 1.0001       # radial range of the level-2 charts
QZMAX = math.cos(1.413716694115407)               # 0.156434...
KAPMAX = 1 - 0.95 * math.sin(1.413716694115407)   # 0.0616958...
EPSMAX = 0.01
MUMAX = 0.0505
J, M0, ROUNDS, CAP = 6, 64, 16, 8192
ROUNDS1, CAP1 = 6, 1024        # first (cheap) pass; ROUNDS/CAP only for small boxes (normalised width < DEEP_W)
DEEP_W = 2e-3
DEEP_ALWAYS = False            # set True by re-certification tools: always try ROUNDS/CAP after a failed cheap pass
# parity lemma: near + intermediate (0) + away (<= 2, run KA) = total; an even total <= 11 is <= 10, so the near target
# may be 9 (Ca) and the zone outer target 5 (inner <= 4 from the level-2 charts).
ZONE_OUTER = 5

DOM = {
    'Ca': (0, NU1, 0, R1**2, 0, R1**2, 0, R1**3), 'CK': (0, NU1, 0, R1, 0, R1**2, 0, R1**3),
    'CE': (0, NU1, 0, R1, 0, R1**2, 0, R1**3), 'Cm': (0, NU1, 0, R1, 0, R1**2, 0, R1**2),
    'PhiP': (0, NU1, 0, R1, 0, R1**2, 0, R1**2), 'PhiM': (0, NU1, 0, R1, 0, R1**2, 0, R1**2),
    'L2aE': (0, NU1, 0, R1**2, 0, ETA1, 0, R2**2), 'L2am': (0, NU1, 0, R1**2, 0, ETA1, 0, R2**2),
    'L2aYp': (0, NU1, 0, R1**2, 0, ETA1, 0, R2**2), 'L2aYm': (0, NU1, 0, R1**2, 0, ETA1, 0, R2**2),
    'L2KE': (0, NU1, 0, R1, 0, ETA1, 0, R2**2), 'L2Km': (0, NU1, 0, R1, 0, ETA1, 0, R2**2),
    'L2KYp': (0, NU1, 0, R1, 0, ETA1, 0, R2**2), 'L2KYm': (0, NU1, 0, R1, 0, ETA1, 0, R2**2),
}
ROOT = {'Ca': (-Y0, Y0), 'CK': (-Y0, Y0), 'CE': (-Y0, Y0), 'Cm': (-Y0, Y0), 'PhiP': (0, R1**3), 'PhiM': (0, R1**3),
        'L2aE': (-R2, R2), 'L2am': (-R2, R2), 'L2KE': (-R2, R2), 'L2Km': (-R2, R2),
        'L2aYp': (0, R2**2), 'L2aYm': (0, R2**2), 'L2KYp': (0, R2**2), 'L2KYm': (0, R2**2)}
TARGET = {'Ca': 9, 'CK': 8, 'CE': 8, 'Cm': 8, 'PhiP': 0, 'PhiM': 0, 'L2aE': 4, 'L2am': 4, 'L2KE': 4, 'L2Km': 4,
          'L2aYp': 0, 'L2aYm': 0, 'L2KYp': 0, 'L2KYm': 0}
# Face charts (e -> 1 and mu -> 0 at general Pz; C3_CORE.md section 6): weighted radius s = max(sqrt(eps), sqrt(mu)) <= SF,
# box variables (s, ratio, Pz, Qz) with Pz, Qz in [0, 1]; near count on |y| <= YF (th = s y), intermediate charts on
# th = +-YF s (no zero for the ratio m in [0, R1]).  The face region also needs the outer count of kq_run --face.
SF = 0.02
YF = 3.0
DOM.update({'FE': (0, SF * R1, 0, R1, 0, 1, 0, 1), 'FM': (0, SF * R1, 0, R1, 0, 1, 0, 1),
            'FIP': (0, SF * R1, 0, R1, 0, 1, 0, 1), 'FIM': (0, SF * R1, 0, R1, 0, 1, 0, 1)})
ROOT.update({'FE': (-YF, YF), 'FM': (-YF, YF), 'FIP': (0, R1), 'FIM': (0, R1)})
TARGET.update({'FE': 5, 'FM': 5, 'FIP': 0, 'FIM': 0})
C3_L, C3_W = 0.95, 9 * math.pi / 20          # the C3 cell: l >= C3_L and w >= C3_W (covered by the C3 core charts)

CH = None
SCALE = None


def engine_module(chart):
    """the second-order engine (make_c3_2.py) where it has been generated, else the make_c3.py engine"""
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    return 'nbeng_c3_%s_2' % chart if os.path.exists(os.path.join(here, 'nbeng_c3_%s_2.py' % chart)) else 'nbeng_c3_' + chart


TW = False          # --tw: tau-winding certificate (code/tau_wind.py) for the charts Ca, CK, CE, Cm, leaf code '7'
TWENG = None


def init(chart, scale=None, tw=False):
    global CH, SCALE, cnt, TW, TWENG
    import importlib
    cnt = importlib.import_module(engine_module(chart)).count
    TW = tw
    if tw:
        import tau_wind
        if chart not in tau_wind.TW_CHARTS:
            raise SystemExit('--tw is only defined for the charts ' + ', '.join(tau_wind.TW_CHARTS))
        TWENG = importlib.import_module('nbeng_c3_%s_2' % chart)
    CH = chart
    d = DOM[chart]
    SCALE = tuple(d[2 * i + 1] - d[2 * i] for i in range(4)) if scale is None else tuple(scale)
    cnt(0.1, 0.11, 0.5, 0.51, 0.3, 0.31, 0.2, 0.21, ROOT[chart][0], ROOT[chart][1], J, M0, 2, 100, 0)


def face_bounds(bx):
    """face charts: (unphysical, outside) for a box (s, ratio, Pz, Qz).  unphysical: Pz^2 + Qz^2 > 1 on the whole box;
    outside: the whole box lies in the C3 cell (l >= C3_L and w >= C3_W), which the C3 core charts cover, or has
    s_lo > SF (only s <= SF is needed: C3_CORE.md section 6)."""
    pz0, pz1, qz0, qz1 = bx[4], bx[5], bx[6], bx[7]
    unphys = pz0 * pz0 + qz0 * qz0 > 1.0 * (1 + 1e-12)
    l_lo = math.sqrt(pz0 * pz0 + qz0 * qz0) * (1 - 1e-12)
    w_lo = math.atan2(pz0, qz1) * (1 - 1e-12) if qz1 > 0 else math.pi / 2
    outside = (l_lo >= C3_L and w_lo >= C3_W) or bx[0] > SF
    return unphys, outside


def phys_bounds(ch, bx):
    """lower bounds of (Qz, kappa, eps, mu) over the box, and the max over the box of the l <= 1 slack
    (Qz^2 <= 2 kappa - kappa^2) expressed so that slack < 0 everywhere means unphysical."""
    n0, n1 = bx[0], bx[1]
    if ch in ('Ca', 'L2aE', 'L2am', 'L2aYp', 'L2aYm'):
        K0, K1 = bx[2], bx[3]
        lo = (n0, n0 * n0 * K0)
        unphys = 2 * K1 - n0 * n0 * K1 * K1 < 1.0          # needs 1 <= 2K - nu^2 K^2 (max at K1, nu0)
    elif ch in ('CK', 'L2KE', 'L2Km', 'L2KYp', 'L2KYm'):
        a0 = bx[2]
        lo = (n0 * a0, n0 * n0)
        unphys = False                                    # a <= R1: a^2 <= 2 - nu^2 always
    else:
        a0, K1 = bx[2], bx[5]
        K0 = bx[4]
        lo = (n0 * a0, n0 * n0 * K0)
        unphys = a0 * a0 > 2 * K1 - n0 * n0 * K1 * K1     # a^2 <= 2K - nu^2 K^2 violated everywhere
    # eps, mu lower bounds
    if ch == 'Ca':
        eps0, mu0 = n0 * n0 * bx[4], n0 ** 3 * bx[6]
    elif ch == 'CK':
        eps0, mu0 = n0 * n0 * bx[4], n0 ** 3 * bx[6]
    elif ch == 'CE':
        eps0, mu0 = n0 * n0, n0 ** 3 * bx[6]
    elif ch == 'Cm':
        eps0, mu0 = n0 * n0 * bx[6], n0 ** 3
    elif ch in ('PhiP', 'PhiM'):
        eps0, mu0 = n0 * n0 * bx[6], 0.0
    elif ch.startswith('L2') and ch[3:] == 'E':
        eps0, mu0 = n0 * n0 * E2 * bx[4] ** 2, 0.0
    elif ch.startswith('L2') and ch[3:] == 'm':
        eps0, mu0 = 0.0, n0 ** 3 * M2 * bx[4] ** 2
    else:
        eps0, mu0 = 0.0, 0.0
    outside = lo[0] > QZMAX * (1 + 1e-9) or lo[1] > KAPMAX * (1 + 1e-9) or eps0 > EPSMAX * (1 + 1e-9) or \
        mu0 > MUMAX * (1 + 1e-9)
    return unphys, outside


def cnt2(bx, r0, r1, target):
    """root count on [r0, r1]: first with ROUNDS1/CAP1, then (small boxes, or DEEP_ALWAYS) with ROUNDS/CAP"""
    ok, c = cnt(*bx, r0, r1, J, M0, ROUNDS1, CAP1, target)
    if ok:
        return ok, c
    if DEEP_ALWAYS or max((bx[2 * i + 1] - bx[2 * i]) / SCALE[i] for i in range(4)) < DEEP_W:
        ok2, c2 = cnt(*bx, r0, r1, J, M0, ROUNDS, CAP, target)
        if ok2 or c2 < c:
            return ok2, c2
    return ok, c


def certify(bx):
    unphys, outside = face_bounds(bx) if CH in ('FE', 'FM', 'FIP', 'FIM') else phys_bounds(CH, bx)
    if unphys:
        return 'u', -1
    if outside:
        return 'o', -1
    r0, r1 = ROOT[CH]
    if CH in ('Ca', 'CK') and bx[5] <= E2 and bx[7] <= M2:
        ok1, c1 = cnt2(bx, -Y0, -Y2, ZONE_OUTER)
        if not ok1:
            return None, int(c1)
        ok2, c2 = cnt2(bx, Y2, Y0, ZONE_OUTER - int(c1))
        c = int(c1) + int(c2)
        return ('Z', c) if ok2 and c <= ZONE_OUTER else (None, c)
    ok, c = cnt2(bx, r0, r1, TARGET[CH])
    return ('a', int(c)) if ok else (None, int(c))


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
    stack = [tuple(bx0)]; code = []; nleaf = 0; hist = {}; fails = []; twl = []
    while stack:
        if time.time() - t0 > tmax:
            return dict(box=list(bx0), done=False, tree=''.join(code), nleaf=nleaf, hist=hist, fails=fails,
                        left=[list(b) for b in stack], tw=twl)
        bx = stack.pop()
        kind, c = certify(bx)
        if kind in ('u', 'o'):
            code.append(kind); nleaf += 1; hist[kind] = hist.get(kind, 0) + 1
            continue
        if kind == 'a':
            code.append(chr(ord('a') + c)); nleaf += 1; hist[c] = hist.get(c, 0) + 1
            continue
        if kind == 'Z':            # level-2 zone leaf: outer count c (<= ZONE_OUTER = 5); inner <= 4 by the level-2 charts
            code.append(chr(ord('A') + c)); nleaf += 1; hist['Z%d' % c] = hist.get('Z%d' % c, 0) + 1
            continue
        if TW:                     # tau-winding certificate: <= 10 real zeros on the whole circle for every point of bx
            import tau_wind
            okt, cre, cim, rr_ = tau_wind.certify_tw_chart(TWENG, bx)
            if okt:
                code.append('7'); twl.append([0, cre, cim, rr_]); nleaf += 1; hist['t'] = hist.get('t', 0) + 1
                continue
        if min((bx[2 * i + 1] - bx[2 * i]) / SCALE[i] for i in range(4)) < minw:
            code.append('!'); fails.append((list(bx), c))
            continue
        code.append('+')
        stack.extend(split(bx))
    return dict(box=list(bx0), done=True, tree=''.join(code), nleaf=nleaf, hist=hist, fails=fails, left=[], tw=twl)


if __name__ == '__main__':
    chart = sys.argv[1]
    minw = 1e-5
    bfile = None; scale = None
    if sys.argv[2] == '--boxes':
        # continuation run: python c3_run.py CHART --boxes BOXES.json --scale s0 s1 s2 s3 TMAX WORKERS OUT.json
        bfile = sys.argv[3]
        scale = tuple(float(v) for v in sys.argv[5:9])
        tmax = float(sys.argv[9]); workers = int(sys.argv[10]); out = sys.argv[11]
        boxes = [tuple(float(v) for v in b) for b in json.load(open(bfile))['boxes']]
        dom = None; n = None
    else:
        n = [int(v) for v in sys.argv[2:6]]
        tmax = float(sys.argv[6]); workers = int(sys.argv[7]); out = sys.argv[8]
        dom = DOM[chart]
        edges = [np.linspace(dom[2 * i], dom[2 * i + 1], n[i] + 1) for i in range(4)]
        boxes = [tuple(float(v) for v in (edges[0][p], edges[0][p + 1], edges[1][q], edges[1][q + 1], edges[2][r],
                                          edges[2][r + 1], edges[3][s], edges[3][s + 1]))
                 for p in range(n[0]) for q in range(n[1]) for r in range(n[2]) for s in range(n[3])]
    header = dict(chart=chart, engine_module=engine_module(chart), domain=dom, n0=n, root=ROOT[chart], target=TARGET[chart], Y0=Y0, Y2=Y2, E2=E2, M2=M2,
                  R1=R1, R2=R2, NU1=NU1, ETA1=ETA1, QZMAX=QZMAX, KAPMAX=KAPMAX, EPSMAX=EPSMAX, MUMAX=MUMAX, J=J, M0=M0,
                  rounds=ROUNDS, cap=CAP, rounds1=ROUNDS1, cap1=CAP1, deep_w=DEEP_W, minw=minw,
                  scale=list(scale) if scale is not None else [dom[2 * i + 1] - dom[2 * i] for i in range(4)],
                  boxes_file=bfile, tw='--tw' in sys.argv, argv=sys.argv)
    t0 = time.time(); res = []
    import os
    if '--resume' in sys.argv and os.path.exists(out + '.records.jsonl'):    # keep completed initial boxes
        keep = {}
        for line in open(out + '.records.jsonl'):
            try:
                rec = json.loads(line)
            except Exception:
                continue
            # keep every record returned by a worker (complete, or capped with its open stack); one per box,
            # a complete record taking precedence
            b_ = tuple(rec['box'])
            if rec['done'] or b_ not in keep:
                keep[b_] = rec
        res = list(keep.values())
        done = set(keep)
        boxes = [b for b in boxes if b not in done]
        print('resume: kept', len(res), 'completed records', flush=True)
    frec = open(out + '.records.jsonl', 'a' if '--resume' in sys.argv else 'w')
    with Pool(workers, initializer=init, initargs=(chart, scale, '--tw' in sys.argv)) as pool:
        for k, r in enumerate(pool.imap_unordered(work, [(b, tmax, minw) for b in boxes])):
            res.append(r); frec.write(json.dumps(r) + chr(10)); frec.flush()
            if (k + 1) % 8 == 0 or k + 1 == len(boxes):
                print(f"[{time.time()-t0:7.0f}s] {k+1}/{len(boxes)} returned, {sum(x['done'] for x in res)} complete, "
                      f"{sum(x['nleaf'] for x in res)} leaves, {sum(len(x['fails']) for x in res)} failed", flush=True)
    frec.close()
    json.dump(dict(header=header, results=res), open(out, 'w'))
    H = {}
    for r in res:
        for kk, v in r['hist'].items():
            H[str(kk)] = H.get(str(kk), 0) + v
    print('done', chart, 'incomplete', sum(not x['done'] for x in res), 'failed', sum(len(x['fails']) for x in res),
          'hist', dict(sorted(H.items())), 'wall', round(time.time() - t0), flush=True)

"""Engine KQ for the near-parabolic region e in [0.95, 1) at moderate-to-large circles (f-chart variables).
Box variables (e, l, w, t) exactly as the f-chart (t = kappa^2/(1 + kappa^2), kappa = R/q), restricted to t >= T_MIN so
that mu = q/R = sqrt((1 - t)/t) is bounded.  With phi = f - pi:
  |phi| <= PHI_0 : real zeros of Q(th) = P(f)(1 + th^2)^6, th = tan(phi/2), written in eps = 1 - e (no cancellation as
                   e -> 1; the aphelion pair sits at th = O(eps^2) as a well-conditioned double cluster), monomial-basis
                   Rolle runs on th in [-tan(PHI_0/2), tan(PHI_0/2)]  (c3/nbeng_c3_KBQe_2: second-order enclosure,
                   variables (eps, Pz, Qz, mu));
  |phi| >= PHI_0 : real zeros of the f-chart polynomial on the away arc (nb_cert2.away_count_t2: the arc count of
                   c3_away.away_count_t with the second-order parameter enclosure of nb_cert2).
  A box may instead get the winding certificate of the f-chart polynomial on the whole circle (nb_cert2.wind_box2 =
  wind_cert.wind_box with the second-order enclosure; PROOF.md §2): W >= 1 roots of the self-inversive polynomial
  inside |z| < e^L give at most 12 - 2W <= 10 real zeros (leaf code '1'..'6', radius L stored in 'wL').
  P = mu^2 k^2 sin^2 f sigma^2 - D^2 B^2 is a positive multiple of the f-chart polynomial, so both parts count the same
  zeros with multiplicity; total <= 11 suffices by the parity lemma (an even count <= 11 is <= 10).
  --peri (C3_CORE.md section 7): the circle is instead covered by four overlapping parts, the aphelion polynomial on
  |f - pi| <= phi0, the perihelion polynomial G(tp) = tp^12 Q(-1/tp), tp = tan(f/2), on |f| <= PHI_P = 1.2 (c3 chart
  KBPe; G = P(f)(1 + tp^2)^6), and the f-chart arcs between them; the box gets the smaller of this count and the plain
  one.
usage: python kq_run.py elo ehi llo lhi wlo whi tlo thi n0 n1 n2 n3 TMAX WORKERS OUT.json [--resume]
       python kq_run.py --boxes BOXES.json --scale s0 s1 s2 s3 TMAX WORKERS OUT.json      (continuation run on the boxes
       of BOXES.json, e.g. the open boxes of a base run; the split rule uses the given scale, recorded in the header)"""
import json
import math
import sys
import time
import numpy as np
from multiprocessing import Pool
import splitrule

sys.path.insert(0, 'c3')
PHI_0 = 0.9
T_HI = math.tan(PHI_0 / 2) + 1e-9
FACE_S = 0.02                  # face neighbourhood: eps = 1 - e <= FACE_S^2 and mu <= FACE_S^2 (c3_run.SF)
FACE_Y = 3.0                   # face charts count |th| <= FACE_Y s; the outer count starts at T_X = FACE_Y FACE_S
T_X = FACE_Y * FACE_S
LIMIT_OUTER = 6                # outer count target (near count of the face charts <= 5): total <= 11, parity
C3_L, C3_W = 0.95, 9 * math.pi / 20
FACE_MODE = False
C3PRUNE = False                # --c3prune: boxes inside the C3 cell are pruned 'o' in the regular mode too
PHI_P = 1.2                    # --peri: perihelion part |f| <= PHI_P (c3 chart KBPe), C3_CORE.md section 7
PERI = None                    # None, or PHI_P when --peri is given
SPLIT_RULE = None               # None: widest normalised width; 'face': splitrule.py (--facesplit)
PHI_BAND = 1.25                # split angle for boxes with l1 >= 0.9 whose axis-crossing angle pi/2 - w meets [0.6, 1.05]
T_MIN = 0.5
J, M0, ROUNDS, CAP = 6, 64, 16, 8192
ROUNDS1, CAP1 = 6, 1024        # first (cheap) pass; ROUNDS/CAP only for small boxes (normalised width < DEEP_W)
DEEP_W = 2e-3
DEEP_ALWAYS = False            # set True by re-certification tools: always try ROUNDS/CAP after a failed cheap pass
LIMIT_P = 11
WIND_M0, WIND_ROUNDS, WIND_CAP = 64, 8, 16384
ORDER = 'wind_first'
SCALE = (1.0, 1.0, math.pi / 2, 1.0)


def dn(x):
    return np.nextafter(np.nextafter(x, -np.inf), -np.inf)


def up(x):
    return np.nextafter(np.nextafter(x, np.inf), np.inf)


def init(scale=None, face=False, rule=None, peri=None, faces=None, tx=None, c3prune=False):
    global cnt, away, wbc, wcl, wbox, SCALE, FACE_MODE, SPLIT_RULE, PERI, cntp, awarc, FACE_S, T_X, C3PRUNE
    if scale is not None:
        SCALE = tuple(scale)
    FACE_MODE = face
    SPLIT_RULE = rule
    PERI = peri
    if faces is not None:
        FACE_S = faces
    if tx is not None:
        T_X = tx
    C3PRUNE = c3prune
    import importlib
    import nb_cert2
    import wind_cert
    cnt = importlib.import_module('nbeng_c3_KBQe_2').count       # second-order enclosure (c3/make_c3_2.py)
    away = nb_cert2.away_count_t2
    awarc = nb_cert2.away_count_arc2
    cntp = importlib.import_module('nbeng_c3_KBPe_2').count     # perihelion polynomial (kq_run --peri)
    wbc = nb_cert2.box_coeffs2; wcl = wind_cert.choose_L; wbox = nb_cert2.wind_box2
    certify((0.96, 0.961, 0.8, 0.801, 0.5, 0.501, 0.9, 0.901))
    wind((0.96, 0.961, 0.8, 0.801, 0.5, 0.501, 0.9, 0.901))


def wind(bx):
    """winding certificate of the f-chart polynomial (chart 1 of nb_cert) on the box: (ok, W, L)"""
    cm = wbc(*bx[:6])
    L = wcl(cm[0], 0.5 * (bx[6] + bx[7]))
    if L is None:
        return False, -1, None
    ok, W = wbox(*cm, bx[6], bx[7], L, WIND_M0, WIND_ROUNDS, WIND_CAP)
    return ok and W >= 1, W, L


def phi_split(bx):
    """split angle phi0 between the near part |f - pi| <= phi0 and the away part, a fixed function of the box"""
    if bx[3] >= 0.9 and (math.pi / 2 - bx[4]) > 0.6 and (math.pi / 2 - bx[5]) < 1.05:
        return PHI_BAND
    return PHI_0


def count_peri(bx, rounds, cap, phi0, peri, t_lo=None):
    """--peri count of the box: f-chart arcs on [pi + phi0, 2 pi - peri] and [peri, pi - phi0], the perihelion
    polynomial (chart KBPe, tp = tan(f/2)) on |tp| <= tan(peri/2), and the aphelion polynomial (KBQe) on
    |th| <= tan(phi0/2), or only on t_lo <= |th| <= tan(phi0/2) (both signs) when t_lo is given (the outer count of
    --face).  The four parts cover the circle with overlaps; returns the sum of the four bounds."""
    t_hi = math.tan(phi0 / 2) + 1e-9
    tp_hi = math.tan(peri / 2) + 1e-9
    e0, e1, l0, l1, w0, w1, t0, t1 = bx
    pz0 = dn(l0 * dn(math.sin(w0))); pz1 = up(l1 * up(math.sin(w1)))
    qz0 = dn(l0 * dn(math.cos(w1))); qz1 = up(l1 * up(math.cos(w0)))
    ep0 = max(dn(1.0 - e1), 0.0); ep1 = up(1.0 - e0)
    t1c = min(t1, 1.0)
    mu0 = max(dn(math.sqrt(max((1.0 - t1c) / t1c, 0.0))), 0.0); mu1 = up(math.sqrt((1.0 - t0) / t0))
    lim = LIMIT_OUTER if t_lo is not None else LIMIT_P
    arcl = math.pi - phi0 - peri + 2e-9
    tot = 0.0
    for stage_target in (2.0, 0.0):
        if phi0 + peri >= math.pi:       # the two polynomial parts already cover the circle: no f-chart arc
            a1 = a2 = 0.0
        else:
            ok1, a1 = awarc(e0, e1, l0, l1, w0, w1, t0, t1, math.pi + phi0 - 1e-9, arcl, 16, rounds, J, stage_target, cap)
            ok2, a2 = awarc(e0, e1, l0, l1, w0, w1, t0, t1, peri - 1e-9, arcl, 16, rounds, J, stage_target, cap)
        if a1 + a2 > lim:
            if stage_target == 0.0:
                return False, 999
            continue
        if t_lo is None:
            okq, q = cnt(ep0, ep1, pz0, pz1, qz0, qz1, mu0, mu1, -t_hi, t_hi, J, M0, rounds, cap, lim - a1 - a2)
        else:
            okq1, q1 = cnt(ep0, ep1, pz0, pz1, qz0, qz1, mu0, mu1, t_lo, t_hi, J, M0, rounds, cap, lim - a1 - a2)
            okq2, q2 = cnt(ep0, ep1, pz0, pz1, qz0, qz1, mu0, mu1, -t_hi, -t_lo, J, M0, rounds, cap,
                           lim - a1 - a2 - q1)
            q = q1 + q2
        if a1 + a2 + q > lim:
            if stage_target == 0.0:
                return False, int(min(a1 + a2 + q, 999))
            continue
        okp, qp = cntp(ep0, ep1, pz0, pz1, qz0, qz1, mu0, mu1, -tp_hi, tp_hi, J, M0, rounds, cap, lim - a1 - a2 - q)
        tot = a1 + a2 + q + qp
        if tot <= lim or stage_target == 0.0:
            return tot <= lim, int(min(tot, 999))
    return False, int(min(tot, 999))


def count_total(bx, rounds, cap, phi0=PHI_0, peri=None):
    """the KQ count (away + near) of the box with the given refinement rounds, arc cap and split angle:
    (ok, total); with peri, the count of count_peri"""
    if peri is not None:
        return count_peri(bx, rounds, cap, phi0, peri)
    t_hi = math.tan(phi0 / 2) + 1e-9
    e0, e1, l0, l1, w0, w1, t0, t1 = bx
    pz0 = dn(l0 * dn(math.sin(w0))); pz1 = up(l1 * up(math.sin(w1)))
    qz0 = dn(l0 * dn(math.cos(w1))); qz1 = up(l1 * up(math.cos(w0)))
    ep0 = max(dn(1.0 - e1), 0.0); ep1 = up(1.0 - e0)
    t1c = min(t1, 1.0)
    mu0 = dn(math.sqrt(max((1.0 - t1c) / t1c, 0.0))); mu1 = up(math.sqrt((1.0 - t0) / t0))
    mu0 = max(mu0, 0.0)
    ok, aw = away(e0, e1, l0, l1, w0, w1, t0, t1, phi0 - 1e-9, 64, rounds, J, 6.0, cap)
    if aw > LIMIT_P:
        return False, 999
    budget = LIMIT_P - aw
    okq, q = cnt(ep0, ep1, pz0, pz1, qz0, qz1, mu0, mu1, -t_hi, t_hi, J, M0, rounds, cap, budget)
    if q > budget:
        # second stage: tighten the away count and retry
        ok, aw = away(e0, e1, l0, l1, w0, w1, t0, t1, phi0 - 1e-9, 64, rounds, J, 0.0, cap)
        budget = LIMIT_P - aw
        if budget < 0:
            return False, 999
        okq, q = cnt(ep0, ep1, pz0, pz1, qz0, qz1, mu0, mu1, -t_hi, t_hi, J, M0, rounds, cap, budget)
    tot = aw + q
    return tot <= LIMIT_P, int(min(tot, 999))


def in_face(bx):
    """True if every point of the box has eps <= FACE_S^2 and mu <= FACE_S^2 (mu = sqrt((1 - t)/t) is decreasing in t)"""
    eps_hi = 1.0 - bx[0]
    mu_hi = math.sqrt(max((1.0 - bx[6]) / bx[6], 0.0))
    lim = FACE_S * FACE_S * (1 - 1e-9)
    return eps_hi <= lim and mu_hi <= lim


def in_c3(bx):
    """True if every point of the box lies in the C3 cell directions l >= C3_L, w >= C3_W"""
    return bx[2] >= C3_L and bx[4] >= C3_W


def count_outer(bx, rounds, cap):
    """--face mode: zeros on T_X <= |th| <= tan(PHI_0/2) (near polynomial, both signs) plus the away arc |f - pi| >= PHI_0"""
    e0, e1, l0, l1, w0, w1, t0, t1 = bx
    pz0 = dn(l0 * dn(math.sin(w0))); pz1 = up(l1 * up(math.sin(w1)))
    qz0 = dn(l0 * dn(math.cos(w1))); qz1 = up(l1 * up(math.cos(w0)))
    ep0 = max(dn(1.0 - e1), 0.0); ep1 = up(1.0 - e0)
    t1c = min(t1, 1.0)
    mu0 = max(dn(math.sqrt(max((1.0 - t1c) / t1c, 0.0))), 0.0); mu1 = up(math.sqrt((1.0 - t0) / t0))
    ok, aw = away(e0, e1, l0, l1, w0, w1, t0, t1, PHI_0 - 1e-9, 64, rounds, J, 0.0, cap)
    if aw > LIMIT_OUTER:
        return False, 999
    ok1, q1 = cnt(ep0, ep1, pz0, pz1, qz0, qz1, mu0, mu1, T_X, T_HI, J, M0, rounds, cap, LIMIT_OUTER - aw)
    if aw + q1 > LIMIT_OUTER:
        return False, int(min(aw + q1, 999))
    ok2, q2 = cnt(ep0, ep1, pz0, pz1, qz0, qz1, mu0, mu1, -T_HI, -T_X, J, M0, rounds, cap, LIMIT_OUTER - aw - q1)
    tot = aw + q1 + q2
    return tot <= LIMIT_OUTER, int(min(tot, 999))


def certify_outer(bx):
    peris = [PERI] if PERI is None else [PERI, None]
    if DEEP_ALWAYS:
        peris += [p_ for p_ in (None, PHI_P) if p_ not in peris]
    best = (False, 999)
    for peri in peris:
        if peri is None:
            ok, tot = count_outer(bx, ROUNDS1, CAP1)
        else:
            ok, tot = count_peri(bx, ROUNDS1, CAP1, PHI_0, peri, T_X)
        if ok:
            return ok, tot
        if DEEP_ALWAYS or max((bx[2 * i + 1] - bx[2 * i]) / SCALE[i] for i in range(4)) < DEEP_W:
            if peri is None:
                ok2, tot2 = count_outer(bx, ROUNDS, CAP)
            else:
                ok2, tot2 = count_peri(bx, ROUNDS, CAP, PHI_0, peri, T_X)
            if ok2:
                return ok2, tot2
            tot = min(tot, tot2)
        if tot < best[1]:
            best = (False, tot)
    return best


def certify(bx):
    """KQ count: first with ROUNDS1/CAP1 (cheap), then, for boxes with all normalised widths below DEEP_W (or always
    when DEEP_ALWAYS, as set by the re-certification tools), with ROUNDS/CAP; split angle phi_split(bx), then the
    default PHI_0; with --peri first the perihelion count (count_peri), then the plain count; in re-certification
    mode every option.  Any of these counts is a valid certificate."""
    if bx[6] < T_MIN:
        return False, 999
    phis = [phi_split(bx)]
    if phis[0] != PHI_0:             # always also the default split angle (any split angle gives a valid count)
        phis.append(PHI_0)
    peris = [PERI] if PERI is None else [PERI, None]   # --peri: the perihelion count, then the plain count
    if DEEP_ALWAYS:                  # re-certification: also the other perihelion option (records of either kind)
        peris += [p_ for p_ in (None, PHI_P) if p_ not in peris]
    best = (False, 999)
    for peri in peris:
        for phi0 in phis:
            ok, tot = count_total(bx, ROUNDS1, CAP1, phi0, peri)
            if ok:
                return ok, tot
            if DEEP_ALWAYS or max(splitrule.rel_widths(bx, SCALE, SPLIT_RULE)) < DEEP_W:
                ok2, tot2 = count_total(bx, ROUNDS, CAP, phi0, peri)
                if ok2:
                    return ok2, tot2
                tot = min(tot, tot2)
            if tot < best[1]:
                best = (False, tot)
    return best


def split(bx):
    if SPLIT_RULE is not None:
        return splitrule.split_box(bx, SCALE, SPLIT_RULE)
    w = [(bx[2 * i + 1] - bx[2 * i]) / SCALE[i] for i in range(4)]
    i = int(np.argmax(w))
    mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx)
    a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)


def work(args):
    bx0, tmax, minw = args
    t0 = time.time()
    stack = [tuple(bx0)]; code = []; nleaf = 0; hist = {}; fails = []; wL = []; nwind = 0
    while stack:
        if time.time() - t0 > tmax:
            return dict(box=list(bx0), done=False, tree=''.join(code), nleaf=nleaf, nwind=nwind, wL=wL, hist=hist,
                        fails=fails, left=[list(b) for b in stack])
        bx = stack.pop()
        if FACE_MODE:
            if in_c3(bx):
                code.append('o'); nleaf += 1; hist['o'] = hist.get('o', 0) + 1
                continue
            good, tot = certify_outer(bx)
            if good:
                code.append(chr(ord('a') + tot)); nleaf += 1; hist[tot] = hist.get(tot, 0) + 1
                continue
            if min((bx[2 * i + 1] - bx[2 * i]) / SCALE[i] for i in range(4)) < minw:
                code.append('!'); fails.append((list(bx), tot))
                continue
            code.append('+')
            stack.extend(split(bx))
            continue
        if C3PRUNE and in_c3(bx):
            code.append('o'); nleaf += 1; hist['o'] = hist.get('o', 0) + 1
            continue
        if in_face(bx):
            code.append('p'); nleaf += 1; hist['p'] = hist.get('p', 0) + 1
            continue
        if ORDER == 'wind_first':
            okw, W, L = wind(bx)
            if okw:
                code.append(chr(ord('0') + W)); wL.append(L); nleaf += 1; nwind += 1
                hist['w%d' % W] = hist.get('w%d' % W, 0) + 1
                continue
        good, tot = certify(bx)
        if good:
            code.append(chr(ord('a') + tot)); nleaf += 1; hist[tot] = hist.get(tot, 0) + 1
            continue
        if ORDER != 'wind_first':
            okw, W, L = wind(bx)
            if okw:
                code.append(chr(ord('0') + W)); wL.append(L); nleaf += 1; nwind += 1
                hist['w%d' % W] = hist.get('w%d' % W, 0) + 1
                continue
        if min(splitrule.rel_widths(bx, SCALE, SPLIT_RULE)) < minw:
            code.append('!'); fails.append((list(bx), tot))
            continue
        code.append('+')
        stack.extend(split(bx))
    return dict(box=list(bx0), done=True, tree=''.join(code), nleaf=nleaf, nwind=nwind, wL=wL, hist=hist, fails=fails,
                left=[])


if __name__ == '__main__':
    resume = '--resume' in sys.argv
    FACE_MODE = '--face' in sys.argv
    C3PRUNE = '--c3prune' in sys.argv
    for a_ in sys.argv:
        if a_.startswith('--faces='):
            FACE_S = float(a_.split('=', 1)[1]); T_X = FACE_Y * FACE_S
    for a_ in sys.argv:
        if a_.startswith('--tx='):
            T_X = float(a_.split('=', 1)[1])
    SPLIT_RULE = 'face' if '--facesplit' in sys.argv else None
    for a_ in sys.argv:
        if a_.startswith('--splitrule='):
            SPLIT_RULE = a_.split('=', 1)[1]
    PERI = PHI_P if '--peri' in sys.argv else None
    minw = 2e-5
    for a_ in sys.argv:            # --minw=X: smaller-width threshold for a failed leaf (continuation levels whose scale
        if a_.startswith('--minw='):   # differs from the previous level's, so inherited thin boxes are not failed at once)
            minw = float(a_.split('=', 1)[1])
    if sys.argv[1] == '--boxes':
        bfile = sys.argv[2]
        SCALE = tuple(float(v) for v in sys.argv[4:8])
        tmax = float(sys.argv[8]); workers = int(sys.argv[9]); out = sys.argv[10]
        boxes = [tuple(float(v) for v in b) for b in json.load(open(bfile))['boxes']]
        dom = None; n = None
    else:
        bfile = None
        dom = tuple(float(v) for v in sys.argv[1:9])
        n = [int(v) for v in sys.argv[9:13]]
        tmax = float(sys.argv[13]); workers = int(sys.argv[14]); out = sys.argv[15]
        edges = [np.linspace(dom[2 * i], dom[2 * i + 1], n[i] + 1) for i in range(4)]
        boxes = [tuple(float(v) for v in (edges[0][p], edges[0][p + 1], edges[1][q], edges[1][q + 1], edges[2][r],
                                          edges[2][r + 1], edges[3][s], edges[3][s + 1]))
                 for p in range(n[0]) for q in range(n[1]) for r in range(n[2]) for s in range(n[3])]
        if '--scale' in sys.argv:                 # optional split scale for a grid run (recorded in the header)
            k_ = sys.argv.index('--scale')
            SCALE = tuple(float(v) for v in sys.argv[k_ + 1:k_ + 5])
    header = dict(engine='KQF' if FACE_MODE else 'KQ', FACE_S=FACE_S, FACE_Y=FACE_Y, T_X=T_X, LIMIT_OUTER=LIMIT_OUTER,
                  domain=dom, n0=n, boxes_file=bfile, split_rule=SPLIT_RULE, peri=PERI, c3prune=C3PRUNE, exclude=[], skipped=[], PHI_0=PHI_0, T_HI=T_HI, PHI_BAND=PHI_BAND, phi_rule='l1 >= 0.9 and pi/2 - w0 > 0.6 and pi/2 - w1 < 1.05 -> PHI_BAND', T_MIN=T_MIN, J=J, M0=M0,
                  rounds=ROUNDS, cap=CAP, rounds1=ROUNDS1, cap1=CAP1, deep_w=DEEP_W, limit=LIMIT_P, minw=minw, scale=SCALE, wind_m0=WIND_M0,
                  wind_rounds=WIND_ROUNDS, wind_cap=WIND_CAP, order=ORDER, enclosure='nb_cert2 (second order)',
                  argv=sys.argv)
    res = []
    import os
    if resume and os.path.exists(out + '.records.jsonl'):
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
    frec = open(out + '.records.jsonl', 'a' if resume else 'w')
    t0 = time.time()
    with Pool(workers, initializer=init, initargs=(SCALE, FACE_MODE, SPLIT_RULE, PERI, FACE_S, T_X, C3PRUNE)) as pool:
        for k, r in enumerate(pool.imap_unordered(work, [(b, tmax, minw) for b in boxes])):
            res.append(r); frec.write(json.dumps(r) + chr(10)); frec.flush()
            if (k + 1) % 16 == 0 or k + 1 == len(boxes):
                print(f"[{time.time()-t0:7.0f}s] {k+1}/{len(boxes)} returned, {sum(x['done'] for x in res)} complete, "
                      f"{sum(x['nleaf'] for x in res)} leaves, {sum(len(x['fails']) for x in res)} failed", flush=True)
    frec.close()
    json.dump(dict(header=header, results=res), open(out, 'w'))
    H = {}
    for r in res:
        for kk, v in r['hist'].items():
            H[str(kk)] = H.get(str(kk), 0) + v
    print('done KQ incomplete', sum(not x['done'] for x in res), 'failed', sum(len(x['fails']) for x in res),
          'hist', dict(sorted(H.items())), 'wall', round(time.time() - t0), flush=True)

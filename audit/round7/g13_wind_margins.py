"""AUDIT-7 g13 (final pass; AUDIT-4 d5 with gzip/records input and engine KQ = f-chart, chart 1 of nb_cert):
for EVERY winding leaf of the given runs, re-run wind_box's refinement with the numba labels(), then
evaluate the winding steps in python exactly as wind_box does and record
  * max over arcs of |step| (a step within ~1e-15 of pi is the branch-flip risk of d4),
  * the minimum over nodes of |labelled component| / |G| (how close node values come to the label boundary),
  * the wrap-around discrepancy |G(left0 + 2pi) - G(left0)| / |G(left0)|,
  * W agreement with the recorded digit, L < 0, max |k L|.
usage: python g13_wind_margins.py NMAX RUN.json[.gz|.records.jsonl] [...]   (.venv python; NMAX < 0 = all winding leaves)"""
import os, sys, json, gzip, math, random, importlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.join(HERE, '..', '..', 'code'); sys.path.insert(0, CODE)
import wind_cert
import nb_cert2
import splitrule
NMAX = int(sys.argv[1]); files = sys.argv[2:]
random.seed(3)


SPLIT_RULE = None                                    # header split_rule of the file being replayed (round 9)


def split(bx, S):
    if SPLIT_RULE:
        return splitrule.split_box(bx, S, SPLIT_RULE)
    w = [(bx[2 * i + 1] - bx[2 * i]) / S[i] for i in range(4)]
    i = int(np.argmax(w)); mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)


def steps(cm, t0, t1, L, M0, rounds, cap):
    second = len(cm) == 8                       # nb_cert2 data (cmid, crad, dcen, dwid, cmax, rr, hcen, rem)
    cmid, crad, dcen, dwid, cmax, rr = cm[:6]
    left = np.arange(M0) * (2 * math.pi / M0); width = np.full(M0, 2 * math.pi / M0)
    for rd in range(rounds + 1):
        if second:
            lab = nb_cert2.labels2(cmid, crad, dcen, dwid, cmax, rr, cm[6], cm[7], t0, t1, L, left + width / 2, width / 2 + wind_cert.ARC_PAD)
        else:
            lab = wind_cert.labels(cmid, crad, dcen, dwid, cmax, rr, t0, t1, L, left + width / 2, width / 2 + wind_cert.ARC_PAD)
        nbad = int(np.sum(lab == 0))
        if nbad == 0:
            M = left.size
            vals = [wind_cert.G_centre(cmid, t0, L, left[i]) for i in range(M)] + [wind_cert.G_centre(cmid, t0, L, left[0] + 2 * math.pi)]
            tot = 0.0; maxstep = 0.0; minrel = 1.0
            for i in range(M):
                (re0, im0), (re1, im1) = vals[i], vals[i + 1]
                lb = lab[i]
                for rv, iv in ((re0, im0), (re1, im1)):
                    comp = {1: rv, 2: iv, 3: -rv, 4: -iv}[int(lb)]
                    if not comp > 0:
                        return None
                    minrel = min(minrel, comp / math.hypot(rv, iv))
                d = math.atan2(im1, re1) - math.atan2(im0, re0)
                while d > math.pi: d -= 2 * math.pi
                while d < -math.pi: d += 2 * math.pi
                tot += d; maxstep = max(maxstep, abs(d))
            wrap = math.hypot(vals[-1][0] - vals[0][0], vals[-1][1] - vals[0][1]) / math.hypot(*vals[0])
            w = int(round(tot / (2 * math.pi)))
            return 6 + w, maxstep, minrel, wrap, abs(tot - 2 * math.pi * w)
        if rd == rounds or left.size > cap:
            return None
        nl = []; nw = []
        for i in range(left.size):
            if lab[i] == 0:
                nl += [left[i], left[i] + width[i] / 2]; nw += [width[i] / 2] * 2
            else:
                nl.append(left[i]); nw.append(width[i])
        left = np.array(nl); width = np.array(nw)
    return None


def load_any(fn):
    if fn.endswith('.records.jsonl'):                 # partial KQ / bulk file: completed records, driver scale
        recs = [r for r in (json.loads(l) for l in open(fn) if l.strip()) if r.get('done')]
        eng = 'KQ' if os.path.basename(fn).startswith('KQ') else 'a'
        S = [1.0, 1.0, math.pi / 2, 1.0] if eng == 'KQ' else [math.pi / 2, 1.0, math.pi / 2, 1.0]
        hdr = dict(engine=eng, scale=S)
        if eng == 'KQ':                                # kq_run since 09:09 (the final JSON header records them)
            hdr.update(wind_m0=64, wind_rounds=8, wind_cap=16384)
        import os as _os
        _r = _os.environ.get('AUDIT_SPLIT_RULE'); _s = _os.environ.get('AUDIT_SCALE')
        if _r is None and _s is None:
            print('   WARNING partial file without header: leaves rebuilt with the default rule and driver scale; '
                  'valid only for runs without a split_rule and with the driver scale (set AUDIT_SPLIT_RULE / AUDIT_SCALE)', flush=True)
        if _r:
            hdr['split_rule'] = _r
        if _s:
            hdr['scale'] = [float(v) for v in _s.split(',')]
        return dict(header=hdr, results=recs)
    p = fn if os.path.exists(fn) else fn + '.gz'
    return json.load(gzip.open(p, 'rt') if p.endswith('.gz') else open(p))


for fn in files:
    if not (os.path.exists(fn) or os.path.exists(fn + '.gz')):
        print(os.path.basename(fn), 'not present', flush=True); continue
    D = load_any(fn); H = D['header']; eng = H['engine']; S = H['scale']
    SPLIT_RULE = H.get('split_rule')
    if eng == 'KQF':                          # the face outer count has no winding leaves
        print(os.path.basename(fn), 'engine KQF: no winding leaves', flush=True); continue
    if eng in ('a', 'f', 'KQ'):
        import nb_cert
        ch = 0 if eng == 'a' else 1
        bc = lambda bx: nb_cert.box_coeffs(ch, *bx[:6])
    elif eng in ('C1', 'C1p', 'C1pp'):
        m = importlib.import_module('nbeng_' + eng); bc = lambda bx, m=m: m.box_coeffs(*bx[:6])
    else:
        m = importlib.import_module('nbeng_' + eng); bc = lambda bx, m=m: m.box_coeffs(*bx)
    leaves = []
    for r in D['results']:
        stack = [tuple(r['box'])]; wi = 0
        for c in r['tree']:
            bx = stack.pop()
            if c == '+':
                stack.extend(split(bx, S))
            elif '1' <= c <= '6':
                leaves.append((bx, int(c), r['wL'][wi])); wi += 1
    sample = leaves if NMAX < 0 or NMAX >= len(leaves) else random.sample(leaves, NMAX)
    worst_step = 0.0; worst_rel = 1.0; worst_wrap = 0.0; worst_res = 0.0; mism = 0; fail = 0; maxkL = 0.0; Lpos = 0; n_second = 0
    for bx, W, L in sample:
        if not L < 0:
            Lpos += 1
        maxkL = max(maxkL, 6 * abs(L))
        res = steps(bc(bx), bx[6], bx[7], L, H.get('wind_m0', 64), H.get('wind_rounds', 3), H.get('wind_cap', 1024))
        if res is None and eng in ('a', 'f', 'KQ'):        # leaf certified by the second-order enclosure (nb_cert2)
            res = steps(nb_cert2.box_data2(0 if eng == 'a' else 1, *bx[:6]), bx[6], bx[7], L, H.get('wind_m0', 64),
                        H.get('wind_rounds', 3), H.get('wind_cap', 1024))
            n_second += res is not None
        if res is None:
            fail += 1; continue
        W2, ms, mr, wr, rsd = res
        mism += W2 != W
        worst_step = max(worst_step, ms); worst_rel = min(worst_rel, mr); worst_wrap = max(worst_wrap, wr); worst_res = max(worst_res, rsd)
    print(f"{os.path.basename(fn)}: winding leaves {len(leaves)}, checked {len(sample)}; W mismatches {mism}; re-run failures {fail}; "
          f"L >= 0: {Lpos}; max |kL| {maxkL:.1f}; max |step| = pi - {math.pi - worst_step:.3e}; "
          f"min node |label component|/|G| {worst_rel:.3e}; max wrap discrepancy {worst_wrap:.1e}; max |tot - 2 pi w| {worst_res:.1e}; "
          f"labelled only by the second-order enclosure {n_second}",
          flush=True)

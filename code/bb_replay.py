"""Replay the split trees recorded by bb_run.py and check coverage; optionally re-certify leaves.

usage: python bb_replay.py OUT.json [--recert N] [--seed S]
  Replays every record from its initial box with the recorded split rule and tree code and checks:
    * the code is consumed exactly, the final stack is empty (done) or equals the stored 'left' list;
    * failed boxes are exactly the stored 'fails' (in order);
    * leaf + fail + left volumes add up to the initial-box volume (relative 1e-9);
    * the initial boxes are the grid of the header domain minus the excluded boxes (if the run used a grid).
  --recert N: re-certifies N random leaves (all leaves if N < 0) with certify_box_deep (regular-certified leaves too)
  and checks that each is certified with a bound <= the recorded bound.
"""
import argparse
import json
from runio import load_json
from splitrule import split_box
import math
import random
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('file')
ap.add_argument('--recert', type=int, default=0)
ap.add_argument('--seed', type=int, default=0)
A = ap.parse_args()
D = load_json(A.file)
H = D['header']
SCALE = H['scale']


# engines whose leaf bound is a whole-circle total: a recorded 11 ('l'/'L') is valid only there (parity lemma)
PARITY_ENGINES = ('a', 'f', 'C1', 'C1p', 'C1pp', 'B1', 'B2', 'B3p', 'B3m', 'L2E', 'L2P', 'L2Ap', 'L2Am', 'L2Tp', 'L2Tm',
                  'D2E', 'D2E2', 'D2P', 'D2K', 'D2Tp', 'D2Tm', 'KQ')       # (KQF counts part of the circle: not here)


def split(bx):
    if H.get('split_rule') is not None:          # kq_run --facesplit (splitrule.py)
        return split_box(bx, SCALE, H['split_rule'])
    w = [(bx[2 * i + 1] - bx[2 * i]) / SCALE[i] for i in range(4)]
    i = int(np.argmax(w))
    mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx)
    a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)


def vol(b):
    return float(np.prod([b[2 * i + 1] - b[2 * i] for i in range(4)]))


leaves = []
bad = 0
nfail = 0
nleft = 0
hist = {}
for r in D['results']:
    root = tuple(r['box'])
    stack = [root]
    fl = [tuple(f[0]) for f in r['fails']]
    fi = 0
    wi = 0
    v = 0.0
    nl = 0
    for ch in r['tree']:
        bx = stack.pop()
        if ch == '+':
            stack.extend(split(bx))
        elif ch == '!':
            if fi >= len(fl) or fl[fi] != bx:
                bad += 1
                print('fail mismatch in', root)
            fi += 1
            v += vol(bx)
        else:
            Lw = None
            if 'a' <= ch <= 'k':
                bound, method = ord(ch) - ord('a'), 0
            elif 'A' <= ch <= 'K':
                bound, method = ord(ch) - ord('A'), 1
            elif ch in 'lL' and H['engine'] in PARITY_ENGINES:
                bound, method = 11, (0 if ch == 'l' else 1)   # parity lemma: total 11 -> at most 10 (whole-circle engines)
            elif ch in 'uop':                     # pruned: unphysical (D2), C3 cell ('o', KQF), face neighbourhood ('p', KQ)
                v += vol(bx); nl += 1
                continue
            elif '1' <= ch <= '6':
                bound, method = 12 - 2 * (ord(ch) - ord('0')), 2
                Lw = r['wL'][wi]; wi += 1
            else:
                raise SystemExit('bad code ' + ch)
            leaves.append((bx, bound, method, Lw))
            hist[bound] = hist.get(bound, 0) + 1
            nl += 1
            v += vol(bx)
    if fi != len(fl):
        bad += 1
        print('unconsumed fails in', root)
    left = [tuple(b) for b in r['left']]
    if r['done']:
        if stack:
            bad += 1
            print('non-empty stack in completed record', root)
    else:
        if stack != left:
            bad += 1
            print('left mismatch in', root)
    v += sum(vol(b) for b in left)
    if nl != r['nleaf']:
        bad += 1
        print('leaf count mismatch', root)
    if abs(v - vol(root)) > 1e-9 * vol(root):
        bad += 1
        print('volume mismatch', root, v, vol(root))
    nfail += len(fl)
    nleft += len(left)

# initial boxes against the header grid
if H['domain'] is not None:
    dom = H['domain']; n0 = H['n0']
    edges = [np.linspace(dom[2 * i], dom[2 * i + 1], n0[i] + 1) for i in range(4)]
    grid = set(tuple(float(x) for x in (edges[0][a], edges[0][a + 1], edges[1][b], edges[1][b + 1], edges[2][c],
                                        edges[2][c + 1], edges[3][d], edges[3][d + 1]))
               for a in range(n0[0]) for b in range(n0[1]) for c in range(n0[2]) for d in range(n0[3]))
    roots = set(tuple(r['box']) for r in D['results'])
    skipped = set(tuple(b) for b in H['skipped'])
    if roots | skipped != grid or roots & skipped:
        bad += 1
        print('initial boxes do not match the grid:', len(grid), len(roots), len(skipped))
print('records', len(D['results']), 'complete', sum(r['done'] for r in D['results']), 'leaves', len(leaves),
      'deep', sum(1 for l in leaves if l[2] == 1), 'wind', sum(1 for l in leaves if l[2] == 2), 'fails', nfail, 'left', nleft, 'bound histogram', dict(sorted(hist.items())),
      'replay problems', bad)

if A.recert:
    eng = H['engine']
    if eng in ('a', 'f'):
        # re-certified with the second-order engine (nb_cert2): its radius is the minimum of the mean-value bound of
        # nb_cert and the second-order bound, so every leaf of either engine version re-certifies
        import nb_cert, nb_cert2
        ch = 0 if eng == 'a' else 1
        cb = lambda bx: nb_cert2.certify_box_deep2(ch, *bx, H['J'], H['M0'], H['deep_rounds'], H['deep_cap'])
    elif eng == 'K':
        import nb_corner
        cb = lambda bx: nb_corner.certify_corner_box(*bx, H['phi0'], H['phi_in'])[:2]
    elif eng == 'KQ':
        import kq_run
        kq_run.init(H['scale'], False, None, None, H.get('FACE_S'), H.get('T_X'))
        kq_run.DEEP_ALWAYS = True          # re-certification: always allow the deep pass
        cb = kq_run.certify
    elif eng == 'KQF':
        import kq_run
        kq_run.init(H['scale'], True, None, None, H.get('FACE_S'), H.get('T_X'))   # the run's own outer start T_X
        kq_run.DEEP_ALWAYS = True
        cb = kq_run.certify_outer
    elif eng == 'KA':
        import c3_away
        cb = lambda bx: c3_away.away_count_t(bx[0], bx[1], bx[2], bx[3], bx[4], bx[5], 1.0 / (1.0 + bx[7] ** 2),
                                             1.0 / (1.0 + bx[6] ** 2), H['phi_a'], 64, 16, H['J'], 2.0, 8192)
    else:
        import importlib
        m = importlib.import_module('nbeng_' + eng)
        cb = lambda bx: m.certify_box_deep(*bx, H['J'], H['M0'], H['deep_rounds'], H['deep_cap'])
    random.seed(A.seed)
    sample = leaves if A.recert < 0 else random.sample(leaves, min(A.recert, len(leaves)))
    import wind_cert
    if eng in ('a', 'f'):
        bcf = lambda bx: nb_cert2.box_data2(ch, *bx[:6])
    elif eng in ('C1', 'C1p', 'C1pp'):
        bcf = lambda bx: m.box_coeffs(*bx[:6])
    elif eng in ('B1', 'B2', 'L2E', 'L2P', 'L2Ap', 'L2Am', 'D2E', 'D2E2', 'D2P', 'D2K'):
        bcf = lambda bx: m.box_coeffs(*bx)
    elif eng == 'KQ':
        import nb_cert2
        bcf = lambda bx: nb_cert2.box_coeffs2(*bx[:6])
    else:
        bcf = None
    worse = 0
    for bx, bound, method, Lw in sample:
        if method == 2:
            wbx = nb_cert2.wind_box2 if eng in ('KQ', 'a', 'f') else wind_cert.wind_box
            okw, W = wbx(*bcf(bx), bx[6], bx[7], Lw, H.get('wind_m0', 64), H.get('wind_rounds', 3),
                         H.get('wind_cap', 1024))
            if not okw or 12 - 2 * W > bound:
                worse += 1
                print('winding re-certification problem', bx, bound, okw, W)
            continue
        ok, tot = cb(bx)
        # whole-circle engines: the final engines stop as soon as the total is <= 11 (parity lemma), so a leaf recorded
        # with a smaller bound by a pre-parity run is re-certified when the total is <= 11 (audit R7-1)
        lim = 11 if eng in PARITY_ENGINES else bound
        if not ok or int(tot) > max(bound, lim):
            worse += 1
            print('re-certification problem', bx, bound, method, ok, tot)
    print('re-certified', len(sample), 'leaves; problems', worse)

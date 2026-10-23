"""Replay and cross-check all certificate runs, and write ../RESULTS.md.

A run family is a base JSON plus optional continuation JSONs (bb_run --boxes ...), named
  runs/NAME.json, runs/NAME.r1.json, runs/NAME.r2.json, ...
Every record is replayed exactly (split rule from the header scale, or the header split_rule: splitrule.py).  The family is COMPLETE when
  * every record of every file replays with no inconsistency,
  * every box left on the stack of an incomplete record, and every failed leaf ('!'), is the root of a record in a
    later continuation file of the same family, and that record is complete,
  * the base file's initial boxes are its grid minus the declared exclusions (bb_run) or its full grid (c3, kb).
Leaf codes: '+' split, '!' failed, 'a'..'k' / 'A'..'K' certified bound 0..10, 'l' / 'L' certified total 11 (whole-circle
engines only; at most 10 by the parity lemma), '1'..'6' winding certificate with W roots inside (bound 12 - 2W),
'u' 'o' 'p' pruned (outside the region), '7' tau-winding certificate (whole circle, <= 10; tau_wind.py).  Run files may be gzip-compressed (OUT.json.gz).
"""
import glob
import json
import math
import os
import sys
import numpy as np
from runio import exists, load_json
from splitrule import split_box
import tau_wind
assert tau_wind.DISC_W == 6, 'disc-winding leaves (code 8) certify the outer count 12 - W <= 6 only with W >= 6'

FAMILIES = [('bulk', 'runs/bulk'), ('C1', 'runs/C1'), ('C1p', 'runs/C1p'), ('C1pp', 'runs/C1pp'), ('B1', 'runs/B1'),
            ('B2', 'runs/B2'), ('B3p', 'runs/B3p'), ('B3m', 'runs/B3m'), ('L2E', 'runs/L2E'), ('L2P', 'runs/L2P'),
            ('L2Ap', 'runs/L2Ap'), ('L2Am', 'runs/L2Am'), ('L2Tp', 'runs/L2Tp'), ('L2Tm', 'runs/L2Tm'),
            ('F0', 'runs/F0'), ('F1', 'runs/F1'), ('KQ0', 'runs/KQ0'), ('KQ0T', 'runs/KQ0T'), ('KQ1', 'runs/KQ1'), ('KQ2', 'runs/KQ2'), ('KQB', 'runs/KQB'),
            ('KQK1', 'runs/KQK1'), ('KQK2', 'runs/KQK2'), ('KA', 'runs/KA'), ('KQF', 'runs/KQF'), ('KQF2', 'runs/KQF2'), ('KQS', 'runs/KQS'),
            ('D2E', 'runs/D2E'), ('D2E2', 'runs/D2E2'), ('D2P', 'runs/D2P'), ('D2K', 'runs/D2K'), ('D2Tp', 'runs/D2Tp'), ('D2Tm', 'runs/D2Tm')] + \
           [('C3 ' + c, 'c3/runs/' + c) for c in ('Ca', 'CK', 'CE', 'Cm', 'PhiP', 'PhiM', 'L2aE', 'L2am', 'L2aYp',
                                                  'L2aYm', 'L2KE', 'L2Km', 'L2KYp', 'L2KYm')] + \
           [('face ' + c, 'c3/runs/' + c) for c in ('FE', 'FM', 'FIP', 'FIM')]


def vol(b):
    return float(np.prod([b[2 * i + 1] - b[2 * i] for i in range(4)]))


PARITY_ENGINES = ('a', 'f', 'C1', 'C1p', 'C1pp', 'B1', 'B2', 'B3p', 'B3m', 'L2E', 'L2P', 'L2Ap', 'L2Am', 'L2Tp', 'L2Tm',
                  'D2E', 'D2E2', 'D2P', 'D2K', 'D2Tp', 'D2Tm', 'KQ')   # whole-circle totals: 11 allowed (parity lemma)


def replay(rec, scale, engine=None, rule=None, face_s=None, c3prune=False):
    """returns (ok, leaves_by_bound, pruned, fails, left, problems)"""
    stack = [tuple(rec['box'])]
    hist = {}; pruned = 0; fails = []; probs = []; nwind = [0]; ntw = [0]
    v = 0.0
    for ch in rec['tree']:
        if not stack:
            probs.append('stack underflow'); break
        bx = stack.pop()
        if ch == '+' and rule is not None:                  # header split_rule (kq_run --facesplit): splitrule.py
            stack.extend(split_box(bx, scale, rule))
            continue
        if ch == '+':
            w = [(bx[2 * i + 1] - bx[2 * i]) / scale[i] for i in range(4)]
            i = int(np.argmax(w))
            mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
            a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
            stack.extend([tuple(a), tuple(c)])
            continue
        v += vol(bx)
        if ch == '!':
            fails.append(bx)
        elif ch in 'uop':
            pruned += 1
            if not prune_ok(ch, engine, bx, face_s, c3prune):
                probs.append('pruning predicate fails for code ' + ch)
        elif 'a' <= ch <= 'k' or 'A' <= ch <= 'K':
            b = ord(ch) - (ord('a') if ch >= 'a' else ord('A'))
            hist[b] = hist.get(b, 0) + 1
        elif ch in 'lL' and engine in PARITY_ENGINES:  # total 11 -> at most 10 by the parity lemma
            hist[11] = hist.get(11, 0) + 1
        elif '1' <= ch <= '6':                     # winding certificate: W roots inside |z| < rho -> <= 12 - 2W real
            b = 12 - 2 * (ord(ch) - ord('0'))
            hist[b] = hist.get(b, 0) + 1
            nwind[0] += 1
        elif ch == '7' and engine in ('KQ', 'KQF', 'Ca', 'FE', 'FM'):  # tau-winding (tau_wind.py): <= 10 real
            hist[10] = hist.get(10, 0) + 1
            ntw[0] += 1
        elif ch == '8' and engine == 'KQF':        # disc winding (tau_wind.certify_disc, radius the pinned T_X):
            hist[6] = hist.get(6, 0) + 1           # outer count <= 12 - W <= 6
        else:
            probs.append('bad code ' + ch)
    left = [tuple(b) for b in rec.get('left', [])]
    if rec['done'] and stack:
        probs.append('non-empty stack in a complete record')
    if not rec['done'] and stack != left:
        probs.append('stored stack differs from replay')
    if [tuple(f[0]) for f in rec.get('fails', [])] != fails:
        probs.append('fails differ from replay')
    v += sum(vol(b) for b in left)
    if nwind[0] != len(rec.get('wL', [])):
        probs.append('number of winding leaves differs from the stored radii')
    if ntw[0] != len(rec.get('tw', [])):
        probs.append('number of tau-winding leaves differs from the stored contours')
    for tw_ in rec.get('tw', []):
        if not (len(tw_) == 4 and tw_[0] in ((0, 1) if engine in ('KQ', 'KQF') else (0,)) and tw_[2] > tw_[3] > 0):
            probs.append('invalid tau-winding contour (must lie in the open upper half-plane)')
    if abs(v - vol(rec['box'])) > 1e-9 * max(vol(rec['box']), 1e-300):
        probs.append('volume mismatch')
    return hist, pruned, fails, left, probs


# per-family limits on the recorded leaf bounds (audit R7-4).  Lower-case 'a'+c / upper-case 'A'+c codes; digits are
# winding leaves (bound 12 - 2W <= 10).  C3 charts: lower case = near count (target of c3/c3_run.py TARGET), upper case =
# level-2 zone leaves (outer count <= ZONE_OUTER = 5; the inner count <= 4 comes from the level-2 charts).
C3_TARGET = {'Ca': 9, 'CK': 8, 'CE': 8, 'Cm': 8, 'PhiP': 0, 'PhiM': 0, 'L2aE': 4, 'L2am': 4, 'L2KE': 4, 'L2Km': 4,
             'L2aYp': 0, 'L2aYm': 0, 'L2KYp': 0, 'L2KYm': 0, 'FE': 5, 'FM': 5, 'FIP': 0, 'FIM': 0}
# face neighbourhood (kq_run FACE_S): eps <= FACE_S^2 and mu <= FACE_S^2; C3 cell directions l >= 0.95, w >= 9 pi/20
FACE_S = 0.02
# corner K of the face construction (C3_CORE.md section 6): face charts only for s <= SF2 there, shell -> run KQS
SF2 = 0.010005
K_PZ, K_QZ0, K_QZ1 = 0.875, 0.125, 0.25
K_L, K_W = 0.88, 1.29
FACE_DEPS = ('c3/runs/FE', 'c3/runs/FM', 'c3/runs/FIP', 'c3/runs/FIM', 'runs/KQF', 'runs/KQF2', 'runs/KQS')
_FAMILY_CACHE = {}


def family_cached(prefix):
    if prefix not in _FAMILY_CACHE:
        _FAMILY_CACHE[prefix] = family(prefix)
    return _FAMILY_CACHE[prefix]
C3_L, C3_W = 0.95, 9 * math.pi / 20


def in_c3_cell(bx):
    """the whole KQ box (e, l, w, t) lies in the C3 cell e >= 0.99, mu <= 0.0505, l >= 0.95, w >= 9 pi/20"""
    return (bx[0] >= 0.99 and math.sqrt(max((1.0 - bx[6]) / bx[6], 0.0)) <= 0.0505 and bx[2] >= C3_L
            and bx[4] >= C3_W)


def prune_ok(code, engine, bx, face_s=None, c3prune=False):
    """verify the pruning predicate of a 'p' leaf of a KQ run (inside the face neighbourhood) or an 'o' leaf of a
    KQF run (inside the C3 cell directions); other pruning codes are not verified here (audited separately)."""
    if code == 'p' and engine == 'KQ':
        fs = FACE_S if face_s is None else face_s          # the run's face neighbourhood (header FACE_S)
        eps_hi = 1.0 - bx[0]; mu_hi = math.sqrt(max((1.0 - bx[6]) / bx[6], 0.0))
        return eps_hi <= fs ** 2 * (1 - 1e-9) and mu_hi <= fs ** 2 * (1 - 1e-9)
    if code == 'o' and engine == 'KQ' and c3prune:          # kq_run --c3prune
        return in_c3_cell(bx)
    if code == 'o' and engine == 'KQF':
        return in_c3_cell(bx)
    if engine in ('FE', 'FM', 'FIP', 'FIM'):   # face charts (c3_run.face_bounds): box (s, ratio, Pz, Qz)
        pz0, qz0, qz1 = bx[4], bx[6], bx[7]
        if code == 'u':
            return pz0 * pz0 + qz0 * qz0 > 1.0
        if code == 'o':
            l_lo = math.sqrt(pz0 * pz0 + qz0 * qz0) * (1 - 1e-12)
            w_lo = math.atan2(pz0, qz1) * (1 - 1e-12) if qz1 > 0 else math.pi / 2
            return (l_lo >= C3_L and w_lo >= C3_W) or bx[0] > FACE_S
        return False
    if code == 'p' or (code == 'o' and engine in ('KQ',)):
        return False
    return True
C3_ZONE_OUTER = 5


def code_limits(prefix, engine):
    """(max lower-case bound, max upper-case bound) allowed for the family"""
    name = os.path.basename(prefix)
    if '/c3/' in prefix.replace(chr(92), '/') or prefix.startswith('c3/'):
        return C3_TARGET[name], (C3_ZONE_OUTER if name in ('Ca', 'CK') else -1)
    if engine == 'KA':
        return 2, 2
    if engine == 'KQF':
        return 6, -1
    if engine in PARITY_ENGINES:
        return 11, 11
    return 10, 10


def family(prefix):
    conts = [s[:-3] if s.endswith('.gz') else s for s in glob.glob(prefix + '.r*.json') + glob.glob(prefix + '.r*.json.gz')]
    files = [prefix + '.json'] + sorted(set(conts), key=lambda s: int(s.split('.r')[-1].split('.')[0]))
    files = [f for f in files if exists(f)]
    if not files:
        return None
    hist = {}; pruned = 0; probs = []; owed = set(); roots_done = set(); nrec = 0; nfiles = len(files)
    incomplete_total = 0; fails_total = 0; uses_p = False
    for k, fn in enumerate(files):
        D = load_json(fn)
        H = D['header']
        scale = H['scale']
        lo_max, up_max = code_limits(prefix, H.get('engine'))
        fs_pin = None
        if H.get('engine') in ('KQ', 'KQF'):
            # the face radius and the outer start are pinned per family, not taken from the header (audit R11-1)
            fs_pin = SF2 if os.path.basename(prefix) in ('KQS', 'KQF2') else FACE_S
            if H.get('FACE_S') is not None and H['FACE_S'] != fs_pin:
                probs.append(fn + ': header FACE_S %r differs from the family value %r' % (H['FACE_S'], fs_pin))
            if H.get('engine') == 'KQF' and H.get('T_X') != 3.0 * fs_pin:
                probs.append(fn + ': header T_X %r differs from the family value %r' % (H.get('T_X'), 3.0 * fs_pin))
            if H.get('engine') == 'KQF' and any(c in r['tree'] for r in D['results'] for c in '123456'):
                probs.append(fn + ': winding leaf in a partial-count (KQF) family')
        uses_p = uses_p or any('p' in r['tree'] for r in D['results'])
        for rec in D['results']:
            nrec += 1
            for ch in set(rec['tree']):
                if 'a' <= ch <= 'l' and ord(ch) - ord('a') > lo_max:
                    probs.append(fn + ': leaf bound %d above the family limit %d' % (ord(ch) - ord('a'), lo_max))
                if 'A' <= ch <= 'L' and ord(ch) - ord('A') > up_max:
                    probs.append(fn + ': leaf code %s above the family limit %d' % (ch, up_max))
            h, p, fl, left, pr = replay(rec, scale, H.get('engine') or H.get('chart'), H.get('split_rule'),
                                        fs_pin, bool(H.get('c3prune')))
            for b, c in h.items():
                hist[b] = hist.get(b, 0) + c
            pruned += p
            probs += [fn + ': ' + x for x in pr]
            root = tuple(rec['box'])
            if k > 0:
                if root in owed:
                    owed.discard(root)
                else:
                    probs.append(fn + ': continuation root not owed by an earlier file')
            if rec['done']:
                roots_done.add(root)
            else:
                incomplete_total += 1
            for b in left:
                owed.add(b)
            for b in fl:
                owed.add(b)
            fails_total += len(fl)
        if k == 0 and os.path.basename(prefix) in ('KQS', 'KQF2'):
            # the corner runs must cover K_lw (l >= K_L, w >= K_W) of the face neighbourhood (KQS: eps, mu <= FACE_S^2
            # up to the 'p' tolerance of the other KQ runs; KQF2: eps, mu <= SF2^2)
            dom = H.get('domain') or [1, 0, 1, 0, 9, 0, 1, 0]
            rad = FACE_S if os.path.basename(prefix) == 'KQS' else SF2
            need_e0 = 1.0 - rad * rad * (1 - 1e-9); need_t0 = 1.0 / (1.0 + (rad * rad) ** 2)
            if not (dom[0] <= need_e0 and dom[1] >= 1.0 and dom[2] <= K_L and dom[3] >= 1.0 and dom[4] <= K_W
                    and dom[5] >= math.pi / 2 and dom[6] <= need_t0 and dom[7] >= 1.0):
                probs.append(fn + ': domain does not cover the corner region it must cover')
        if k == 0 and H.get('domain') is not None and H.get('n0') is not None:   # bb_run, kq_run and c3_run base files
            dom = H['domain']; n0 = H['n0']
            edges = [np.linspace(dom[2 * i], dom[2 * i + 1], n0[i] + 1) for i in range(4)]
            grid = set(tuple(float(x) for x in (edges[0][a], edges[0][a + 1], edges[1][b], edges[1][b + 1], edges[2][c],
                                                edges[2][c + 1], edges[3][d], edges[3][d + 1]))
                       for a in range(n0[0]) for b in range(n0[1]) for c in range(n0[2]) for d in range(n0[3]))
            roots = set(tuple(r['box']) for r in D['results'])
            skipped = set(tuple(b) for b in H.get('skipped', []))
            if roots | skipped != grid or roots & skipped:
                probs.append(fn + ': initial boxes do not match grid minus exclusions')
    # C1p: failed leaves at the D2b corner are covered by the D2 cone charts, which cover the cube
    # {e <= SD2, p <= SD2, kappa = 1 - q <= SD2, |tau - 1/2| <= SD2} (code/gen_bu.py D2*, runs/D2*).
    if prefix.endswith('C1p') and owed:
        SD2 = 5e-3
        keep = set()
        for b in owed:
            e_hi, l_lo, w_hi, t_lo, t_hi = b[1], b[2], b[5], b[6], b[7]
            p_hi = math.sin(w_hi); k_hi = 1 - l_lo * math.cos(w_hi)
            if not (e_hi <= SD2 and p_hi <= SD2 and k_hi <= SD2 and t_lo >= 0.5 - SD2 and t_hi <= 0.5 + SD2):
                keep.add(b)
        d2_ok = all(exists('runs/' + c + '.json') and (family('runs/' + c) or {}).get('complete')
                    for c in ('D2E', 'D2P', 'D2K', 'D2Tp', 'D2Tm'))
        if d2_ok:
            owed = keep
    # D2 charts: they are only needed on the failed (or unprocessed) leaves of C1p (rule above).  A failed or unprocessed
    # D2 box whose smallest physical e exceeds the largest e of every owed C1p box contains no point of those boxes, so
    # it is not owed.  Physical e = SD2 * nu (D2E) or SD2 * E * nu (D2P, D2K, D2T+-), nu = coordinate 0, E = coordinate 1.
    if os.path.basename(prefix).startswith('D2') and owed:
        SD2 = 5e-3
        C = load_json('runs/C1p.json')
        c1p_owed = [f[0] for r in C['results'] for f in r['fails']] + [b for r in C['results'] for b in r['left']]
        emax = max(b[1] for b in c1p_owed) if c1p_owed else -1.0
        def e_min(b):
            v = SD2 * b[0] if os.path.basename(prefix) in ('D2E', 'D2E2') else SD2 * b[0] * b[2]
            return v * (1 - 1e-12)
        owed = set(b for b in owed if not e_min(b) > emax)
    # D2Tp: a box (nu, E, P, K) whose image in the D2E variables (nu_E = E nu, P_E = P/E, K_E = K/E, T = tau'/e = 1/E)
    # lies in D2E's domain (nu_E <= 1.01, P_E, K_E <= 1.05, T in [-1.05, 1.05]) for T <= 1.05 and in D2E2's domain
    # (nu_E <= 0.14, P_E <= 0.35, K_E in [0.24, 0.68], T in [1.04, 1.28]) for T >= 1.04 is covered by those two charts
    # (same C1p function up to positive factors) when both families are complete.
    if os.path.basename(prefix) == 'D2Tp' and owed:
        d2e = family('runs/D2E') or {}; d2e2 = family('runs/D2E2') or {}
        if d2e.get('complete') and d2e2.get('complete'):
            def covered(b):
                n0, n1, E0, E1, P0, P1, K0, K1 = b
                if not E0 > 0:
                    return False
                nuE1 = E1 * n1 * (1 + 1e-12); PE1 = P1 / E0 * (1 + 1e-12)
                KE0 = K0 / E1 * (1 - 1e-12); KE1 = K1 / E0 * (1 + 1e-12)
                T0 = 1 / E1 * (1 - 1e-12); T1 = 1 / E0 * (1 + 1e-12)
                ok_e = nuE1 <= 1.01 and PE1 <= 1.05 and KE1 <= 1.05
                ok_e2 = nuE1 <= 0.14 and PE1 <= 0.35 and KE0 >= 0.24 and KE1 <= 0.68
                if T0 >= -1.05 and T1 <= 1.05:
                    return ok_e
                if T0 >= 1.04 and T1 <= 1.28:
                    return ok_e2
                return T0 >= -1.05 and T1 <= 1.28 and ok_e and ok_e2
            owed = set(b for b in owed if not covered(b))
    # face charts (c3/c3_run.py FE, FM, FIP, FIM; C3_CORE.md section 6): a point of the face neighbourhood has
    # s = max(sqrt(eps), sqrt(mu)) <= FACE_S (FE, FM), and a gap point tau in [YF s, T_X], T_X = YF FACE_S, has
    # sigma = tau/YF <= FACE_S (FIP, FIM).  The chart domains run to FACE_S R1 (R1 = 1.001); open or failed boxes with
    # s_lo > FACE_S contain no needed point.
    if os.path.basename(prefix) in ('FE', 'FM', 'FIP', 'FIM') and owed:
        owed = set(b for b in owed if not b[0] > FACE_S)
    # Corner K (C3_CORE.md section 6): for P_z >= 0.875, 0.125 <= Q_z <= 0.25 the face construction is used only for
    # s <= SF2 (outer count KQF2 from T_X2 = 3 SF2), and the shell SF2 < s <= FACE_S is covered by the whole-circle run
    # KQS.  So a face-chart box inside K with s_lo > SF2 is not needed, and neither is a KQF box inside
    # K_lw = {l >= 0.88, w >= 1.29} (which contains K, by l = |(P_z, Q_z)|, w = atan2(P_z, Q_z)).
    if os.path.basename(prefix) in ('FE', 'FM', 'FIP', 'FIM') and owed:
        owed = set(b for b in owed if not (b[4] >= K_PZ and b[6] >= K_QZ0 and b[7] <= K_QZ1 and b[0] > SF2))
    if os.path.basename(prefix) == 'KQF' and owed:
        owed = set(b for b in owed if not (b[2] >= K_L and b[4] >= K_W))
    if uses_p and not all((family_cached(fp) or {}).get('complete') for fp in FACE_DEPS if fp != prefix):
        probs.append('face-neighbourhood boxes pruned (p) but the face families are not all complete')
    complete = not probs and not owed
    return dict(files=nfiles, records=nrec, leaves=sum(hist.values()), pruned=pruned, hist=dict(sorted(hist.items())),
                fails=fails_total, incomplete=incomplete_total, owed=len(owed), problems=probs[:5],
                nproblems=len(probs), complete=complete)


if __name__ == '__main__':
    rows = []
    allok = True
    for name, prefix in FAMILIES:
        r = family(prefix)
        if r is None:
            rows.append((name, None)); allok = False; continue
        rows.append((name, r)); allok = allok and r['complete']
        print(name, r, flush=True)
    L = ['# Results of the certificate runs', '',
         'Generated by `code/check_all.py`. A family is complete when every split tree replays exactly and every',
         'incomplete stack or failed leaf is re-covered by a complete continuation run (see check_all.py).', '',
         '| run | files | records | certified leaves | pruned | failed leaves (rechecked) | open boxes | replay problems | bound histogram | status |',
         '|---|---|---|---|---|---|---|---|---|---|']
    for name, r in rows:
        if r is None:
            L.append(f'| {name} | – | – | – | – | – | – | – | – | **missing** |')
            continue
        L.append(f"| {name} | {r['files']} | {r['records']} | {r['leaves']} | {r['pruned']} | {r['fails']} | {r['owed']} | "
                 f"{r['nproblems']} | {r['hist']} | {'complete' if r['complete'] else '**INCOMPLETE**'} |")
    L += ['', '**All families complete: ' + ('YES' if allok else 'NO') + '**', '']
    open('../RESULTS.md', 'w').write('\n'.join(L) + '\n')
    print('ALL COMPLETE' if allok else 'NOT COMPLETE')

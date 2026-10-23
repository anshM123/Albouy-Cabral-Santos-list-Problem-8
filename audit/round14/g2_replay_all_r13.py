"""AUDIT-7 g2 (final pass): independent replay of EVERY recorded run family, without check_all / bb_replay code.

For each family (bb_run, kq_run and c3_run outputs; NAME.json or NAME.json.gz, plus continuations NAME.r1.json, ...):
  * header: engine/chart, domain, n0 and exclusions equal the values of code/launch_runs.sh (hard-coded here from my
    own reading of launch_runs.sh, 2026-10-23), split scale consistent;
  * initial boxes: exactly the np.linspace grid of the domain minus the skipped boxes, no duplicates, every skipped box
    inside a declared exclusion block;
  * every tree: own replay with its file's recorded scale (continuations may use another scale) (split at the
    float midpoint of the largest width/scale, ties -> lowest index, children
    pushed lower then upper), code consumed exactly, done -> empty stack, not done -> stack == 'left', '!' leaves
    == 'fails' in order, leaf count == nleaf, number of winding leaves == len(wL) (== nwind when stored), volume sum;
  * leaf codes allowed per family, and the per-family bound limits of the proof:
      whole-circle engines (a, f, C1*, B*, L2*, D2*, KQ): 'a'..'l' / 'A'..'L' (11 = parity), '1'..'6' (W, bound 12-2W),
          'u' (D2 only: unphysical);
      KA (away arc of the C3 cell): bound <= 2;
      C3 charts: near/partial count <= target (Ca 9; CK, CE, Cm 8; L2*E, L2*m 4; Phi*, L2*Y* 0), zone leaves
          ('A'+c, Ca and CK only) c <= 5, 'u'/'o' pruned;
    winding radii: -15 <= L < 0;
  * owed boxes ('!' leaves and open stacks) must be the roots of complete continuation records, except
      C1p: owed boxes must lie in the D2 cube {e <= 5e-3, p <= 5e-3, 1 - q <= 5e-3, |tau - 1/2| <= 5e-3}
           (checked with outward-rounded bounds) and the D2 families must be complete;
      D2*: an owed D2 box is excused only if its smallest e exceeds the largest e of every owed C1p box;
  * pruned boxes ('u', 'o'): sampled points (corners + random) must be unphysical / outside the C3 cell.
  * round 11 (corner K): families KQF2 (outer count from 3 SF2) and KQS (whole circle, SF2 < s <= SF); their header
    fields (FACE_S, T_X, c3prune, engine) are checked in EVERY file and the base domains against launch_runs.sh; the 'p'
    radius is the family's role (SF2 for KQS, SF = 0.02 otherwise), not the header value; KQ 'o' only in KQS (exact
    C3-cell predicate); face-chart boxes inside K with s_lo > SF2 and KQF boxes inside K_lw are excused; a family with
    'p' leaves needs FE, FM, FIP, FIM, KQF, KQF2 and KQS complete.
usage: python g2_replay_all.py [--prune-samples K] [--only NAME,...]"""
import os, sys, json, gzip, math, random, itertools, collections
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', 'code'))
sys.path.insert(0, ROOT)
import splitrule                                  # round 9: header split_rule 'face' (kq_run --facesplit)
args = sys.argv[1:]
KPR = int(args[args.index('--prune-samples') + 1]) if '--prune-samples' in args else 6
ONLY = set(args[args.index('--only') + 1].split(',')) if '--only' in args else None
MAXLEV = {kv.split(':')[0]: int(kv.split(':')[1]) for kv in os.environ.get('AUDIT_MAXLEV', '').split(';') if kv}
random.seed(77)
HP = 1.5707963267948966
W3, W9 = 0.47123889803846897, 1.413716694115407
NU1, R1, R2, ETA1 = 0.37, 1.001, 1.001, 1.0001
SD2 = 5e-3
WHOLE = ('a', 'f', 'C1', 'C1p', 'C1pp', 'B1', 'B2', 'B3p', 'B3m', 'L2E', 'L2P', 'L2Ap', 'L2Am', 'L2Tp', 'L2Tm',
         'D2E', 'D2E2', 'D2P', 'D2K', 'D2Tp', 'D2Tm', 'KQ')
C3DOM = {'Ca': (0, NU1, 0, R1**2, 0, R1**2, 0, R1**3), 'CK': (0, NU1, 0, R1, 0, R1**2, 0, R1**3),
         'CE': (0, NU1, 0, R1, 0, R1**2, 0, R1**3), 'Cm': (0, NU1, 0, R1, 0, R1**2, 0, R1**2),
         'PhiP': (0, NU1, 0, R1, 0, R1**2, 0, R1**2), 'PhiM': (0, NU1, 0, R1, 0, R1**2, 0, R1**2)}
for c in ('L2aE', 'L2am', 'L2aYp', 'L2aYm'):
    C3DOM[c] = (0, NU1, 0, R1**2, 0, ETA1, 0, R2**2)
for c in ('L2KE', 'L2Km', 'L2KYp', 'L2KYm'):
    C3DOM[c] = (0, NU1, 0, R1, 0, ETA1, 0, R2**2)
C3TARGET = {'FE': 5, 'FM': 5, 'FIP': 0, 'FIM': 0, 'Ca': 9, 'CK': 8, 'CE': 8, 'Cm': 8, 'PhiP': 0, 'PhiM': 0, 'L2aE': 4, 'L2am': 4, 'L2KE': 4, 'L2Km': 4,
            'L2aYp': 0, 'L2aYm': 0, 'L2KYp': 0, 'L2KYm': 0}
ZONE_OUTER = 5
L2D = (0, 1.01, 0, 1.05, -1.05, 1.05, 0, 1)
FAM = [  # name, path prefix (relative to code/), kind, engine/chart, domain, n0, exclusions
    ('bulk', 'runs/bulk', 'bb', 'a', (0.20033740617957885, 1.2535563099909826, 0, 1, 0, HP, 0, 1), (8, 10, 10, 10),
     [(0.9902515840381316, 1.2535563099909826, 0.0, 1.0, 0.0, HP, 0.0, 1.0)]),          # round 9: eta >= 0.99025 excluded
    ('F1', 'runs/F1', 'bb', 'f', (0.835, 0.9501, 0, 1, 0, HP, 0, 0.5), (2, 10, 10, 2), []),        # round 9
    ('KQ0', 'runs/KQ0', 'kq', 'KQ', None, None, []),                                               # round 9 (header)
    ('KQ0T', 'runs/KQ0T', 'kq', 'KQ', None, None, []),                                             # round 9 (header)
    ('KQF', 'runs/KQF', 'kqf', 'KQF', None, None, []),                                             # round 9 (header)
    ('KQF2', 'runs/KQF2', 'kqf', 'KQF', None, None, []),                                           # round 11 (header, R11)
    ('KQS', 'runs/KQS', 'kq', 'KQ', None, None, []),                                               # round 11 (header, R11)
    ('C1', 'runs/C1', 'bb', 'C1', (0, 0.26, 0, 1.1071487177940904, 0, HP, 0, 1), (8, 12, 4, 16),
     [(0, 0.0650001, 0, 0.1846, 0, 1.5708, 0.375, 0.625)]),
    ('C1p', 'runs/C1p', 'bb', 'C1p', (0, 0.2001, 0.4, 1, 0, HP, 0, 1), (4, 4, 4, 4), []),
    ('C1pp', 'runs/C1pp', 'bb', 'C1pp', (0, 0.51, 0, 0.52, 0, HP, 0, 1), (4, 4, 4, 4), []),
    ('B1', 'runs/B1', 'bb', 'B1', (0, 0.07, 0, 1.75, 0, 1.75, 0, 1), (7, 7, 7, 8),
     [(0, 0.0100001, 0, 0.2500001, 0.7499999, 1.2500001, 0.3749999, 0.6250001)]),
    ('B2', 'runs/B2', 'bb', 'B2', (0, 0.2, 0, 0.3675, 0, HP, 0, 1), (4, 4, 4, 4), []),
    ('B3p', 'runs/B3p', 'bb', 'B3p', (0, 0.13, 0, 0.5654, 0, 1.6154, 0, HP), (4, 4, 4, 4), []),
    ('B3m', 'runs/B3m', 'bb', 'B3m', (0, 0.13, 0, 0.5654, 0, 1.6154, 0, HP), (4, 4, 4, 4), []),
    ('L2E', 'runs/L2E', 'bb', 'L2E', L2D, (4, 4, 4, 4), []),
    ('L2P', 'runs/L2P', 'bb', 'L2P', L2D, (4, 4, 4, 4), []),
    ('L2Ap', 'runs/L2Ap', 'bb', 'L2Ap', (0, 1.01, 0, 1.05, 0, 1.05, 0, 1), (4, 4, 4, 4), []),
    ('L2Am', 'runs/L2Am', 'bb', 'L2Am', (0, 1.01, 0, 1.05, 0, 1.05, 0, 1), (4, 4, 4, 4), []),
    ('L2Tp', 'runs/L2Tp', 'bb', 'L2Tp', (0, 1.01, 0, 1.05, 0, 1.05, -1.05, 1.05), (4, 4, 4, 4), []),
    ('L2Tm', 'runs/L2Tm', 'bb', 'L2Tm', (0, 1.01, 0, 1.05, 0, 1.05, -1.05, 1.05), (4, 4, 4, 4), []),
    ('D2E', 'runs/D2E', 'bb', 'D2E', (0, 1.01, 0, 1.05, 0, 1.05, 0, 1), (4, 4, 4, 4), []),
    ('D2E2', 'runs/D2E2', 'bb', 'D2E2', (0, 0.14, 0, 0.35, 0.24, 0.68, 0, 1), (4, 4, 4, 4), []),   # round-7 addendum
    ('D2P', 'runs/D2P', 'bb', 'D2P', (0, 1.01, 0, 1.05, 0, 1.05, 0, 1), (4, 4, 4, 4), []),
    ('D2K', 'runs/D2K', 'bb', 'D2K', (0, 1.01, 0, 1.05, 0, 1.05, 0, 1), (4, 4, 4, 4), []),
    ('D2Tp', 'runs/D2Tp', 'bb', 'D2Tp', (0, 1.01, 0, 1.05, 0, 1.05, 0, 1.05), (4, 4, 4, 4), []),
    ('D2Tm', 'runs/D2Tm', 'bb', 'D2Tm', (0, 1.01, 0, 1.05, 0, 1.05, 0, 1.05), (4, 4, 4, 4), []),
    ('F0', 'runs/F0', 'bb', 'f', (0.9499, 1.0, 0, 1, 0, HP, 0, 0.5), (4, 5, 5, 4), []),
    ('KQ1', 'runs/KQ1', 'kq', 'KQ', (0.9499, 1.0, 0, 1, 0, HP, 0.5, 0.9975), (4, 10, 10, 5), []),
    ('KQ2', 'runs/KQ2', 'kq', 'KQ', (0.9499, 0.99, 0, 1, 0, HP, 0.9975, 1.0), (4, 10, 10, 1), []),
    ('KQK1', 'runs/KQK1', 'kq', 'KQ', (0.99, 1.0, 0, 0.95, 0, HP, 0.997456, 1.0), (2, 19, 10, 2), []),
    ('KQK2', 'runs/KQK2', 'kq', 'KQ', (0.99, 1.0, 0.95, 1.0, 0, W3, 0.997456, 1.0), (2, 1, 3, 2), []),
    ('KQB', 'runs/KQB', 'kq', 'KQ', (0.99, 1.0, 0.95, 1.0, W3, W9, 0.997456, 1.0), (2, 5, 6, 2), []),
    ('KA', 'runs/KA', 'bb', 'KA', (0.99, 1.0, 0.95, 1.0, W9, HP, 0, 0.0505), (2, 2, 4, 4), []),
] + [('face ' + c, 'c3/runs/' + c, 'c3', c, (0, 0.02 * R1, 0, R1, 0, 1, 0, 1), None, []) for c in ('FE', 'FM', 'FIP', 'FIM')] + \
    [('C3 ' + c, 'c3/runs/' + c, 'c3', c, C3DOM[c], (8, 4, 4, 4), []) for c in
     ('Ca', 'CK', 'CE', 'Cm', 'PhiP', 'PhiM', 'L2aE', 'L2am', 'L2aYp', 'L2aYm', 'L2KE', 'L2Km', 'L2KYp', 'L2KYm')]


def load(path):
    for p in (path, path + '.gz'):
        if os.path.exists(p):
            return json.load(gzip.open(p, 'rt') if p.endswith('.gz') else open(p))
    return None


def vol(b):
    v = 1.0
    for i in range(4):
        v *= b[2 * i + 1] - b[2 * i]
    return v


def grid(dom, n0):
    edges = [np.linspace(dom[2 * i], dom[2 * i + 1], n0[i] + 1) for i in range(4)]
    return [tuple(float(v) for v in (edges[0][a], edges[0][a + 1], edges[1][b], edges[1][b + 1], edges[2][c],
                                     edges[2][c + 1], edges[3][d], edges[3][d + 1]))
            for a in range(n0[0]) for b in range(n0[1]) for c in range(n0[2]) for d in range(n0[3])]


def inside(bx, ex):
    return all(ex[2 * i] <= bx[2 * i] and bx[2 * i + 1] <= ex[2 * i + 1] for i in range(4))


# ---- pruning predicates at a point (my own maps; C3 maps as in audit round3/c2, D2 maps as in round6/f3) ----
def d2_lq(eng, x):
    n = x[0]; SD = mp.mpf(5) / 1000
    if eng in ('D2Tp', 'D2Tm'):
        E, P, K = x[1], x[2], x[3]; p, k = SD * P * n, SD * K * n
    elif eng in ('D2E', 'D2E2'):
        p, k = SD * x[1] * n, SD * x[2] * n
    elif eng == 'D2P':
        p, k = SD * n, SD * x[2] * n
    else:
        p, k = SD * x[2] * n, SD * n
    return p * p + (1 - k) ** 2            # l^2


def c3_point(ch, x):
    """(Qz, kappa, eps, mu) of a C3 chart box point (root variable excluded)."""
    nu = x[0]
    if ch == 'Ca':
        return nu, nu * nu * x[1], nu * nu * x[2], nu ** 3 * x[3]
    if ch == 'CK':
        return nu * x[1], nu * nu, nu * nu * x[2], nu ** 3 * x[3]
    if ch == 'CE':
        return nu * x[1], nu * nu * x[2], nu * nu, nu ** 3 * x[3]
    if ch == 'Cm':
        return nu * x[1], nu * nu * x[2], nu * nu * x[3], nu ** 3
    if ch in ('PhiP', 'PhiM'):
        return nu * x[1], nu * nu * x[2], nu * nu * x[3], None          # mu = nu^3 * root (root >= 0): min 0
    rv, eta, w4 = x[1], x[2], x[3]
    qz, kap = (nu, nu * nu * rv) if ch.startswith('L2a') else (nu * rv, nu * nu)
    kind = ch[3:]
    if kind == 'E':
        return qz, kap, nu * nu * mp.mpf('0.01') * eta * eta, nu ** 3 * mp.mpf('0.01') * eta * eta * w4
    if kind == 'm':
        return qz, kap, nu * nu * mp.mpf('0.01') * eta * eta * w4, nu ** 3 * mp.mpf('0.01') * eta * eta
    return qz, kap, nu * nu * mp.mpf('0.01') * eta * eta * w4, None              # *Yp/*Ym: mu = ... * root, min 0


mp.mp.dps = 40
QZMAX = mp.cos(mp.mpf(W9)); KAPMAX = 1 - mp.mpf(0.95) * mp.sin(mp.mpf(W9)); EPSMAX = mp.mpf('0.01'); MUMAX = mp.mpf('0.0505')


FACE_S2 = mp.mpf('0.02') ** 2 if False else None
# round 11 (corner K, C3_CORE.md section 6): the face construction is used only for s <= SF2 inside K; the whole-circle
# run KQS prunes 'p' only into {eps, mu <= SF2^2} (covered by FE/FM/FIP/FIM with s <= SF2 and KQF2 from 3 SF2) and 'o'
# into the C3 cell.  The 'p' radius is fixed by the ROLE of the family (not taken from the run-written header):
SF, SF2 = 0.02, 0.010005
P_RADIUS = {'KQS': SF2}                        # every other KQ family: SF
R11 = {'KQF2': dict(domain=[0.99989989, 1.0, 0.88, 1.0, 1.29, HP, 0.9999999899, 1.0], every=dict(engine='KQF', FACE_S=SF2, T_X=0.030015, LIMIT_OUTER=6)),
       'KQS': dict(domain=[0.9996, 1.0, 0.88, 1.0, 1.29, HP, 0.99999983, 1.0], every=dict(engine='KQ', FACE_S=SF2, c3prune=True)),
       'KQF': dict(domain=None, every=dict(engine='KQF', FACE_S=SF, T_X=0.06, LIMIT_OUTER=6))}
K_PZ, K_QZ0, K_QZ1, K_L, K_W = 0.875, 0.125, 0.25, 0.88, 1.29


def face_prune_ok(code, eng, bx, fs=SF, c3p=False):
    """exact checks (mpmath): KQ 'p' = whole box in {eps <= fs^2, mu <= fs^2} (fs = the family's role radius); KQF 'o'
    (and KQ 'o' of a --c3prune run) = whole box in the C3 cell directions l >= 0.95, w >= 9pi/20 (float literal of the
    C3 core); face-chart 'u' = Pz^2 + Qz^2 > 1 on the box, face-chart 'o' = l >= 0.95 and w >= 9pi/20 on the box
    (Pz, Qz >= 0)."""
    S2 = mp.mpf(fs) ** 2
    if eng == 'KQ' and code == 'p':
        eps_hi = 1 - mp.mpf(bx[0]); t0 = mp.mpf(bx[6]); mu_hi = mp.sqrt((1 - t0) / t0)
        return eps_hi <= S2 and mu_hi <= S2
    if (eng == 'KQF' or (eng == 'KQ' and c3p)) and code == 'o':
        # round 11: also inside the region where the C3 core covers the cell (e >= 0.99, mu <= 0.0505; C3_CORE.md)
        t0 = mp.mpf(bx[6]); in_core = 1 - mp.mpf(bx[0]) <= EPSMAX and mp.sqrt((1 - t0) / t0) <= MUMAX
        return mp.mpf(bx[2]) >= mp.mpf(0.95) and mp.mpf(bx[4]) >= mp.mpf(W9) and in_core
    if eng in ('FE', 'FM', 'FIP', 'FIM'):
        pz0, qz0, qz1 = mp.mpf(bx[4]), mp.mpf(bx[6]), mp.mpf(bx[7])
        if code == 'u':
            return pz0 * pz0 + qz0 * qz0 > 1
        if code == 'o':                       # C3 cell, or (round 10b) s_lo > SF: no needed point (AUDIT_9 excuse rule)
            in_c3 = mp.sqrt(pz0 * pz0 + qz0 * qz0) >= mp.mpf(0.95) and (qz1 <= 0 or mp.atan2(pz0, qz1) >= mp.mpf(W9))
            return in_c3 or bx[0] > 0.02
    return False


def pruned_ok(kind, code, eng, bx):
    pts = list(itertools.product(*[bx[2 * i:2 * i + 2] for i in range(4)]))
    pts = random.sample(pts, min(KPR, len(pts))) + [tuple(random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(4))
                                                   for _ in range(KPR)]
    for x in pts:
        if x[0] <= 0:                         # nu = 0: the blow-up face (e = 0 for D2, e = 1 for C3), not physical
            continue
        x = [mp.mpf(v) for v in x]            # exact predicate (float rounding of l^2 - 1 gave a false alarm)
        if kind == 'bb':                      # D2 'u': l > 1
            if not d2_lq(eng, x) > 1:
                return False
        else:
            qz, kap, eps, mu = c3_point(eng, x)
            l2 = (1 - kap) ** 2 + qz * qz
            if code == 'u' and not l2 > 1:
                return False
            if code == 'o' and not (qz > QZMAX or kap > KAPMAX or eps > EPSMAX or (mu is not None and mu > MUMAX)):
                return False
    return True


TW_ENG = ('KQ', 'KQF', 'Ca', 'CK', 'CE', 'Cm')     # round 13: '7' = tau-winding certificate, whole circle <= 10


def allowed(kind, eng, ch):
    """bound of a leaf code, or None if the code is not allowed for this family"""
    if ch == '7':
        return 10 if eng in TW_ENG else None
    if kind == 'kqf':
        return ord(ch) - ord('a') if 'a' <= ch <= 'g' else None      # outer count <= LIMIT_OUTER = 6
    if kind in ('bb', 'kq'):
        if eng == 'KA':
            return ord(ch) - ord('a') if 'a' <= ch <= 'c' else None
        if 'a' <= ch <= 'l':
            return ord(ch) - ord('a')
        if 'A' <= ch <= 'L' and kind == 'bb':
            return ord(ch) - ord('A')
        if '1' <= ch <= '6' and eng in ('a', 'f', 'C1', 'C1p', 'C1pp', 'B1', 'B2', 'L2E', 'L2P', 'L2Ap', 'L2Am', 'D2E',
                                         'D2E2', 'D2P', 'D2K', 'KQ'):
            return 12 - 2 * int(ch)
        return None
    t = C3TARGET[eng]
    if 'a' <= ch <= 'z' and ord(ch) - ord('a') <= t:
        return ord(ch) - ord('a')
    if 'A' <= ch <= 'Z' and eng in ('Ca', 'CK') and ord(ch) - ord('A') <= ZONE_OUTER:
        return 100 + ord(ch) - ord('A')        # zone leaf: outer count c (+ inner <= 4)
    return None


def replay_family(name, prefix, kind, eng, dom, n0, excl):
    base = os.path.join(ROOT, prefix)
    D = load(base + '.json')
    if D is None:
        return dict(status='MISSING'), [], []
    files = [D]; k = 1
    maxk = MAXLEV.get(name)                       # AUDIT-12 snapshot: continuation levels up to K only (AUDIT_MAXLEV)
    while maxk is None or k <= maxk:
        Dk = load(base + '.r%d.json' % k)
        if Dk is None:
            break
        files.append(Dk); k += 1
    H = D['header']
    probs = []
    if name in R11:                               # round 11: header fields of every file, base domain
        exp = R11[name]
        if exp['domain'] is not None and [float(x) for x in H.get('domain') or []] != [float(x) for x in exp['domain']]:
            probs.append('round-11 base domain %s != %s' % (H.get('domain'), exp['domain']))
        for Df in files:
            for kk, vv in exp['every'].items():
                if Df['header'].get(kk) != vv:
                    probs.append('round-11 header field %s = %r (expected %r)' % (kk, Df['header'].get(kk), vv))
    p_fs = P_RADIUS.get(name, SF)                  # 'p' radius of this family (role)
    c3p_ok = name == 'KQS'                         # KQ 'o' leaves allowed only in the --c3prune run KQS
    heng = H.get('engine') if kind != 'c3' else H.get('chart')
    if heng != eng:
        probs.append('engine/chart %s != %s' % (heng, eng))
    hdr_only = dom is None                       # round 9 families: domain and n0 taken from the base header (reported)
    if n0 is None:
        n0 = H['n0']
    if dom is None:
        dom = H['domain']
    elif [float(v) for v in H['domain']] != [float(v) for v in dom] or list(H['n0']) != list(n0):
        probs.append('domain/n0 differ from launch_runs.sh: %s %s' % (H['domain'], H['n0']))
    if kind == 'bb' and [list(map(float, e)) for e in H.get('exclude', [])] != [list(map(float, e)) for e in excl]:
        probs.append('exclusions differ from launch_runs.sh')
    S = H['scale']
    if hdr_only or kind == 'kqf':
        exp_scale = S                             # header scale (reported)
    elif kind == 'kq' or eng == 'f':
        exp_scale = [1.0, 1.0, math.pi / 2, 1.0]
    elif eng == 'a':
        exp_scale = [math.pi / 2, 1.0, math.pi / 2, 1.0]
    else:
        exp_scale = [dom[2 * i + 1] - dom[2 * i] for i in range(4)]
    if any(abs(a - b) > 1e-15 * max(1, abs(b)) for a, b in zip(S, exp_scale)):
        probs.append('scale %s unexpected (%s)' % (S, exp_scale))
    G = grid(dom, n0)
    roots0 = [tuple(r['box']) for r in D['results']]
    skipped = [tuple(b) for b in H.get('skipped', [])]
    if len(set(roots0)) != len(roots0):
        probs.append('duplicate initial boxes: %d' % (len(roots0) - len(set(roots0))))
    if set(roots0) | set(skipped) != set(G) or set(roots0) & set(skipped):
        probs.append('initial boxes != grid minus skipped (%d grid, %d roots, %d skipped, %d missing)' %
                     (len(G), len(set(roots0)), len(skipped), len(set(G) - set(roots0) - set(skipped))))
    for b in skipped:
        if not any(inside(b, e) for e in excl):
            probs.append('skipped box outside every exclusion'); break
    uses_p = [False]; ntw_ = [0]                     # round 13: number of tau-winding leaves
    hist = {}; nleaf = 0; nw = 0; Lmin, Lmax = 0.0, -99.0; owed = set(); nprune = 0; prune_bad = 0; nfail = 0; ninc = 0
    leaves_out = []
    for fi, Df in enumerate(files):
        Sf = Df['header'].get('scale') if fi > 0 else S        # a continuation replays with its own recorded scale
        rule_f = Df['header'].get('split_rule')
        if fi > 0 and (Sf is None or len(Sf) != 4 or not all(v > 0 for v in Sf)):
            probs.append('continuation without a valid scale'); continue
        if fi > 0 and Df['header'].get('engine', Df['header'].get('chart')) != heng:
            probs.append('continuation engine/chart differs')
        for r in Df['results']:
            root = tuple(r['box'])
            if fi > 0:
                if root in owed:
                    owed.discard(root)
                else:
                    probs.append('continuation root not owed')
            stack = [root]; fl = [tuple(f[0]) for f in r.get('fails', [])]; fpos = 0; wl = r.get('wL', []); wpos = 0
            twl_ = r.get('tw', []); twpos = 0                # round 13: tau-winding contours, in tree order
            v = 0.0; cnt = 0; ndig = 0
            rc = collections.Counter()                     # round 11: replayed code histogram vs the record's 'hist'
            for ch in r['tree']:
                if not stack:
                    probs.append('stack underflow'); break
                bx = stack.pop()
                if ch == '+':
                    if rule_f is not None:                         # header split_rule (round 9): splitrule.split_box
                        a, c = splitrule.split_box(bx, Sf, rule_f)
                        diff = [i for i in range(8) if a[i] != bx[i] or c[i] != bx[i]]
                        if not (len(diff) == 2 and diff[1] == diff[0] + 1 and diff[0] % 2 == 0 and
                                a[diff[1]] == c[diff[0]] and bx[diff[0]] < a[diff[1]] < bx[diff[1]] and
                                a[diff[0]] == bx[diff[0]] and c[diff[1]] == bx[diff[1]]):
                            probs.append('split does not tile the parent box')
                        stack.append(tuple(a)); stack.append(tuple(c))
                        continue
                    w = [(bx[2 * i + 1] - bx[2 * i]) / Sf[i] for i in range(4)]
                    i = max(range(4), key=lambda j: w[j])          # first maximum, as np.argmax
                    mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
                    a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
                    stack.append(tuple(a)); stack.append(tuple(c))
                    continue
                v += vol(bx)
                if ch == '!':
                    if fpos >= len(fl) or fl[fpos] != bx:
                        probs.append('fail mismatch')
                    fpos += 1
                    continue
                cnt += 1
                if (ch == 'p' and eng == 'KQ') or (ch == 'o' and (kind == 'kqf' or (eng == 'KQ' and c3p_ok))) or \
                        (ch in 'uo' and eng in ('FE', 'FM', 'FIP', 'FIM')):
                    nprune += 1; uses_p[0] = uses_p[0] or ch == 'p'; rc[ch] += 1
                    if not face_prune_ok(ch, eng, bx, p_fs, c3p_ok):
                        prune_bad += 1
                    continue
                if ch in 'uo' and (kind == 'c3' or (ch == 'u' and eng.startswith('D2'))):
                    nprune += 1
                    if KPR > 0 and (nprune <= 4000 or random.random() < 0.02):
                        if not pruned_ok(kind, ch, eng, bx):
                            prune_bad += 1
                    continue
                b = allowed(kind, eng, ch)
                if b is None:
                    probs.append('code %r not allowed for %s' % (ch, name)); continue
                if ch == '7':                                  # round 13: tau-winding leaf, contour [which, cre, cim, r]
                    if twpos >= len(twl_):
                        probs.append('missing tw contour')
                    else:
                        t_ = twl_[twpos]
                        if not (len(t_) == 4 and t_[0] in ((0, 1) if kind in ('kq', 'kqf') else (0,)) and t_[2] > t_[3] > 0):
                            probs.append('invalid tw contour %r' % (t_,))
                    twpos += 1; ntw_[0] += 1
                    hist['tw'] = hist.get('tw', 0) + 1; rc['t'] += 1; nleaf += 1
                    continue
                if '1' <= ch <= '6':
                    if wpos >= len(wl):
                        probs.append('missing wL'); L = None
                    else:
                        L = wl[wpos]
                        Lmin = min(Lmin, L); Lmax = max(Lmax, L)
                        if not (-15.0 <= L < 0.0):
                            probs.append('winding radius out of range %r' % L)
                    wpos += 1; ndig += 1; nw += 1
                    key = 'w%s' % ch
                else:
                    L = None
                    key = b
                hist[key] = hist.get(key, 0) + 1
                rc[key if isinstance(key, str) else str(key)] += 1
                nleaf += 1
            if r['done'] and stack:
                probs.append('stack left in a done record')
            left = [tuple(x) for x in r.get('left', [])]
            if not r['done']:
                ninc += 1
                if stack != left:
                    probs.append('left differs from replay')
            v += sum(vol(x) for x in left)
            if fpos != len(fl):
                probs.append('unconsumed fails')
            if twpos != len(twl_):
                probs.append('number of 7 leaves differs from len(tw)')
            if wpos != len(wl) or ('nwind' in r and r['nwind'] != ndig):
                probs.append('winding count differs from wL/nwind')
            if cnt != r['nleaf']:
                probs.append('nleaf differs (%d vs %d)' % (cnt, r['nleaf']))
            # round 11 (R11-3, concurrent writers): a line spliced from two records would not reproduce its own 'hist'
            if kind in ('kq', 'kqf') and 'hist' in r and dict(rc) != {str(k): v_ for k, v_ in r['hist'].items()}:
                probs.append('record hist differs from the replayed leaf codes')
            if abs(v - vol(root)) > 1e-9 * max(vol(root), 1e-300):
                probs.append('volume mismatch')
            for x in left + fl:
                owed.add(x)
            nfail += len(fl)
    st = dict(files=len(files), records=sum(len(Df['results']) for Df in files), leaves=nleaf, winding=nw,
              L_range=(Lmin, Lmax) if nw else None, pruned=nprune, prune_violations=prune_bad, fails=nfail,
              incomplete=ninc, owed=len(owed), hist=dict(sorted(hist.items(), key=lambda kv: str(kv[0]))),
              domain=H['domain'], n0=H['n0'], scale=S, split_rules=[Df['header'].get('split_rule') for Df in files],
              uses_p=uses_p[0], tw_leaves=ntw_[0],
              problems=probs)
    return st, sorted(owed), leaves_out


def c1p_owed_in_cube(b):
    """C1p box (e, l, w, tau): e <= SD2, p = l sin w <= SD2, kappa = 1 - l cos w <= SD2, |tau - 1/2| <= SD2,
    with outward-rounded upper bounds (mpmath, 50 digits)"""
    mp.mp.dps = 50
    e_hi = mp.mpf(b[1]); p_hi = mp.mpf(b[3]) * mp.sin(mp.mpf(b[5])); k_hi = 1 - mp.mpf(b[2]) * mp.cos(mp.mpf(b[5]))
    return e_hi <= SD2 and p_hi <= SD2 and k_hi <= SD2 and mp.mpf(b[6]) >= mp.mpf(0.5) - mp.mpf(SD2) and \
        mp.mpf(b[7]) <= mp.mpf(0.5) + mp.mpf(SD2)


if __name__ == '__main__':
    res = {}; owed_all = {}
    for fam in FAM:
        if ONLY and fam[0] not in ONLY:
            continue
        st, owed, _ = replay_family(*fam)
        res[fam[0]] = st; owed_all[fam[0]] = owed
        print(fam[0], json.dumps(st, default=str), flush=True)
    # owed-box rules
    verdict = {}
    c1p_owed = owed_all.get('C1p', [])
    if 'C1p' in res and res['C1p'].get('status') != 'MISSING':
        bad = [b for b in c1p_owed if not c1p_owed_in_cube(b)]
        emax = max((b[1] for b in c1p_owed), default=-1.0)
        print(f'C1p owed boxes {len(c1p_owed)}; outside the D2 cube {len(bad)}; max e of owed boxes {emax:.3e}')
        mp.mp.dps = 50
        emax_x = max((mp.mpf(b[1]) for b in c1p_owed), default=mp.mpf(-1))
        for nm in ('D2E', 'D2E2', 'D2P', 'D2K', 'D2Tp', 'D2Tm'):
            if nm not in owed_all:
                continue
            ob = owed_all[nm]
            # rule (i), exact: smallest physical e of the box (SD2 nu, or SD2 E nu) > largest e of the owed C1p boxes
            need = [b for b in ob if not ((mp.mpf(SD2) * mp.mpf(b[0]) * (1 if nm in ('D2E', 'D2E2') else mp.mpf(b[2]))) > emax_x)]
            print(f'{nm}: owed {len(ob)}, owed boxes that may meet an owed C1p box {len(need)}')
            owed_all[nm] = need
        # rule (ii), exact: a D2Tp box (nu, E, P, K), E_lo > 0, whose D2E image (E nu, P/E, K/E, T = 1/E) lies in D2E's
        # domain for T <= 1.05 and in D2E2's domain for T >= 1.04, both families complete (checked below)
        def covered_ii(b):
            n0, n1, E0, E1, P0, P1, K0, K1 = [mp.mpf(v) for v in b]
            if not E0 > 0:
                return False
            nu1, PE1, KE0, KE1, T0, T1 = E1 * n1, P1 / E0, K0 / E1, K1 / E0, 1 / E1, 1 / E0
            inE = lambda: nu1 <= mp.mpf('1.01') and PE1 <= mp.mpf('1.05') and KE1 <= mp.mpf('1.05') and T0 >= mp.mpf('-1.05')
            inE2 = lambda: nu1 <= mp.mpf('0.14') and PE1 <= mp.mpf('0.35') and KE0 >= mp.mpf('0.24') and KE1 <= mp.mpf('0.68')                 and T1 <= mp.mpf('1.28')
            if T1 <= mp.mpf('1.05'):
                return inE()
            if T0 >= mp.mpf('1.04'):
                return inE2()
            return inE() and inE2()
        if 'D2Tp' in owed_all and owed_all['D2Tp']:
            def fam_ok(nm):
                st = res.get(nm)
                return bool(st) and st.get('status') != 'MISSING' and not st['problems'] and not owed_all.get(nm)                     and st['prune_violations'] == 0
            okE, okE2 = fam_ok('D2E'), fam_ok('D2E2')
            left_ii = [b for b in owed_all['D2Tp'] if not (okE and okE2 and covered_ii(b))]
            print(f'D2Tp rule (ii): {len(owed_all["D2Tp"]) - len(left_ii)} owed boxes covered by D2E u D2E2 '
                  f'(D2E complete: {okE}, D2E2 complete: {okE2}); still open {len(left_ii)}')
            owed_all['D2Tp'] = left_ii
        if not bad:
            owed_all['C1p'] = []
    # face charts: open or failed boxes with s_lo > SF = 0.02 contain no needed point (round 9 reduction, AUDIT_9)
    for nm in ('face FE', 'face FM', 'face FIP', 'face FIM'):
        if nm in owed_all and owed_all[nm]:
            keep = [b for b in owed_all[nm] if not b[0] > 0.02]
            print(f'{nm}: owed {len(owed_all[nm])}, of which with s_lo > 0.02 (not needed) {len(owed_all[nm]) - len(keep)}')
            # round 11: inside K = {Pz >= 0.875, 0.125 <= Qz <= 0.25} only s <= SF2 is needed (KQF2 from 3 SF2, KQS above
            # SF2; K lies in K_lw = {l >= 0.88, w >= 1.29}, checked in audit/round11/i1_seams.py)
            keep2 = [b for b in keep if not (b[4] >= K_PZ and b[6] >= K_QZ0 and b[7] <= K_QZ1 and b[0] > SF2)]
            print(f'{nm}: of the rest, inside K with s_lo > SF2 (not needed) {len(keep) - len(keep2)}; still owed {len(keep2)}')
            owed_all[nm] = keep2
    # round 11: KQF boxes inside K_lw are not needed (KQF2 for s <= SF2, KQS for SF2 < s <= SF)
    if owed_all.get('KQF'):
        keep = [b for b in owed_all['KQF'] if not (b[2] >= K_L and b[4] >= K_W)]
        print(f"KQF: owed {len(owed_all['KQF'])}, inside K_lw (not needed) {len(owed_all['KQF']) - len(keep)}; still owed {len(keep)}")
        owed_all['KQF'] = keep
    def fam_complete(nm):
        st = res.get(nm)
        return bool(st) and st.get('status') != 'MISSING' and not st['problems'] and not owed_all.get(nm)             and st['prune_violations'] == 0
    face_names = ('face FE', 'face FM', 'face FIP', 'face FIM', 'KQF', 'KQF2', 'KQS')     # round 11: + KQF2, KQS
    for nm, st in res.items():
        if st.get('status') != 'MISSING' and st.get('uses_p'):
            deps = [f for f in face_names if f != nm]
            missing = [f for f in deps if f not in res]
            if missing:
                st['problems'].append('uses p: face families not checked in this run: %s' % missing)
            elif not all(fam_complete(f) for f in deps):
                st['problems'].append('uses p but the face families (FE, FM, FIP, FIM, KQF, KQF2, KQS) are not all complete')
    allok = True
    for nm, st in res.items():
        ok = st.get('status') != 'MISSING' and not st['problems'] and not owed_all[nm] and st['prune_violations'] == 0
        allok = allok and ok
        print(f"{nm:10s} {'OK' if ok else 'NOT OK'}  leaves {st.get('leaves')}  winding {st.get('winding')}  "
              f"fails {st.get('fails')}  incomplete {st.get('incomplete')}  open {len(owed_all.get(nm, []))}  "
              f"problems {len(st.get('problems', [])) if st.get('status') != 'MISSING' else 'MISSING'}")
    print('ALL FAMILIES REPLAYED, COMPLETE AND WITHIN THEIR BOUNDS:', 'YES' if allok else 'NO')

"""AUDIT-7 g1: engine KQ with the winding fallback (kq_run.py, ORDER = 'wind_first').

Input: a kq_run output (final JSON, JSON.gz, or the partial .records.jsonl; only completed records are used for the
partial file).  For every record: own replay (split rule with the KQ scale (1, 1, pi/2, 1)), exact consumption of the
tree code, empty stack (done) or stack == 'left', fails == '!' leaves in order, leaf volume sum, nleaf, number of
winding leaves == len(wL) == nwind (when present).  Leaf codes: 'a'..'l' KQ count total 0..11, '1'..'6' winding W.
Then:
  (i)   winding leaves (all, or a sample of NW): re-certify with wind_cert.wind_box on nb_cert.box_coeffs(1, box)
        (the f-chart polynomial, 07:50 enclosure) with the stored L and the header's wind parameters: must give ok and
        the same W; if that fails, nb_cert2.wind_box2 (since 08:40) is tried and counted separately;
  (ii)  count leaves (sample of NC): kq_run_v2_backup.certify (07:50 code) must give ok and total <= recorded bound;
        if that fails, the current kq_run.certify (nb_cert2 away count) is tried and counted separately;
  (iii) 100-digit spot checks with MY OWN formula P(f) = mu^2 k^2 sin^2 f sigma^2 - D^2 B^2 (R-normalised,
        mu = sqrt((1-t)/t)), a positive multiple of the f-chart polynomial: at NPT points of NSPOT sampled leaves,
        N_mult(P) (zeros on |z| = 1 with multiplicity) <= min(bound, 10) and even; for winding leaves also
        #{zeros of the self-inversive Q(z) = sum c_k z^(k+6) in |z| < e^L} == W, and the relative distance of the
        nearest zero to the test circle (margin).
usage: python g1_kq_wind.py NW NC NSPOT NPT SEED FILE [FILE ...]      (NW/NC < 0: all leaves)"""
import os, sys, json, gzip, math, random
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.normpath(os.path.join(HERE, '..', '..', 'code'))
NW, NC, NSPOT, NPT, SEED = [int(v) for v in sys.argv[1:6]]
files = [os.path.abspath(f) for f in sys.argv[6:]]
sys.path.insert(0, CODE); sys.path.insert(0, os.path.join(CODE, 'c3'))
os.chdir(CODE)                      # kq_run inserts the relative path 'c3'
import importlib
import splitrule
import nb_cert
import wind_cert
kq_run = importlib.import_module('kq_run')          # current code (since 08:40: nb_cert2 enclosure)
kq_run.init()
if hasattr(kq_run, 'DEEP_ALWAYS'):
    kq_run.DEEP_ALWAYS = True                        # round 8: two-tier counts; re-certification allows the deep pass
try:                                                 # the audited 07:50 code (nb_cert / c3_away), if still present
    kq_old = importlib.import_module('kq_run_v2_backup'); kq_old.init()
except Exception:
    kq_old = None
try:
    import nb_cert2
except Exception:
    nb_cert2 = None
mp.mp.dps = 100
random.seed(SEED)
KQ_SCALE = (1.0, 1.0, math.pi / 2, 1.0)


def load(fn):
    if fn.endswith('.jsonl'):
        recs = []
        for line in open(fn):
            try:
                r = json.loads(line)
            except Exception:
                continue
            if r.get('done'):
                recs.append(r)
        hdr = dict(engine='KQ', scale=list(KQ_SCALE), wind_m0=kq_run.WIND_M0, wind_rounds=kq_run.WIND_ROUNDS,
                   wind_cap=kq_run.WIND_CAP)
        import os as _os
        _r = _os.environ.get('AUDIT_SPLIT_RULE'); _s = _os.environ.get('AUDIT_SCALE')
        if _r is None and _s is None:
            print('   WARNING partial file without header: leaves rebuilt with the default rule and driver scale; '
                  'valid only for runs without a split_rule and with the driver scale (set AUDIT_SPLIT_RULE / AUDIT_SCALE)', flush=True)
        if _r:
            hdr['split_rule'] = _r
        if _s:
            hdr['scale'] = [float(v) for v in _s.split(',')]
        return hdr, recs, True
    p = fn if os.path.exists(fn) else fn + '.gz'
    D = json.load(gzip.open(p, 'rt') if p.endswith('.gz') else open(p))
    return D['header'], D['results'], False


SPLIT_RULE = None                                    # header split_rule of the file being replayed (round 9)


def split(bx, S):
    if SPLIT_RULE:
        return splitrule.split_box(bx, S, SPLIT_RULE)
    w = [(bx[2 * i + 1] - bx[2 * i]) / S[i] for i in range(4)]
    i = int(np.argmax(w)); mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)


def vol(b):
    return float(np.prod([b[2 * i + 1] - b[2 * i] for i in range(4)]))


def coeffs(e, l, w, t):
    e, l, w, t = [mp.mpf(v) for v in (e, l, w, t)]
    pz, qz = l * mp.sin(w), l * mp.cos(w)
    mu = mp.sqrt((1 - t) / t)
    N = 24; vals = []
    for n in range(N):
        f = 2 * mp.pi * n / N
        zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
        B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta); D = 1 + e * mp.cos(f); k = 2 * e * (1 + e)
        vals.append(mu ** 2 * k ** 2 * mp.sin(f) ** 2 * (1 - zeta ** 2) - D ** 2 * B ** 2)
    return [mp.fsum(vals[n] * mp.expj(-kk * 2 * mp.pi * n / N) for n in range(N)) / N for kk in range(-6, 7)]


def zeros(c):
    """zeros of Q(z) = sum_{k=-6..6} c_k z^(k+6): list of finite zeros (zeros at 0 included), number at infinity"""
    sc = max(abs(v) for v in c)
    lo = 0
    while abs(c[lo]) < mp.mpf(10) ** -80 * sc:
        lo += 1
    hi = 12
    while abs(c[hi]) < mp.mpf(10) ** -80 * sc:
        hi -= 1
    z = list(mp.polyroots(c[lo:hi + 1][::-1], maxsteps=800, extraprec=600)) if hi > lo else []
    return z + [mp.mpf(0)] * lo, 12 - hi


for fn in files:
    if not (os.path.exists(fn) or os.path.exists(fn + '.gz')):
        print(os.path.basename(fn), 'not present', flush=True); continue
    H, recs, partial = load(fn)
    SPLIT_RULE = H.get('split_rule')
    KQF = H.get('engine') == 'KQF'
    IS_KQS = os.path.basename(fn).split('.')[0] == 'KQS'     # round 11: the --c3prune shell run
    TXF = H.get('T_X', 0.06)                                  # round 11: outer count from the header T_X (KQF2: 3 SF2)
    if KQF and not (TXF <= 3 * (0.010005 if os.path.basename(fn).startswith('KQF2') else 0.02)):
        print('   T_X of the header too large for the gap charts:', TXF, flush=True)
    kq_run.T_X = TXF; kq_run.FACE_S = H.get('FACE_S', 0.02)
    S = H['scale']
    m0, rnds, cap = H.get('wind_m0', 64), H.get('wind_rounds', 3), H.get('wind_cap', 1024)
    leaves = []; bad = []; nfail = 0; npr = [0]
    for r in recs:
        root = tuple(r['box']); stack = [root]; fl = [tuple(f[0]) for f in r.get('fails', [])]; fi = 0; wi = 0
        wl = r.get('wL', []); v = 0.0; nl = 0; ndig = 0
        for ch in r['tree']:
            if not stack:
                bad.append(('underflow', root)); break
            bx = stack.pop()
            if ch == '+':
                stack.extend(split(bx, S)); continue
            v += vol(bx)
            if ch == '!':
                if fi >= len(fl) or fl[fi] != bx:
                    bad.append(('fail mismatch', root))
                fi += 1
            elif 'a' <= ch <= 'l':
                leaves.append((bx, ord(ch) - ord('a'), None, None)); nl += 1
            elif ch == 'p' and not KQF:               # round 9: face neighbourhood, verified exactly
                nl += 1; npr[0] += 1                  # round 11: radius by role (SF2 for KQS, 0.02 otherwise)
                t0_ = mp.mpf(bx[6]); fs_ = mp.mpf(0.010005 if IS_KQS else 0.02)
                if not (1 - mp.mpf(bx[0]) <= fs_ ** 2 and mp.sqrt((1 - t0_) / t0_) <= fs_ ** 2):
                    bad.append(('p not in the face neighbourhood', bx))
            elif ch == 'o' and (KQF or IS_KQS):       # round 9: C3 cell directions, verified exactly (round 11: + KQS)
                nl += 1; npr[0] += 1
                if not (mp.mpf(bx[2]) >= mp.mpf(0.95) and mp.mpf(bx[4]) >= mp.mpf(1.413716694115407)):
                    bad.append(('o not in the C3 cell', bx))
            elif '1' <= ch <= '6':
                if wi >= len(wl):
                    bad.append(('missing wL', root)); L = None
                else:
                    L = wl[wi]
                wi += 1; ndig += 1
                leaves.append((bx, 12 - 2 * int(ch), int(ch), L)); nl += 1
            else:
                bad.append(('bad code ' + ch, root))
        if r['done'] and stack:
            bad.append(('stack left in done record', root))
        if not r['done'] and stack != [tuple(b) for b in r.get('left', [])]:
            bad.append(('left mismatch', root))
        v += sum(vol(b) for b in r.get('left', []))
        if fi != len(fl) or wi != len(wl) or nl != r['nleaf'] or ('nwind' in r and r['nwind'] != ndig):
            bad.append(('count mismatch', root, fi, len(fl), wi, len(wl), nl, r['nleaf'], r.get('nwind'), ndig))
        if abs(v - vol(root)) > 1e-9 * vol(root):
            bad.append(('volume', root))
        nfail += len(fl)
    wleaves = [l for l in leaves if l[2] is not None]
    cleaves = [l for l in leaves if l[2] is None]
    hist = {}
    for l in leaves:
        k = ('w%d' % l[2]) if l[2] is not None else l[1]
        hist[k] = hist.get(k, 0) + 1
    print(f"{os.path.basename(fn)}: records {len(recs)}{' (completed, partial file)' if partial else ''}, leaves {len(leaves)}"
          f" (winding {len(wleaves)}), fails {nfail}, replay problems {len(bad)} {bad[:3]}; "
          f"histogram {dict(sorted(hist.items(), key=lambda kv: str(kv[0])))}", flush=True)
    # (i) winding re-certification
    ws = wleaves if NW < 0 else random.sample(wleaves, min(NW, len(wleaves)))
    wbad = 0; wnew = 0
    for bx, bound, W, L in ws:
        ok, W2 = wind_cert.wind_box(*nb_cert.box_coeffs(1, *bx[:6]), bx[6], bx[7], L, m0, rnds, cap)   # 07:50 enclosure
        if (not ok or W2 != W) and nb_cert2 is not None:
            ok, W2 = nb_cert2.wind_box2(*nb_cert2.box_coeffs2(*bx[:6]), bx[6], bx[7], L, m0, rnds, cap)
            wnew += ok and W2 == W
        if not ok or W2 != W:
            wbad += 1; print('   WINDING RECERT PROBLEM', bx, W, L, ok, W2, flush=True)
    # (ii) count re-certification
    cs = cleaves if NC < 0 else random.sample(cleaves, min(NC, len(cleaves)))
    cbad = 0; cnew = 0
    for bx, bound, _, _ in cs:
        if KQF:
            kq_run.FACE_MODE = True
            ok, tot = kq_run.certify_outer(bx)
            kq_run.FACE_MODE = False
        else:
            ok, tot = (kq_old or kq_run).certify(bx)
            if (not ok or tot > bound) and kq_old is not None:
                ok, tot = kq_run.certify(bx)
                cnew += ok and tot <= bound
        if not ok or tot > bound:
            cbad += 1; print('   COUNT RECERT PROBLEM', bx, bound, ok, tot, flush=True)
    print(f"    re-certified winding leaves {len(ws)} (problems {wbad}; certified only by the nb_cert2 enclosure: {wnew}); "
          f"count leaves {len(cs)} (problems {cbad}; only by the current kq_run: {cnew})", flush=True)
    # (iii) 100-digit spot checks (winding leaves first, then count leaves)
    sw = random.sample(wleaves, min(NSPOT, len(wleaves)))
    sc_ = random.sample(cleaves, min(NSPOT, len(cleaves)))
    viol = 0; odd = 0; npts = 0; maxN = 0; minmarg = mp.inf; wchecked = 0
    for bx, bound, W, L in sw + sc_:
        for k in range(NPT):
            x = [0.5 * (bx[2 * i] + bx[2 * i + 1]) if k == 0 else random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(4)]
            if x[0] >= 1.0 or x[3] >= 1.0:          # e < 1 and t < 1 (physical)
                x[0] = min(x[0], 1 - 1e-9); x[3] = min(x[3], 1 - 1e-12)
            z, ninf = zeros(coeffs(*x))
            Nm = sum(1 for r in z if r != 0 and abs(mp.log(abs(r))) < mp.mpf(10) ** -30)
            if KQF:                                   # outer count: zeros with |f - pi| >= 2 atan(T_X)
                Nm = sum(1 for r in z if r != 0 and abs(mp.log(abs(r))) < mp.mpf(10) ** -30 and
                         mp.pi - abs(mp.arg(r)) >= 2 * mp.atan(mp.mpf(TXF)) - mp.mpf(10) ** -12)
            npts += 1; maxN = max(maxN, Nm)
            if Nm % 2 and not KQF:
                odd += 1
            if Nm > min(bound, 10):
                viol += 1; print('   SPOT VIOLATION', bx, x, Nm, bound, flush=True)
            if W is not None:
                rho = mp.e ** mp.mpf(L)
                inside = sum(1 for r in z if abs(r) < rho)
                marg = min(abs(abs(r) - rho) / rho for r in z) if z else mp.inf
                minmarg = min(minmarg, marg); wchecked += 1
                if inside != W:
                    viol += 1; print('   WINDING SPOT VIOLATION', bx, x, 'inside', inside, 'W', W, 'L', L, flush=True)
    print(f"    100-digit: {len(sw)} winding + {len(sc_)} count leaves, {npts} points ({wchecked} winding); violations {viol};"
          f" odd N_mult {odd}; max N_mult {maxN}; min relative distance of a zero to |z| = e^L {mp.nstr(minmarg, 3)}",
          flush=True)

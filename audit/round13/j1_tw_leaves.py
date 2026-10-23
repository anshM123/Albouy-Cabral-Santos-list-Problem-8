"""AUDIT-13 j1: tau-winding (TW, code '7') leaves of the running levels, checked independently.
For each sampled TW leaf (box bx, stored contour [which, c_re, c_im, r]):
 (1) replay: '7' count == len(tw) per record, every contour with which in {0,1} (c3: 0) and c_im > r > 0;
 (2) tau_wind.check_tw / check_tw_chart re-certifies the stored contour;
 (3) my own 60-digit zeros (my P, FFT -> degree-12 polynomial in z = e^{if}, Moebius map to th or tp; c3 Ca: th = nu*y):
     at the box centre and NPT random points the number of zeros strictly inside |w - c| < r is >= 1 and the same at
     every point, and N_P (zeros on |z| = 1) <= 10;
 (4) enclosure test of tau_wind.arc_radius at complex points: for NARC arcs (random centre angle, half-width pi/64 or
     pi/16), at parameter points x of the leaf's parameter box (corners + random) and w in the disc |w - z_j| <= D (on
     its boundary and inside): |F(w; x) - b0| <= R with my own F (normalisation checked against the engine's cmid).
usage: python j1_tw_leaves.py RECORDS.jsonl|RUN.json KIND RULE S0 S1 S2 S3 NLEAF NPT NARC SEED
       KIND: kq or Ca;  RULE: split rule or - (plain)"""
import sys, os, json, math, random, cmath
import numpy as np
import mpmath as mp
CODE = r'<repo>\code'
sys.path.insert(0, CODE); sys.path.insert(0, os.path.join(CODE, 'c3')); os.chdir(CODE)
import splitrule
import tau_wind
fn, KIND, RULE = sys.argv[1], sys.argv[2], (None if sys.argv[3] == '-' else sys.argv[3])
S = [float(v) for v in sys.argv[4:8]] if RULE != '@' else None
if RULE == '@':                                   # AUDIT-15: take split rule and scale from the run file header
    _H = json.load(open(fn))['header']; RULE = _H.get('split_rule'); S = _H['scale']
NLEAF, NPT, NARC, SEED = int(sys.argv[8]), int(sys.argv[9]), int(sys.argv[10]), int(sys.argv[11])
random.seed(SEED)
mp.mp.dps = 60


def load(fn):
    if fn.endswith('.json'):
        D = json.load(open(fn))
        return D['results'], 0
    keep = {}; bad = 0
    for line in open(fn):
        try:
            r = json.loads(line)
        except Exception:
            bad += 1
            continue
        b = tuple(r['box'])
        if r['done'] or b not in keep:
            keep[b] = r
    return list(keep.values()), bad


def split(bx):
    if RULE is not None:
        return splitrule.split_box(bx, S, RULE)
    w = [(bx[2 * i + 1] - bx[2 * i]) / S[i] for i in range(4)]
    i = int(np.argmax(w)); mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)


recs, nbad = load(fn)
tw_leaves = []; probs = []
for r in recs:
    stack = [tuple(r['box'])]; tw = r.get('tw', []); k = 0
    for ch in r['tree']:
        if not stack:
            probs.append('underflow'); break
        bx = stack.pop()
        if ch == '+':
            stack.extend(split(bx)); continue
        if ch == '7':
            if k < len(tw):
                t = tw[k]
                if not (len(t) == 4 and t[0] in ((0, 1) if KIND == 'kq' else (0,)) and t[2] > t[3] > 0):
                    probs.append('bad contour %r' % (t,))
                tw_leaves.append((bx, t))
            k += 1
    if k != len(tw):
        probs.append('number of 7 leaves %d != len(tw) %d' % (k, len(tw)))
    if not r['done'] and stack != [tuple(b) for b in r.get('left', [])]:
        probs.append('left mismatch')
print(f'{os.path.basename(fn)}: records {len(recs)} (unparsable lines {nbad}), TW leaves {len(tw_leaves)}, '
      f'replay problems {len(probs)} {probs[:3]}', flush=True)
if os.environ.get('J1_WHICH') is not None:          # restrict to contours of one chart (0: Q, 1: G)
    tw_leaves = [l for l in tw_leaves if l[1][0] == int(os.environ['J1_WHICH'])]
    print('   restricted to which =', os.environ['J1_WHICH'], ':', len(tw_leaves), 'leaves', flush=True)
sample = random.sample(tw_leaves, min(NLEAF, len(tw_leaves)))


def Pf(e, pz, qz, mu, f):
    zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
    B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta); D = 1 + e * mp.cos(f); k = 2 * e * (1 + e)
    return mu ** 2 * k ** 2 * mp.sin(f) ** 2 * (1 - zeta ** 2) - D ** 2 * B ** 2


def fourier(e, pz, qz, mu):
    N = 32
    vals = [Pf(e, pz, qz, mu, 2 * mp.pi * n / N) for n in range(N)]
    return {k: mp.fsum(vals[n] * mp.expj(-k * 2 * mp.pi * n / N) for n in range(N)) / N for k in range(-6, 7)}


def F_th(c, th):                       # Q(th) = P(f)(1 + th^2)^6, th = tan((f - pi)/2)
    return mp.fsum(c[k] * (-1) ** k * (1 + 1j * th) ** (6 + k) * (1 - 1j * th) ** (6 - k) for k in range(-6, 7))


def F_tp(c, tp):                       # G(tp) = P(f)(1 + tp^2)^6, tp = tan(f/2)
    return mp.fsum(c[k] * (1 + 1j * tp) ** (6 + k) * (1 - 1j * tp) ** (6 - k) for k in range(-6, 7))


def zroots(c):
    co = [c[k] for k in range(6, -7, -1)]          # z^12 ... z^0 coefficients of z^6 P
    sc = max(abs(v) for v in co)
    while abs(co[0]) < mp.mpf(10) ** -50 * sc:
        co = co[1:]
    return mp.polyroots(co, maxsteps=3000, extraprec=2000)


def phys_kq(x):                       # parameter-box point (eps, Pz, Qz, mu) -> my P arguments
    return 1 - x[0], x[1], x[2], x[3]


def phys_ca(x):                       # chart Ca: (nu, K, E, M) -> qz = nu, kap = nu^2 K, eps = nu^2 E, mu = nu^3 M
    nu = x[0]
    return 1 - nu ** 2 * x[2], 1 - nu ** 2 * x[1], nu, nu ** 3 * x[3]


def chart_value(c, w, which, nu):
    if KIND == 'kq':
        return F_th(c, w) if which == 0 else F_tp(c, w)
    return F_th(c, nu * w) / nu ** 10                # Ca: G(y) = Q(nu y)/nu^10


def zeros_chart(c, which, nu):
    out = []
    for z in zroots(c):
        if abs(z - 1) < mp.mpf(10) ** -40 or abs(z + 1) < mp.mpf(10) ** -40:
            out.append((z, None)); continue
        if KIND == 'kq':
            w = -1j * (z + 1) / (z - 1) if which == 0 else -1j * (z - 1) / (z + 1)
        else:
            w = (-1j * (z + 1) / (z - 1)) / nu
        out.append((z, w))
    return out


tot_bad = 0; worst_ratio = 0.0; worst_norm = 0.0; min_margin = mp.inf; nchk = 0; inside_counts = []
for bx, t in sample:
    which, cre, cim, rr = t
    c0 = complex(cre, cim)
    ok_tw = tau_wind.check_tw(bx, which, cre, cim, rr) if KIND == 'kq' else tau_wind.check_tw_chart(KIND, bx, cre, cim, rr)
    pbox = tau_wind.params(bx) if KIND == 'kq' else tuple(bx)
    phys = phys_kq if KIND == 'kq' else phys_ca
    pts = [[0.5 * (pbox[2 * i] + pbox[2 * i + 1]) for i in range(4)]] + \
          [[random.uniform(pbox[2 * i], pbox[2 * i + 1]) for i in range(4)] for _ in range(NPT)]
    counts = []; nps = []
    for x in pts:
        xm = [mp.mpf(v) for v in x]
        if KIND != 'kq' and xm[0] == 0:
            continue
        c = fourier(*phys(xm))
        zs = zeros_chart(c, which, xm[0])
        NP = sum(1 for z, w in zs if abs(abs(z) - 1) < mp.mpf(10) ** -20)
        ins = [w for z, w in zs if w is not None and abs(w - c0) < rr]
        counts.append(len(ins)); nps.append(NP)
        for z, w in zs:
            if w is not None:
                min_margin = min(min_margin, abs(abs(w - c0) - rr) / rr)
    good = ok_tw and min(counts) >= 1 and len(set(counts)) == 1 and max(nps) <= 10
    tot_bad += not good
    inside_counts.append(counts[0])
    if not good:
        print('   PROBLEM', bx, t, 'check_tw', ok_tw, 'inside', counts, 'N_P', nps, flush=True)
    # (4) enclosure test
    eng = tau_wind.engine(which) if KIND == 'kq' else __import__('nbeng_c3_%s_2' % KIND)
    data = eng.box_coeffs2(*pbox)
    xc = [mp.mpf(0.5 * (pbox[2 * i] + pbox[2 * i + 1])) for i in range(4)]
    cc = fourier(*phys(xc))
    for thz in (0.0, 1.0, 2.5, 4.0):                 # points ON the contour (F bounded away from 0 there)
        zt = c0 + rr * cmath.exp(1j * thz)
        ref = chart_value(cc, mp.mpc(zt.real, zt.imag), which, xc[0])
        eng_v = complex(tau_wind.centre_value(data[0], zt))
        worst_norm = max(worst_norm, float(abs(ref - eng_v) / max(abs(ref), mp.mpf(10) ** -300)))
    corners = [[pbox[2 * i + ((m >> i) & 1)] for i in range(4)] for m in range(16)]
    xs = random.sample(corners, 3) + [[random.uniform(pbox[2 * i], pbox[2 * i + 1]) for i in range(4)] for _ in range(2)]
    for a_ in range(NARC):
        th0 = random.uniform(0, 2 * math.pi); h = random.choice((math.pi / 64, math.pi / 16))
        zj = c0 + rr * cmath.exp(1j * th0); D = rr * h * (1 + 1e-12) + 1e-12 * (abs(c0) + rr)
        v0, b0, R = tau_wind.arc_radius(*data, zj, D, tau_wind.BINOM)
        for x in xs:
            xm = [mp.mpf(v) for v in x]
            if KIND != 'kq' and xm[0] == 0:
                continue
            cx = fourier(*phys(xm))
            for q in range(4):
                ang = random.uniform(0, 2 * math.pi); rad = D if q < 2 else D * random.random()
                w = zj + rad * cmath.exp(1j * ang)
                Fw = chart_value(cx, mp.mpc(w.real, w.imag), which, xm[0])
                ratio = float(abs(Fw - mp.mpc(b0.real, b0.imag)) / R)
                worst_ratio = max(worst_ratio, ratio); nchk += 1
                if ratio > 1:
                    tot_bad += 1
                    print('   ENCLOSURE VIOLATION', bx, t, 'arc', th0, h, 'x', x, 'w', w, 'ratio', ratio, flush=True)
print(f'   sampled {len(sample)} TW leaves x {NPT + 1} points: problems {tot_bad}; zeros inside the contour per leaf '
      f'{sorted(set(inside_counts))}; min relative distance of a zero to the contour {mp.nstr(min_margin, 3)}; '
      f'engine vs my F at x_c (normalisation) max rel diff {worst_norm:.2e}; enclosure checks {nchk}, '
      f'max |F(w;x) - b0|/R = {worst_ratio:.3f}', flush=True)

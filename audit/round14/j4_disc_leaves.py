"""AUDIT-14 j4: disc-winding outer-count leaves (code '8', KQF-engine files), checked independently.
For each sampled '8' leaf (box bx; contour |th| = T_X of the file header, centre 0):
 (1) tau_wind.check_disc(bx, T_X) re-certifies;
 (2) my own 60-digit zeros of Q (my P, FFT -> degree-12 polynomial in z = e^{if}, th = -i(z+1)/(z-1)) at the centre and
     NPT random physical points: W = #zeros with |th| < T_X is >= 6 and the same at every point; the outer real count
     (zeros on |z| = 1 with |th| >= T_X, f = 0 included) is <= 12 - W <= 6; margin min | |th| - T_X | / T_X;
 (3) enclosure test of tau_wind.arc_radius on the disc contour (c = 0, r = T_X): |F(w; x) - b0| <= R at corners and
     random points x of the leaf's (eps, Pz, Qz, mu) box and w on / inside the arc discs.
Also: the replay of every record (no '8' outside KQF engine files; '8' count consistent with hist['d']).
usage: python j4_disc_leaves.py RUN.json NLEAF NPT NARC SEED"""
import sys, os, json, math, random, cmath
import numpy as np
import mpmath as mp
CODE = r'<repo>\code'
sys.path.insert(0, CODE); sys.path.insert(0, os.path.join(CODE, 'c3')); os.chdir(CODE)
import splitrule
import tau_wind
fn = sys.argv[1]; NLEAF, NPT, NARC, SEED = [int(v) for v in sys.argv[2:6]]
random.seed(SEED)
mp.mp.dps = 60
D = json.load(open(fn)); H = D['header']; S = H['scale']; RULE = H.get('split_rule'); TX = H['T_X']
assert H['engine'] == 'KQF', H['engine']


def split(bx):
    if RULE is not None:
        return splitrule.split_box(bx, S, RULE)
    w = [(bx[2 * i + 1] - bx[2 * i]) / S[i] for i in range(4)]
    i = int(np.argmax(w)); mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)


dl = []; probs = []; nd_tot = 0
for r in D['results']:
    stack = [tuple(r['box'])]; nd = 0
    for ch in r['tree']:
        bx = stack.pop()
        if ch == '+':
            stack.extend(split(bx)); continue
        if ch == '8':
            dl.append(bx); nd += 1
    if nd != r.get('hist', {}).get('d', 0):
        probs.append(('8 count != hist[d]', nd, r.get('hist', {}).get('d')))
    nd_tot += nd
print(f"{os.path.basename(fn)}: header T_X {TX!r}, FACE_S {H.get('FACE_S')!r}, disc flag {H.get('disc')!r}; records "
      f"{len(D['results'])}, '8' leaves {nd_tot}, problems {len(probs)} {probs[:3]}", flush=True)
sample = random.sample(dl, min(NLEAF, len(dl)))


def Pf(e, pz, qz, mu, f):
    zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
    B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta); Dd = 1 + e * mp.cos(f); k = 2 * e * (1 + e)
    return mu ** 2 * k ** 2 * mp.sin(f) ** 2 * (1 - zeta ** 2) - Dd ** 2 * B ** 2


def fourier(e, pz, qz, mu):
    N = 32
    vals = [Pf(e, pz, qz, mu, 2 * mp.pi * n / N) for n in range(N)]
    return {k: mp.fsum(vals[n] * mp.expj(-k * 2 * mp.pi * n / N) for n in range(N)) / N for k in range(-6, 7)}


def F_th(c, th):
    return mp.fsum(c[k] * (-1) ** k * (1 + 1j * th) ** (6 + k) * (1 - 1j * th) ** (6 - k) for k in range(-6, 7))


def zroots(c):
    co = [c[k] for k in range(6, -7, -1)]
    sc = max(abs(v) for v in co)
    while abs(co[0]) < mp.mpf(10) ** -50 * sc:
        co = co[1:]
    return mp.polyroots(co, maxsteps=3000, extraprec=2000), 12 - (len(co) - 1)   # (roots, zeros at z = infinity)


bad = 0; worst = 0.0; nchk = 0; margin = mp.inf; Wset = set(); maxout = 0
txm = mp.mpf(TX)
for bx in sample:
    ok = tau_wind.check_disc(bx, TX)
    pb = tau_wind.params(bx)
    pts = [[0.5 * (bx[2 * i] + bx[2 * i + 1]) for i in range(4)]] + \
          [[random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(4)] for _ in range(NPT)]
    Ws = set(); outs = []
    for x in pts:
        x[0] = min(x[0], 1 - 1e-12); x[3] = min(x[3], 1 - 1e-15)             # physical: e < 1, t < 1
        e, l, w, t = [mp.mpf(v) for v in x]
        c = fourier(e, l * mp.sin(w), l * mp.cos(w), mp.sqrt((1 - t) / t))
        zs, ninf = zroots(c)
        Win = 0; outer = 0
        for z in zs:
            if abs(z - 1) < mp.mpf(10) ** -30:           # f = 0 (th = infinity): outer, real
                outer += 1; continue
            th = -1j * (z + 1) / (z - 1)
            margin = min(margin, abs(abs(th) - txm) / txm)
            if abs(th) < txm:
                Win += 1
            elif abs(abs(z) - 1) < mp.mpf(10) ** -20:
                outer += 1
        Ws.add(Win); outs.append(outer); maxout = max(maxout, outer)
    Wset |= Ws
    good = ok and len(Ws) == 1 and min(Ws) >= 6 and max(outs) <= 12 - min(Ws)
    if not good:
        bad += 1
        print('   PROBLEM', bx, 'check_disc', ok, 'W', Ws, 'outer', outs, flush=True)
    # (3) enclosure on the disc contour
    data = tau_wind.engine(0).box_coeffs2(*pb)
    corners = [[pb[2 * i + ((m >> i) & 1)] for i in range(4)] for m in range(16)]
    xs = random.sample(corners, 3) + [[random.uniform(pb[2 * i], pb[2 * i + 1]) for i in range(4)] for _ in range(2)]
    for a_ in range(NARC):
        th0 = random.uniform(0, 2 * math.pi); h = random.choice((math.pi / 64, math.pi / 16))
        zj = TX * cmath.exp(1j * th0); Dd = TX * h * (1 + 1e-12) + 1e-12 * TX
        v0, b0, R = tau_wind.arc_radius(*data, zj, Dd, tau_wind.BINOM)
        for x in xs:
            cx = fourier(1 - mp.mpf(x[0]), mp.mpf(x[1]), mp.mpf(x[2]), mp.mpf(x[3]))
            for q in range(4):
                wv = zj + (Dd if q < 2 else Dd * random.random()) * cmath.exp(1j * random.uniform(0, 2 * math.pi))
                ratio = float(abs(F_th(cx, mp.mpc(wv.real, wv.imag)) - mp.mpc(b0.real, b0.imag)) / R)
                worst = max(worst, ratio); nchk += 1
                if ratio > 1:
                    bad += 1; print('   ENCLOSURE VIOLATION', bx, x, wv, ratio, flush=True)
print(f"   sampled {len(sample)} '8' leaves x {NPT + 1} points: problems {bad}; W (zeros in |th| < T_X) {sorted(Wset)}; "
      f"max outer real count {maxout}; min relative distance of a zero to |th| = T_X {mp.nstr(margin, 3)}; enclosure "
      f"checks {nchk}, max |F(w;x) - b0|/R = {worst:.12f}", flush=True)

"""AUDIT-13 j2: (a) the counting identity N_P = 12 - #nonreal zeros of Q (and of G), with degree drops;
(b) negative tests: tau_wind must NOT certify a box that contains a point with 12 real zeros (e = 1, mu = 0,
    l >= l*(w)), nor a contour that encloses no zero, nor one that passes through a zero;
(c) adversarial enclosure test on enlarged boxes (larger parameter variation T) and large discs.
Q is built independently here: its th-coefficients by exact interpolation of P(f(th))(1 + th^2)^6 at 13 nodes."""
import sys, os, math, random, cmath
import numpy as np
import mpmath as mp
CODE = r'<repo>\code'
sys.path.insert(0, CODE); sys.path.insert(0, os.path.join(CODE, 'c3')); os.chdir(CODE)
import tau_wind
mp.mp.dps = 60
random.seed(13)


def Pf(e, pz, qz, mu, f):
    zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
    B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta); D = 1 + e * mp.cos(f); k = 2 * e * (1 + e)
    return mu ** 2 * k ** 2 * mp.sin(f) ** 2 * (1 - zeta ** 2) - D ** 2 * B ** 2


def coeffs_interp(e, pz, qz, mu, which):
    """coefficients (low to high) of Q(th) (which 0) or G(tp) (which 1) by interpolation at 13 real nodes"""
    nodes = [mp.mpf(j - 6) / 3 for j in range(13)]
    vals = []
    for s in nodes:
        f = mp.pi + 2 * mp.atan(s) if which == 0 else 2 * mp.atan(s)
        vals.append(Pf(e, pz, qz, mu, f) * (1 + s * s) ** 6)
    V = mp.matrix([[s ** k for k in range(13)] for s in nodes])
    return list(mp.lu_solve(V, mp.matrix(vals)))


def poly_roots(c):
    c = list(c)
    sc = max(abs(v) for v in c)
    while abs(c[-1]) < mp.mpf(10) ** -40 * sc:
        c.pop()
    return (mp.polyroots(c[::-1], maxsteps=3000, extraprec=2000) if len(c) > 1 else []), len(c) - 1


def N_P(e, pz, qz, mu, tol=mp.mpf(10) ** -15):
    N = 32
    vals = [Pf(e, pz, qz, mu, 2 * mp.pi * n / N) for n in range(N)]
    cf = [mp.fsum(vals[n] * mp.expj(-k * 2 * mp.pi * n / N) for n in range(N)) / N for k in range(-6, 7)]
    co = cf[::-1]
    sc = max(abs(v) for v in co)
    while abs(co[0]) < mp.mpf(10) ** -40 * sc:
        co = co[1:]
    zz = mp.polyroots(co, maxsteps=3000, extraprec=2000)
    return sum(1 for z in zz if abs(abs(z) - 1) < tol)


print('(a) counting identity  N_P = 12 - #nonreal(Q) = 12 - #nonreal(G)')
bad = 0; rows = 0
cases = []
for _ in range(25):                                      # random physical points
    e = random.uniform(0.0, 0.999); l = random.uniform(0, 1); w = random.uniform(0, math.pi / 2); mu = 10 ** random.uniform(-4, 0.5)
    cases.append(('random', e, l * math.sin(w), l * math.cos(w), mu))
for _ in range(8):                                       # Qz = 0: zero of P at f = 0 (deg Q <= 10)
    e = random.uniform(0.5, 0.999); l = random.uniform(0.3, 1); mu = 10 ** random.uniform(-4, -1)
    cases.append(('Qz=0', e, l, 0.0, mu))
for _ in range(8):                                       # e = 1: zero of P at f = pi (deg G drops)
    l = random.uniform(0.3, 1); w = random.uniform(0, math.pi / 2); mu = 10 ** random.uniform(-4, -1)
    cases.append(('e=1', 1.0, l * math.sin(w), l * math.cos(w), mu))
for name, e, pz, qz, mu in cases:
    e, pz, qz, mu = [mp.mpf(v) for v in (e, pz, qz, mu)]
    n = N_P(e, pz, qz, mu)
    out = []
    for which in (0, 1):
        r, deg = poly_roots(coeffs_interp(e, pz, qz, mu, which))
        nonreal = sum(1 for z in r if abs(mp.im(z)) > mp.mpf(10) ** -15 * max(1, abs(z)))
        nreal = deg - nonreal
        out.append((deg, nreal, nonreal))
        if 12 - nonreal != n or nreal + (12 - deg) != n:
            bad += 1
            print('   IDENTITY FAILS', name, [mp.nstr(v, 6) for v in (e, pz, qz, mu)], 'N_P', n, 'which', which, 'deg', deg,
                  'real', nreal, 'nonreal', nonreal, flush=True)
    rows += 1
    if name != 'random' and rows % 4 == 0:
        print(f'   {name}: N_P {n}; (deg, real, nonreal) of Q {out[0]}, of G {out[1]}', flush=True)
print(f'   {rows} points (25 random, 8 with Qz = 0, 8 with e = 1): identity failures {bad}', flush=True)

print('(b) negative tests')
# (b1) boxes containing a degenerate point e = 1, mu = 0 (t = 1) with l >= l*(w): 12 real zeros there
HP = math.pi / 2
neg = 0; neg_bad = 0
for (l0, l1, w0, w1) in ((0.99, 1.0, 1.0, 1.05), (0.995, 1.0, 0.6, 0.62), (0.999, 1.0, 1.3, 1.31), (0.98, 0.99, 0.2, 0.22)):
    for (e0, t0) in ((0.99999, 0.9999999999), (0.9999999, 0.999999999999)):
        bx = (e0, 1.0, l0, l1, w0, w1, t0, 1.0)
        n12 = N_P(mp.mpf(1), mp.mpf(l1) * mp.sin(mp.mpf(w1)), mp.mpf(l1) * mp.cos(mp.mpf(w1)), mp.mpf(0), tol=mp.mpf(10) ** -6)
        ok, which, cre, cim, r = tau_wind.certify_tw(bx, rounds=10, cap=65536, tries=6)
        neg += 1; neg_bad += bool(ok)
        print(f'   box {bx}: N_P at the corner (e = 1, mu = 0, l = {l1}, w = {w1}) = {n12}; certify_tw -> {ok}'
              + (f'  !!! contour {which, cre, cim, r}' if ok else ''), flush=True)
# (b2) wrong contours on a TW leaf box of KQK1.r7 (if a certified one is given): far from every zero -> W = 0
bx = (0.99, 0.990625, 0.9359375, 0.9375, 0.574322406984384, 0.579231145505618, 0.999860875, 0.99988075)
okc = tau_wind.certify_tw(bx)
print('   sample box', bx, 'certify_tw', okc, flush=True)
if okc[0]:
    which, cre, cim, r = okc[1:]
    far = tau_wind.check_tw(bx, which, cre + 10 * r, cim + 10 * r, r)
    big = tau_wind.check_tw(bx, which, cre, cim + 0.0, min(cim * 0.999, 50 * r))
    lower = tau_wind.check_tw(bx, which, cre, -cim, r)
    print(f'   contour moved away from the zero (W = 0 expected): {far};  contour in the lower half-plane: {lower}; '
          f'larger contour (may still be valid): {big}', flush=True)
    neg += 2; neg_bad += bool(far) + bool(lower)
print(f'   negative tests {neg}: wrongly certified {neg_bad}', flush=True)

print('(c) enclosure on enlarged boxes and large discs (my own F vs tau_wind.arc_radius)')


def F(c, w, which):
    if which == 0:
        return mp.fsum(c[k] * (-1) ** k * (1 + 1j * w) ** (6 + k) * (1 - 1j * w) ** (6 - k) for k in range(-6, 7))
    return mp.fsum(c[k] * (1 + 1j * w) ** (6 + k) * (1 - 1j * w) ** (6 - k) for k in range(-6, 7))


def fourier(e, pz, qz, mu):
    N = 32
    vals = [Pf(e, pz, qz, mu, 2 * mp.pi * n / N) for n in range(N)]
    return {k: mp.fsum(vals[n] * mp.expj(-k * 2 * mp.pi * n / N) for n in range(N)) / N for k in range(-6, 7)}


worst = 0.0; nchk = 0; viol = 0
for trial in range(12):
    # random parameter box in (eps, Pz, Qz, mu), widths up to 0.05, and a random disc
    c_ = [random.uniform(0, 0.05), random.uniform(0.3, 0.98), random.uniform(0.0, 0.7), random.uniform(0, 0.05)]
    wd = [10 ** random.uniform(-6, math.log10(0.05)) for _ in range(4)]
    pb = []
    for i in range(4):
        lo = max(c_[i] - wd[i], 0.0); pb += [lo, c_[i] + wd[i]]
    which = trial % 2
    data = tau_wind.engine(which).box_coeffs2(*pb)
    z0 = complex(random.uniform(-1, 1), random.uniform(0.05, 1.0)); D = random.choice((0.01, 0.1, 0.4))
    v0, b0, R = tau_wind.arc_radius(*data, z0, D, tau_wind.BINOM)
    corners = [[pb[2 * i + ((m >> i) & 1)] for i in range(4)] for m in range(16)]
    for x in random.sample(corners, 4) + [[random.uniform(pb[2 * i], pb[2 * i + 1]) for i in range(4)] for _ in range(2)]:
        cx = fourier(1 - mp.mpf(x[0]), mp.mpf(x[1]), mp.mpf(x[2]), mp.mpf(x[3]))
        for q in range(5):
            w = z0 + (D if q < 3 else D * random.random()) * cmath.exp(1j * random.uniform(0, 2 * math.pi))
            ratio = float(abs(F(cx, mp.mpc(w.real, w.imag), which) - mp.mpc(b0.real, b0.imag)) / R)
            worst = max(worst, ratio); nchk += 1
            if ratio > 1:
                viol += 1; print('   ENCLOSURE VIOLATION', pb, which, z0, D, x, w, ratio, flush=True)
print(f'   {nchk} checks on 12 random boxes (widths 1e-6 .. 0.05) and discs D in (0.01, 0.1, 0.4): violations {viol}, '
      f'max |F(w;x) - b0|/R = {worst:.3f}', flush=True)

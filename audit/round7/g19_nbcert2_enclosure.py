"""AUDIT-7 g19 (round-7 addendum A): independent end-to-end tests of the second-order f-chart enclosure nb_cert2.
  (H) Hessian in (e, l, w) at the box centre (box_coeffs2 hcen) against 50-digit numerical second derivatives of the
      exact coefficient c(e, l sin w, l cos w) (plain nbcoef_f polynomials evaluated in mpmath);
  (E) enclose_point2 containment, modes 1 (second-order bound only), 2 (mean-value only) and 0 (min, used): random
      f-chart boxes (e in [0.9, 1] and [0.2, 1], all l, w; widths 1e-5 .. 2e-1 per variable), t in {t0, t1}, random u and
      j = 0..7; exact P^(j)(u; x, t) = sum_k [(1-t) PhiF_k - t PsiF_k](x) d^j/du^j (cos, sin)(k u) at 50 digits at
      random x of the box (corners included) must lie in [vc - rad, vc + rad];
  (W) comp_bounds2 containment (winding components): weights fl(cosh(kL))(1+1e-12), fl(sinh(kL))(1+1e-12) as in
      labels2, phases 0, -pi/2, random L in [-3, -0.05], random arcs [u - h, u + h]: the exact S(theta; x, t) (exact
      coefficients, the given float weights) at random theta of the arc and x of the box must lie in [lo, hi].
Reports violations and the worst ratio |exact - centre value| / radius per mode.
usage: python g19_nbcert2_enclosure.py NBOX NPT SEED"""
import os, re, sys, types, random, math
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.normpath(os.path.join(HERE, '..', '..', 'code'))
sys.path.insert(0, CODE)
NBOX, NPT, SEED = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
random.seed(SEED)
mp.mp.dps = 50
import nb_cert2
J = 6


def load_plain(path):
    src = open(path).read().replace('from numba import njit', '')
    src = re.sub(r'@njit\([^)]*\)\n', '', src)
    mod = types.ModuleType(os.path.basename(path)); exec(compile(src, path, 'exec'), mod.__dict__)
    return mod


F1 = load_plain(os.path.join(CODE, 'nbcoef_f.py'))


def exact_coeffs(e, l, w):
    """[f][part][k] exact (mpf) coefficients of PhiF (f=0) and PsiF (f=1) at (e, l sin w, l cos w)"""
    out = np.zeros((4, 2, 2, 7), dtype=object); out[...] = mp.mpf(0)
    F1.coef_val(mp.mpf(e), mp.mpf(l) * mp.sin(mp.mpf(w)), mp.mpf(l) * mp.cos(mp.mpf(w)), out)
    return out[0]


def exact_P(cf, t, u, j):
    t = mp.mpf(t); u = mp.mpf(u); s = mp.mpf(0)
    for k in range(7):
        if j > 0 and k == 0:
            continue
        a = (1 - t) * cf[0, 0, k] - t * cf[1, 0, k]
        b = (1 - t) * cf[0, 1, k] - t * cf[1, 1, k]
        kj = mp.mpf(k) ** j
        s += kj * (a * mp.cos(k * u + j * mp.pi / 2) + b * mp.sin(k * u + j * mp.pi / 2))
    return s


def tcomb(data, t):
    cmid, crad, dcen, dwid, cmax, rr, hcen, rem = data
    c = (1 - t) * cmid[0] - t * cmid[1]
    cr = (1 - t) * crad[0] + t * crad[1] + 1e-15 * (np.abs(cmid[0]) + np.abs(cmid[1]))
    dc = (1 - t) * dcen[:, 0] - t * dcen[:, 1]
    dw = (1 - t) * dwid[:, 0] + t * dwid[:, 1] + 1e-15 * (np.abs(dcen[:, 0]) + np.abs(dcen[:, 1]))
    hc = (1 - t) * hcen[:, 0] - t * hcen[:, 1]
    rm = (1 - t) * rem[0] + t * rem[1]
    cm = (1 - t) * cmax[0] + t * cmax[1]
    qw = np.array([0.5 * rr[0] * rr[0], rr[0] * rr[1], rr[0] * rr[2], 0.5 * rr[1] * rr[1], rr[1] * rr[2], 0.5 * rr[2] * rr[2]])
    DA = np.abs(dcen[:, 0]) + np.abs(dcen[:, 1]); HA = np.abs(hcen[:, 0]) + np.abs(hcen[:, 1]); RA = rem[0] + rem[1]
    return c, cr, dc, dw, hc, rm, cm, qw, DA, HA, RA


def rand_box():
    lo_e = 0.9 if random.random() < 0.6 else 0.2
    wd = [10 ** random.uniform(-5, math.log10(0.2)) for _ in range(3)]
    e0 = random.uniform(lo_e, 1.0 - wd[0]); l0 = random.uniform(0.0, 1.0 - min(wd[1], 0.99)); w0 = random.uniform(0.0, max(0.0, math.pi / 2 - wd[2]))
    return (e0, e0 + wd[0], l0, min(l0 + wd[1], 1.0), w0, w0 + wd[2])


def box_points(bx, n):
    pts = [[bx[2 * i + ((m >> i) & 1)] for i in range(3)] for m in random.sample(range(8), 3)]
    pts += [[random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(3)] for _ in range(n)]
    return pts


# (H) Hessian check
worstH = mp.mpf(0)
for _ in range(max(5, NBOX // 20)):
    bx = rand_box()
    data = nb_cert2.box_coeffs2(*bx); hcen = data[6]
    xc = [mp.mpf(0.5 * (bx[2 * i] + bx[2 * i + 1])) for i in range(3)]
    h = mp.mpf(10) ** -12
    pairs = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]

    def cval(x):
        return exact_coeffs(*x)
    for m, (a, b) in enumerate(pairs):
        def shifted(sa, sb):
            x = list(xc); x[a] += sa * h; x[b] += sb * h
            return cval(x)
        num = (shifted(1, 1) - shifted(1, -1) - shifted(-1, 1) + shifted(-1, -1)) / (4 * h * h)
        for f in range(2):
            for pt in range(2):
                for k in range(7):
                    ref = num[f, pt, k]
                    worstH = max(worstH, abs(ref - hcen[m, f, pt, k]) / (1 + abs(ref)))
print(f'(H) Hessian (e, l, w) vs 50-digit numerical second derivatives: max rel err {mp.nstr(worstH, 3)}', flush=True)

# (E) enclose_point2 containment
stats = {0: [0, 0, 0.0], 1: [0, 0, 0.0], 2: [0, 0, 0.0]}      # mode -> [checks, violations, worst ratio]
near = []
for _ in range(NBOX):
    bx = rand_box()
    t0 = random.uniform(0.0, 1.0); t1 = min(1.0, t0 + 10 ** random.uniform(-6, -0.5))
    data = nb_cert2.box_coeffs2(*bx)
    pts = box_points(bx, NPT)
    cfs = [exact_coeffs(*x) for x in pts]
    for t in (t0, t1):
        c, cr, dc, dw, hc, rm, cm, qw, DA, HA, RA = tcomb(data, t)
        for _u in range(2):
            u = random.uniform(-3 * math.pi, 3 * math.pi)
            for mode in (0, 1, 2):
                vc, rad = nb_cert2.enclose_point2(c, cr, dc, dw, hc, rm, DA, HA, RA, qw, data[5], u, J, mode)
                for j in range(J + 2):
                    for cf in cfs:
                        ex = exact_P(cf, t, u, j)
                        d = abs(ex - mp.mpf(float(vc[j])))
                        st = stats[mode]; st[0] += 1
                        if d > rad[j]:
                            st[1] += 1
                            if st[1] <= 3:
                                print('   VIOLATION (E) mode', mode, bx, t, u, j, float(d), rad[j], flush=True)
                        if rad[j] > 0:
                            st[2] = max(st[2], float(d / rad[j]))
                            if float(d / rad[j]) > 0.999:
                                near.append((mode, j, float(1 - d / rad[j]), float(rad[j]), float(abs(mp.mpf(float(vc[j])))), bx))
if near:
    near.sort(key=lambda z: z[2])
    print('   tightest checks (mode, j, 1 - |exact - vc|/rad, rad, |vc|, box):', flush=True)
    for z in near[:6]:
        print('     ', z[0], z[1], '%.3e' % z[2], '%.3e' % z[3], '%.3e' % z[4], z[5], flush=True)
for mode, nm in ((1, 'second-order only'), (2, 'mean-value only'), (0, 'min (used)')):
    st = stats[mode]
    print(f'(E) enclose_point2 mode {mode} ({nm}): {st[0]} checks, violations {st[1]}, worst |exact - vc| / rad {st[2]:.3f}', flush=True)

# (W) comp_bounds2 containment
nW = 0; vW = 0; wr = 0.0
for _ in range(NBOX):
    bx = rand_box()
    t0 = random.uniform(0.0, 1.0); t1 = min(1.0, t0 + 10 ** random.uniform(-6, -0.5))
    data = nb_cert2.box_coeffs2(*bx)
    L = -random.uniform(0.05, 3.0)
    wR = np.array([math.cosh(k * L) * (1 + 1e-12) for k in range(7)]); wI = np.array([math.sinh(k * L) * (1 + 1e-12) for k in range(7)])
    pts = box_points(bx, NPT)
    cfs = [exact_coeffs(*x) for x in pts]
    for t in (t0, t1):
        c, cr, dc, dw, hc, rm, cm, qw, DA, HA, RA = tcomb(data, t)
        um = random.uniform(0, 2 * math.pi); hh = 10 ** random.uniform(-4, -0.7)
        for wgt, phase in ((wR, 0.0), (wI, -math.pi / 2)):
            lo, hi = nb_cert2.comp_bounds2(c, cr, dc, dw, hc, rm, DA, HA, RA, qw, cm, data[5], wgt, phase, um, hh)
            mid = 0.5 * (lo + hi); half = 0.5 * (hi - lo)
            for cf in cfs:
                for _ in range(3):
                    th = mp.mpf(um + random.uniform(-hh, hh))
                    s = mp.mpf(0)
                    for k in range(7):
                        a = (1 - mp.mpf(t)) * cf[0, 0, k] - mp.mpf(t) * cf[1, 0, k]
                        b = (1 - mp.mpf(t)) * cf[0, 1, k] - mp.mpf(t) * cf[1, 1, k]
                        wk = mp.mpf(float(wgt[k]))          # the weights the enclosure was computed for
                        s += wk * (a * mp.cos(k * th + phase) + b * mp.sin(k * th + phase))
                    nW += 1
                    if not (lo <= s <= hi):
                        vW += 1
                        if vW <= 3:
                            print('   VIOLATION (W)', bx, t, L, float(s), lo, hi, flush=True)
                    if half > 0:
                        wr = max(wr, float(abs(s - mp.mpf(mid)) / mp.mpf(half)))
print(f'(W) comp_bounds2 containment: {nW} checks, violations {vW}, worst |S - mid| / half-width {wr:.3f}', flush=True)

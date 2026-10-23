"""AUDIT-7 g21 (round-7 addendum 2): the second-order enclosure for the a-chart (gen_nb2a.py -> nbcoef_a2.py,
nb_cert2.box_coeffs2a / box_data2 / certify_box2 / certify_box_deep2; bb_run --order2).
  (A) nbcoef_a2.coef_val2 slots 0..4 == nbcoef_a.coef_val slots 0..4 (the audited first-order module);
  (B) slots 5..14 (second derivatives in (e, b, p, q), sorted pairs) == 50-digit numerical derivatives of slots 1..4;
  (C) coef_abs2 slots 0..14 term-wise |coef_val2| (textual); slots 15..34 == abs-coefficient polynomials of the sympy
      third derivatives of slot 0;
  (H) hcen of box_coeffs2a == 50-digit numerical Hessian in (eta, l, w) of C(sin eta, cos eta, l sin w, l cos w);
  (D) per-coefficient Taylor inequality |c(x) - cmid - dcen.d - 1/2 d'Hd| <= crad + rem (d in (eta, l, w));
  (E) enclose_point2 containment with the chart-0 data, modes 1, 2, 0 (as g19);
  (S) end-to-end: random boxes of the bulk domain (eta in [0.2003, 1.2536], l, w, t); every box certified by
      certify_box2(0, ...) is checked at 100 digits (own a-chart formula) at 6 points: N_mult <= total (<= 11) and
      even.  The first-order certify_box is run on the same boxes for comparison.
usage: python g21_nbcert2_achart.py NBOX NS SEED"""
import os, re, sys, types, random, math
import numpy as np
import mpmath as mp
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.normpath(os.path.join(HERE, '..', '..', 'code'))
sys.path.insert(0, CODE)
NBOX, NS, SEED = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
random.seed(SEED)
mp.mp.dps = 50
import nb_cert, nb_cert2
J = 6
PAIRS = [(0, 0), (0, 1), (0, 2), (0, 3), (1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3)]
TRIPLES = [(a, b, c) for a in range(4) for b in range(a, 4) for c in range(b, 4)]
ETA0, ETA1 = 0.20033740617957885, 1.2535563099909826


def load_plain(path):
    src = open(path).read().replace('from numba import njit', '')
    src = re.sub(r'@njit\([^)]*\)\n', '', src)
    mod = types.ModuleType(os.path.basename(path)); exec(compile(src, path, 'exec'), mod.__dict__)
    return mod, src


A1, _ = load_plain(os.path.join(CODE, 'nbcoef_a.py'))
A2, src2 = load_plain(os.path.join(CODE, 'nbcoef_a2.py'))


def ev(fn, nslot, x, obj=False):
    out = np.zeros((nslot, 2, 2, 7), dtype=object if obj else float)
    if obj:
        out[...] = mp.mpf(0)
    fn(*x, out)
    return out


def exact_c(eta, l, w):
    eta, l, w = mp.mpf(eta), mp.mpf(l), mp.mpf(w)
    return ev(A1.coef_val, 5, [mp.sin(eta), mp.cos(eta), l * mp.sin(w), l * mp.cos(w)], obj=True)[0]


# (A), (B)
wA = 0.0
for _ in range(40):
    x = [random.uniform(0, 1), random.uniform(0, 1), random.uniform(-1, 1), random.uniform(-1, 1)]
    o1 = ev(A1.coef_val, 5, x); o2 = ev(A2.coef_val2, 15, x)
    wA = max(wA, float(np.abs(o1 - o2[:5]).max() / (np.abs(o1).max() + 1e-300)))
wB = mp.mpf(0)
for _ in range(4):
    x0 = [mp.mpf(random.uniform(0, 1)) for _ in range(2)] + [mp.mpf(random.uniform(-1, 1)) for _ in range(2)]
    o2 = ev(A2.coef_val2, 15, x0, obj=True)
    h = mp.mpf(10) ** -20
    for m, (a, b) in enumerate(PAIRS):
        xp = list(x0); xm = list(x0); xp[b] += h; xm[b] -= h
        num = (ev(A2.coef_val2, 15, xp, obj=True)[1 + a] - ev(A2.coef_val2, 15, xm, obj=True)[1 + a]) / (2 * h)
        for f in range(2):
            for pt in range(2):
                for k in range(7):
                    wB = max(wB, abs(num[f, pt, k] - o2[5 + m, f, pt, k]) / (1 + abs(o2[5 + m, f, pt, k])))
# (C)
bv = dict(re.findall(r'out\[([^\]]+)\] = (.*)', src2.split('def coef_val2')[1].split('def coef_abs2')[0]))
ba = dict(re.findall(r'out\[([^\]]+)\] = (.*)', src2.split('def coef_abs2')[1]))
bad = nterm = 0
for key, exv in bv.items():
    tv = exv.split(' + '); ta = ba[key].split(' + ')
    if len(tv) != len(ta):
        bad += 1; continue
    for x_, y_ in zip(tv, ta):
        nterm += 1
        cx, mx = x_.split(' * ', 1) if ' * ' in x_ else (x_, '')
        cy, my = y_.split(' * ', 1) if ' * ' in y_ else (y_, '')
        if mx != my or abs(float(cx)) != float(cy):
            bad += 1
E, B, P, Q = sp.symbols('e b p q'); syms = (E, B, P, Q)
loc = {}
for n, s_ in zip('ebpq', syms):
    for d in range(1, 20):
        loc[f'{n}{d}'] = s_ ** d
wC = 0.0
pts = [[random.uniform(0, 1.3) for _ in range(4)] for _ in range(2)]
for f in range(2):
    for pt in range(2):
        for k in range(7):
            exs = bv[f'0, {f}, {pt}, {k}']
            ex = sp.Integer(0) if exs == '0.0' else sp.sympify(exs.replace('1.0 * ', ''), locals=loc)
            for m, (a, b, c) in enumerate(TRIPLES):
                d3 = sp.expand(sp.diff(ex, syms[a], syms[b], syms[c]))
                for p in pts:
                    got = ev(A2.coef_abs2, 35, p)[15 + m, f, pt, k]
                    want = 0.0 if d3 == 0 else float(sum(abs(cf) * sp.Mul(*[sp.Float(p[i]) ** mon[i] for i in range(4)])
                                                       for mon, cf in sp.Poly(d3, *syms).terms()))
                    wC = max(wC, abs(got - want) / (1 + abs(want)))
print(f'(A) slots 0..4 vs nbcoef_a: max rel diff {wA:.2e}; (B) second derivatives vs numerical: max rel err {mp.nstr(wB, 3)}; '
      f'(C) abs terms {nterm} (mismatches {bad}), third derivatives vs sympy max rel diff {wC:.2e}', flush=True)


def rand_box():
    wd = [10 ** random.uniform(-5, math.log10(0.2)) for _ in range(3)]
    e0 = random.uniform(0.0, math.pi / 2 - wd[0]); l0 = random.uniform(0.0, 1.0 - min(wd[1], 0.99)); w0 = random.uniform(0.0, max(0.0, math.pi / 2 - wd[2]))
    return (e0, e0 + wd[0], l0, min(l0 + wd[1], 1.0), w0, w0 + wd[2])


# (H), (D)
wH = mp.mpf(0); nD = vD = 0; rD = 0.0
for it in range(max(20, NBOX // 4)):
    bx = rand_box()
    cmid, crad, dcen, dwid, cmax, rr, hcen, rem = nb_cert2.box_coeffs2a(*bx)
    xc = [mp.mpf(0.5 * (bx[2 * i] + bx[2 * i + 1])) for i in range(3)]
    if it < 8:
        h = mp.mpf(10) ** -12
        for m, (a, b) in enumerate([(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]):
            def sh(sa, sb):
                x = list(xc); x[a] += sa * h; x[b] += sb * h
                return exact_c(*x)
            num = (sh(1, 1) - sh(1, -1) - sh(-1, 1) + sh(-1, -1)) / (4 * h * h)
            for f in range(2):
                for pt in range(2):
                    for k in range(7):
                        wH = max(wH, abs(num[f, pt, k] - hcen[m, f, pt, k]) / (1 + abs(num[f, pt, k])))
    for _ in range(4):
        x = [mp.mpf(random.uniform(bx[2 * i], bx[2 * i + 1])) for i in range(3)]
        d = [x[i] - xc[i] for i in range(3)]
        ce = exact_c(*x)
        for f in range(2):
            for pt in range(2):
                for k in range(7):
                    lin = sum(mp.mpf(float(dcen[v, f, pt, k])) * d[v] for v in range(3))
                    Hm = [mp.mpf(float(hcen[m, f, pt, k])) for m in range(6)]
                    quad = (Hm[0] * d[0] ** 2 + Hm[3] * d[1] ** 2 + Hm[5] * d[2] ** 2) / 2 + Hm[1] * d[0] * d[1] + \
                        Hm[2] * d[0] * d[2] + Hm[4] * d[1] * d[2]
                    r = abs(ce[f, pt, k] - mp.mpf(float(cmid[f, pt, k])) - lin - quad)
                    allowed = crad[f, pt, k] + rem[f, pt, k]
                    nD += 1
                    if r > allowed:
                        vD += 1
                        if vD <= 3:
                            print('   VIOLATION (D)', bx, f, pt, k, float(r), allowed, flush=True)
                    if allowed > 0:
                        rD = max(rD, float(r / allowed))
print(f'(H) Hessian (eta, l, w) vs numerical: max rel err {mp.nstr(wH, 3)}; (D) Taylor inequality: {nD} checks, '
      f'violations {vD}, worst ratio {rD:.4f}', flush=True)


# (E) containment
def tcomb(data, t):
    cmid, crad, dcen, dwid, cmax, rr, hcen, rem = data
    c = (1 - t) * cmid[0] - t * cmid[1]
    cr = (1 - t) * crad[0] + t * crad[1] + 1e-15 * (np.abs(cmid[0]) + np.abs(cmid[1]))
    dc = (1 - t) * dcen[:, 0] - t * dcen[:, 1]
    dw = (1 - t) * dwid[:, 0] + t * dwid[:, 1] + 1e-15 * (np.abs(dcen[:, 0]) + np.abs(dcen[:, 1]))
    hc = (1 - t) * hcen[:, 0] - t * hcen[:, 1]
    rm = (1 - t) * rem[0] + t * rem[1]
    qw = np.array([0.5 * rr[0] * rr[0], rr[0] * rr[1], rr[0] * rr[2], 0.5 * rr[1] * rr[1], rr[1] * rr[2], 0.5 * rr[2] * rr[2]])
    DA = np.abs(dcen[:, 0]) + np.abs(dcen[:, 1]); HA = np.abs(hcen[:, 0]) + np.abs(hcen[:, 1]); RA = rem[0] + rem[1]
    return c, cr, dc, dw, hc, rm, qw, DA, HA, RA


def exact_P(cf, t, u, j):
    t = mp.mpf(t); u = mp.mpf(u); s = mp.mpf(0)
    for k in range(7):
        if j > 0 and k == 0:
            continue
        a = (1 - t) * cf[0, 0, k] - t * cf[1, 0, k]; b = (1 - t) * cf[0, 1, k] - t * cf[1, 1, k]
        s += mp.mpf(k) ** j * (a * mp.cos(k * u + j * mp.pi / 2) + b * mp.sin(k * u + j * mp.pi / 2))
    return s


st = {0: [0, 0, 0.0], 1: [0, 0, 0.0], 2: [0, 0, 0.0]}
for _ in range(NBOX):
    bx = rand_box()
    t0 = random.uniform(0.0, 1.0); t1 = min(1.0, t0 + 10 ** random.uniform(-6, -0.5))
    data = nb_cert2.box_coeffs2a(*bx)
    pts = [[bx[2 * i + ((m >> i) & 1)] for i in range(3)] for m in random.sample(range(8), 3)]
    pts += [[random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(3)] for _ in range(3)]
    cfs = [exact_c(*x) for x in pts]
    for t in (t0, t1):
        c, cr, dc, dw, hc, rm, qw, DA, HA, RA = tcomb(data, t)
        u = random.uniform(-3 * math.pi, 3 * math.pi)
        for mode in (0, 1, 2):
            vc, rad = nb_cert2.enclose_point2(c, cr, dc, dw, hc, rm, DA, HA, RA, qw, data[5], u, J, mode)
            for j in range(J + 2):
                for cf in cfs:
                    dd = abs(exact_P(cf, t, u, j) - mp.mpf(float(vc[j])))
                    s_ = st[mode]; s_[0] += 1
                    if dd > rad[j]:
                        s_[1] += 1
                        if s_[1] <= 3:
                            print('   VIOLATION (E) mode', mode, bx, t, u, j, float(dd), rad[j], flush=True)
                    if rad[j] > 0:
                        s_[2] = max(s_[2], float(dd / rad[j]))
for mode, nm in ((1, 'second-order only'), (2, 'mean-value only'), (0, 'min (used)')):
    s_ = st[mode]
    print(f'(E) enclose_point2 (chart 0) mode {mode} ({nm}): {s_[0]} checks, violations {s_[1]}, worst ratio {s_[2]:.7f}', flush=True)


# (S) end-to-end
def pol_mult(eta, l, w, t):
    e = mp.sin(mp.mpf(eta)); b = mp.cos(mp.mpf(eta)); p = mp.mpf(l) * mp.sin(mp.mpf(w)); q = mp.mpf(l) * mp.cos(mp.mpf(w)); t = mp.mpf(t)
    N = 32; vals = []
    for n in range(N):
        u = 2 * mp.pi * n / N
        r = 1 - e * mp.cos(u); Z = p * (mp.cos(u) - e) + q * b * mp.sin(u); Zp = -p * mp.sin(u) + q * b * mp.cos(u)
        A = 2 * e * mp.sin(u) * r; rho2 = r * r - Z * Z; d = A - 2 * Z * Zp
        vals.append((1 - t) * A * A * rho2 - t * d * d)
    c = [mp.fsum(vals[n] * mp.expj(-k * 2 * mp.pi * n / N) for n in range(N)) / N for k in range(-6, 7)]
    sc = max(abs(v) for v in c); lo = 0
    while abs(c[lo]) < mp.mpf(10) ** -40 * sc:
        lo += 1
    hi = 12
    while abs(c[hi]) < mp.mpf(10) ** -40 * sc:
        hi -= 1
    z = mp.polyroots(c[lo:hi + 1][::-1], maxsteps=600, extraprec=400) if hi > lo else []
    return sum(1 for r in z if abs(mp.log(abs(r))) < mp.mpf(10) ** -25)


nS = cert2 = cert1 = only2 = vS = odd = 0
for _ in range(NS):
    wd = [10 ** random.uniform(-4, -1.3) for _ in range(4)]
    x0 = random.uniform(ETA0, ETA1 - wd[0]); l0 = random.uniform(0, 1 - wd[1]); w0 = random.uniform(0, math.pi / 2 - wd[2])
    t0 = random.uniform(0, 1 - wd[3])
    bx = (x0, x0 + wd[0], l0, l0 + wd[1], w0, w0 + wd[2], t0, t0 + wd[3])
    ok2, T2 = nb_cert2.certify_box2(0, *bx, 6, 64, 8)
    ok1, T1 = nb_cert.certify_box(0, *bx, 6, 64, 8)
    nS += 1; cert2 += ok2; cert1 += ok1; only2 += ok2 and not ok1
    if not ok2:
        continue
    for k in range(6):
        x = [0.5 * (bx[2 * i] + bx[2 * i + 1]) if k == 0 else random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(4)]
        n = pol_mult(*x)
        if n % 2:
            odd += 1
        if n > T2:
            vS += 1; print('   VIOLATION (S)', bx, x, n, T2, flush=True)
print(f'(S) end-to-end: {nS} random bulk boxes; certified by certify_box2 {cert2} (first-order certify_box {cert1}; '
      f'second-order only {only2}); 100-digit violations {vS}; odd N_mult {odd}', flush=True)

"""AUDIT-4 d6: soundness of wind_cert.comp_bounds / labels against an independent evaluator.
Charts 'a' (nb_cert.box_coeffs chart 0), 'B1' and 'C1pp' (generated engines).  For random boxes, t-intervals and test
circles L in [-6, -0.05]: for each arc of a 64-arc cover, the intervals [lo, hi] of Re G and Im G at each t-end must
contain the true values at sampled (x in box, theta in padded arc); and every label certified by labels() over
[t0, t1] must hold at sampled (x, t, theta).  True G from MY chart maps + a-chart physics (double precision DFT).
usage: python d6_comp_bounds.py NBOX SEED  (.venv python)"""
import os, sys, math, random, importlib, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.join(HERE, '..', '..', 'code'); sys.path.insert(0, CODE)
import wind_cert, nb_cert
NB = int(sys.argv[1]) if len(sys.argv) > 1 else 40
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 6)


def coeffs_phys(e, p, q, t):
    N = 32; u = 2 * np.pi * np.arange(N) / N; b = math.sqrt(1 - e * e)
    r = 1 - e * np.cos(u); Z = p * (np.cos(u) - e) + q * b * np.sin(u); Zp = -p * np.sin(u) + q * b * np.cos(u)
    A = 2 * e * np.sin(u) * r; rho2 = r * r - Z * Z; d = A - 2 * Z * Zp
    F = (1 - t) * A * A * rho2 - t * d * d
    cf = np.fft.rfft(F) / N
    a = np.zeros(7); s = np.zeros(7); a[0] = cf[0].real; a[1:] = 2 * cf[1:7].real; s[1:] = -2 * cf[1:7].imag
    return a, s


def chart(eng, x, x3):
    """returns (a, s) Fourier coefficients of the chart function up to a positive factor at box point x, x3=t-variable"""
    if eng == 'a':
        eta, l, w = x; e = math.sin(eta); return coeffs_phys(e, l * math.sin(w), l * math.cos(w), x3), 1.0
    if eng == 'B1':
        e, ph, qh = x; tauh = -2 + 4 * x3
        a, s = coeffs_phys(e, e * ph, e * qh, 0.5 + e * tauh); return (a / e ** 3, s / e ** 3), 1.0
    if eng == 'C1pp':
        epp, l, w = x; e = epp * l * l; tt = epp ** 2 * x3 / (1 - x3)
        a, s = coeffs_phys(e, l * math.sin(w), l * math.cos(w), tt); return (a * (1 - x3) / e ** 2, s * (1 - x3) / e ** 2), 1.0


def G(a, s, L, th):
    k = np.arange(7)
    re = np.sum(np.cosh(k * L) * (a * np.cos(k * th) + s * np.sin(k * th)))
    im = np.sum(np.sinh(k * L) * (a * np.sin(k * th) - s * np.cos(k * th)))
    return re, im, np.sum(np.cosh(k * L) * (np.abs(a) + np.abs(s)))


DOMS = {'a': [(math.asin(0.2), math.asin(0.95)), (0, 1), (0, math.pi / 2), (0, 1)],
        'B1': [(0.005, 0.07), (0, 1.75), (0, 1.75), (0, 1)],
        'C1pp': [(0.01, 0.51), (0.05, 0.52), (0, math.pi / 2), (0, 0.98)]}
st = {}
for eng, dom in DOMS.items():
    if eng == 'a':
        bcf = lambda bx: nb_cert.box_coeffs(0, *bx[:6])
    elif eng == 'C1pp':
        m = importlib.import_module('nbeng_C1pp'); bcf = lambda bx, m=m: m.box_coeffs(*bx[:6])
    else:
        m = importlib.import_module('nbeng_B1'); bcf = lambda bx, m=m: m.box_coeffs(*bx)
    nchk = 0; viol = 0; lab_chk = 0; lab_viol = 0; nlab = 0
    for ib in range(NB):
        c = [rng.uniform(*d) for d in dom]
        r = np.array([d[1] - d[0] for d in dom]) * 0.01 * np.exp(rng.uniform(-3, 0, 4))
        bx = []
        for i in range(4):
            bx += [max(c[i] - r[i], dom[i][0]), min(c[i] + r[i], dom[i][1])]
        cm = bcf(tuple(bx)); cmid, crad, dcen, dwid, cmax, rr = cm
        L = -10 ** rng.uniform(-1.3, math.log10(6))
        M = 64; left = np.arange(M) * (2 * math.pi / M); width = np.full(M, 2 * math.pi / M)
        um = left + width / 2; hh = width / 2 + wind_cert.ARC_PAD
        wR = np.array([math.cosh(k * L) * (1 + 1e-12) for k in range(7)]); wI = np.array([math.sinh(k * L) * (1 + 1e-12) for k in range(7)])
        pts = [tuple(bx[2 * i + b_] for i, b_ in enumerate(bits)) for bits in itertools.product((0, 1), repeat=3)]
        pts += [tuple(rng.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(3)) for _ in range(6)]
        for t in (bx[6], bx[7]):
            cc = (1 - t) * cmid[0] - t * cmid[1]
            cr = (1 - t) * crad[0] + t * crad[1] + 1e-15 * (np.abs(cmid[0]) + np.abs(cmid[1]))
            dc = (1 - t) * dcen[:, 0] - t * dcen[:, 1]
            dw = (1 - t) * dwid[:, 0] + t * dwid[:, 1] + 1e-15 * (np.abs(dcen[:, 0]) + np.abs(dcen[:, 1]))
            cmx = (1 - t) * cmax[0] + t * cmax[1]
            bounds = [(wind_cert.comp_bounds(cc, cr, dc, dw, cmx, rr, wR, 0.0, um[i], hh[i]),
                       wind_cert.comp_bounds(cc, cr, dc, dw, cmx, rr, wI, -math.pi / 2, um[i], hh[i])) for i in range(M)]
            for x in pts:
                (a, s), _ = chart(eng, x, t)
                for i in range(M):
                    for th in (um[i] - hh[i], um[i] + hh[i], um[i] + rng.uniform(-1, 1) * hh[i]):
                        re, im, sc = G(a, s, L, th)
                        noise = 1e-12 * sc
                        (lR, hR), (lI, hI) = bounds[i]
                        nchk += 2
                        if re < lR - noise or re > hR + noise or im < lI - noise or im > hI + noise:
                            viol += 1
        lab = wind_cert.labels(cmid, crad, dcen, dwid, cmax, rr, bx[6], bx[7], L, um, hh)
        nlab += int(np.sum(lab > 0))
        for x in pts[::2]:
            t = rng.uniform(bx[6], bx[7])
            (a, s), _ = chart(eng, x, t)
            for i in range(M):
                if lab[i] == 0:
                    continue
                th = um[i] + rng.uniform(-1, 1) * hh[i]
                re, im, sc = G(a, s, L, th)
                comp = {1: re, 2: im, 3: -re, 4: -im}[int(lab[i])]
                lab_chk += 1
                if comp < -1e-12 * sc:
                    lab_viol += 1
    st[eng] = dict(enclosure_checks=nchk, enclosure_violations=viol, labelled_arcs=nlab, label_checks=lab_chk, label_violations=lab_viol)
    print(eng, st[eng], flush=True)

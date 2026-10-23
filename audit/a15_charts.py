"""AUDIT a15: near-circular blow-up charts C1, C1p, C1pp (gen_chart.py -> nbcoef_<chart>.py).

1. Values: coef_val(atoms) out[0] vs an independent evaluation built from the a-chart closed forms (r, Z, rho^2, A)
   at the mapped (e, l, w), scaled exactly as the chart claims:
      C1  : F_A = Phi/eps^2,          F_B = Psi/eps^2              (e = eps cos th, l^2 = eps sin th)
      C1p : F_A = Phi/e^2,            F_B = eps^2 Phi/e^2 + Psi    (e = eps)
      C1pp: F_A = Phi/e^2,            F_B = epp^2 Phi/e^2 + Psi/l^4 (e = epp l^2)
   Fourier coefficients by projection sums (linear, hence complex-analytic in the parameters).
2. Gradients out[1+v] (d/d box variable, chain rule through DER) vs complex-step derivatives of the evaluator.
3. Absolute-value polynomials: coef_abs == sum |c_a| atoms^a of the coef_val polynomial (sympy, exact).
4. Identity check of the t-reparametrisation: (1-tau) f_t / e^2 == (1-tau) F_A - tau F_B, t = eps^2 tau/(1-tau)
   (C1p) and t = epp^2 tau/(1-tau) (C1pp), at random points.
5. Coverage of {0 < e < 0.2} by the three chart domains, sampled, and the seam with the bulk a-chart (eta >= 0.20136).
"""
import os
import re
import sys
import types
import math
import numpy as np
import sympy as sp

CODE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'code')
rng = np.random.default_rng(1515)


def load_plain(path):
    src = open(path).read()
    src = src.replace('from numba import njit', '')
    src = re.sub(r'@njit\([^)]*\)\n', '', src)
    mod = types.ModuleType(os.path.basename(path))
    exec(compile(src, path, 'exec'), mod.__dict__)
    return mod


def phipsi(e, l, w, N=32):
    """complex-safe Phi, Psi samples from the a-chart closed forms"""
    u = 2 * np.pi * np.arange(N) / N
    b = np.sqrt(1 - e * e + 0j) if np.iscomplexobj(e) else np.sqrt(1 - e * e)
    p, q = l * np.sin(w), l * np.cos(w)
    r = 1 - e * np.cos(u)
    Z = p * (np.cos(u) - e) + b * q * np.sin(u)
    Zp = -p * np.sin(u) + b * q * np.cos(u)
    A = 2 * e * np.sin(u) * r
    rho2 = r * r - Z * Z
    d = A - 2 * Z * Zp
    return A * A * rho2, d * d, u


def proj(F, u):
    N = len(u)
    a = np.zeros(7, dtype=F.dtype); s = np.zeros(7, dtype=F.dtype)
    a[0] = F.sum() / N
    for k in range(1, 7):
        a[k] = 2 * (F * np.cos(k * u)).sum() / N
        s[k] = 2 * (F * np.sin(k * u)).sum() / N
    return a, s


def chart_funcs(chart, x):
    x0, x1, x2 = x
    if chart == 'C1':
        eps, th, w = x0, x1, x2
        e = eps * np.cos(th); l = np.sqrt(eps * np.sin(th))
        Phi, Psi, u = phipsi(e, l, w)
        FA, FB = Phi / eps ** 2, Psi / eps ** 2
    elif chart == 'C1p':
        eps, l, w = x0, x1, x2
        e = eps
        Phi, Psi, u = phipsi(e, l, w)
        FA, FB = Phi / e ** 2, eps ** 2 * Phi / e ** 2 + Psi
    else:
        epp, l, w = x0, x1, x2
        e = epp * l ** 2
        Phi, Psi, u = phipsi(e, l, w)
        FA, FB = Phi / e ** 2, epp ** 2 * Phi / e ** 2 + Psi / l ** 4
    return proj(FA, u), proj(FB, u)


def atoms_of(chart, x):
    x0, x1, x2 = x
    if chart == 'C1':
        e = x0 * math.cos(x1); b = math.sqrt(1 - e * e)
        return np.array([x0, math.cos(x1), math.sin(x1), b, 1 / b, math.sin(x2), math.cos(x2)])
    if chart == 'C1p':
        b = math.sqrt(1 - x0 * x0)
        return np.array([x0, x1, b, 1 / b, math.sin(x2), math.cos(x2)])
    e = x0 * x1 * x1; b = math.sqrt(1 - e * e)
    return np.array([x0, x1, b, 1 / b, math.sin(x2), math.cos(x2)])


DOM = {'C1': ((1e-3, 0.26), (0.0, math.atan(2.0)), (0, math.pi / 2)),
       'C1p': ((1e-3, 0.2), (0.4, 1.0), (0, math.pi / 2)),
       'C1pp': ((1e-3, 0.5), (0.05, 0.52), (0, math.pi / 2))}

for chart in ('C1', 'C1p', 'C1pp'):
    mod = load_plain(os.path.join(CODE, f'nbcoef_{chart}.py'))
    errv = 0.0; errd = 0.0; magv = 0.0
    for it in range(60):
        x = np.array([rng.uniform(*DOM[chart][i]) for i in range(3)])
        a = atoms_of(chart, x)
        out = np.zeros((4, 2, 2, 7))
        mod.coef_val(a, out)
        (aA, sA), (aB, sB) = chart_funcs(chart, x)
        want = np.array([[aA, sA], [aB, sB]]).real
        errv = max(errv, float(np.max(np.abs(out[0] - want))))
        magv = max(magv, float(np.max(np.abs(want))))
        for v in range(3):
            h = 1e-20
            xc = x.astype(complex); xc[v] += 1j * h
            (aA, sA), (aB, sB) = chart_funcs(chart, xc)
            dwant = np.array([[aA, sA], [aB, sB]]).imag / h
            errd = max(errd, float(np.max(np.abs(out[1 + v] - dwant)) / (1 + np.max(np.abs(dwant)))))
    print(f"{chart}: values max abs err {errv:.2e} (max |coef| {magv:.2e}); gradients max rel err {errd:.2e}")
    # 3. absolute-value polynomials (exact)
    names = mod.ATOMS
    syms = sp.symbols(' '.join(names), real=True)
    ys = sp.symbols(' '.join('y_' + n for n in names), positive=True)
    ov = np.empty((4, 2, 2, 7), dtype=object); oa = np.empty((4, 2, 2, 7), dtype=object)
    mod.coef_val(list(syms), ov); mod.coef_abs(list(ys), oa)
    worst = 0.0
    for idx in np.ndindex(4, 2, 2, 7):
        ex = sp.expand(sp.sympify(ov[idx]))
        if ex == 0:
            want = sp.Integer(0)
        else:
            P = sp.Poly(ex, *syms)
            want = sum(abs(c) * sp.Mul(*[y ** d for y, d in zip(ys, m)]) for m, c in P.terms())
        diff = sp.expand(sp.sympify(oa[idx]) - want)
        if diff != 0:
            worst = max(worst, max(abs(float(c)) for c in sp.Poly(diff, *ys).coeffs()))
    print(f"{chart}: coef_abs vs sum|c| monomials: max coefficient difference {worst:.2e}")

# 4. t-reparametrisation identity
worst = 0.0
for it in range(200):
    e = rng.uniform(1e-3, 0.2); l = rng.uniform(0.05, 1); w = rng.uniform(0, np.pi / 2); tau = rng.uniform(0, 0.99)
    for chart in ('C1p', 'C1pp'):
        if chart == 'C1p':
            x = (e, l, w); scale = e
        else:
            x = (e / l ** 2, l, w); scale = e / l ** 2
        t = scale ** 2 * tau / (1 - tau)
        Phi, Psi, u = phipsi(e, l, w)
        lhs = (1 - tau) * ((1 - t) * Phi - t * Psi) / e ** 2
        (aA, sA), (aB, sB) = chart_funcs(chart, x)
        FA = aA[0] + sum(aA[k] * np.cos(k * u) + sA[k] * np.sin(k * u) for k in range(1, 7))
        FB = aB[0] + sum(aB[k] * np.cos(k * u) + sB[k] * np.sin(k * u) for k in range(1, 7))
        rhs = (1 - tau) * FA - tau * FB
        worst = max(worst, float(np.max(np.abs(lhs - rhs)) / (1 + np.max(np.abs(lhs)))))
print(f"t-reparametrisation identity (C1p, C1pp): max rel err {worst:.2e}")

# 5. coverage
TH = math.atan(2.0)          # arccot 0.5
miss = []; n = 0
for it in range(400000):
    e = (rng.uniform(0, 0.2000020372) if it % 2 else 10 ** rng.uniform(-8, math.log10(0.2000020372))) if it % 3 else rng.uniform(0.199999, 0.2000020372)
    l = rng.uniform(0, 1) if it % 3 else 10 ** rng.uniform(-6, 0)
    n += 1
    inC1p = (0.0 <= e <= 0.2) and (0.4 <= l <= 1.0)
    eps = math.hypot(e, l * l); th = math.atan2(l * l, e)
    inC1 = (eps <= 0.26) and (th <= TH)
    inC1pp = (l > 0) and (e / l ** 2 <= 0.5) and (l <= 0.52)
    inBulk = e >= math.sin(0.20136)
    if not (inC1p or inC1 or inC1pp or inBulk):
        miss.append((e, l))
print(f"coverage sample of 0 < e <= sin(0.20136) = {math.sin(0.20136):.9f}, l in [0,1]: {len(miss)} uncovered of {n}")
if miss:
    m = np.array(miss)
    print("   uncovered e range [%.9f, %.9f], l range [%.6f, %.6f]" % (m[:, 0].min(), m[:, 0].max(), m[:, 1].min(), m[:, 1].max()))

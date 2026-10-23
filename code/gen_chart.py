"""Generic chart generator: a chart is f_x3 = (1 - x3) F_A - x3 F_B, with F_A, F_B trigonometric polynomials of degree
<= 6 in u whose coefficients are polynomials in 'atoms'.  Atoms are functions of the box variables (x0, x1, x2) whose
derivatives are again polynomials in atoms.

Generates nbcoef_<chart>.py with
  coef_val(a, out): out[0, F, part, k] = coefficient; out[1 + v, F, part, k] = d/dx_v of it (chain rule);
  coef_abs(a, out): the same polynomials with all monomial coefficients replaced by their absolute values.
F in {0: F_A, 1: F_B}; part in {0: cos k u, 1: sin k u}; k = 0..6.

Charts near e = 0 (all built from the a-chart: r = 1 - e cos u, Z = l zt, zt = sw (cos u - e) + b cw sin u,
rho^2 = r^2 - l^2 zt^2, A = (r^2)' = 2 e sin u r, (rho^2)' = A - l^2 (zt^2)', Phi = A^2 rho^2, Psi = ((rho^2)')^2):
  C1   : box (eps, th, w), x3 = t.        e = eps cos th, l^2 = eps sin th.     F_A = Phi/eps^2, F_B = Psi/eps^2.
  C1p  : box (eps, l, w),  x3 = tau.      e = eps, t = eps^2 tau/(1-tau).        F_A = Phi/e^2,  F_B = eps^2 Phi/e^2 + Psi.
  C1pp : box (epp, l, w),  x3 = tau.      e = epp l^2, t = epp^2 tau/(1-tau).   F_A = Phi/e^2,  F_B = epp^2 Phi/e^2 + Psi/l^4.
(For C1p/C1pp: (1 - tau) * f/e^2 = (1 - tau) F_A - tau F_B with the stated F_A, F_B.)
"""
import sys
import sympy as sp

chart = sys.argv[1]
z = sp.symbols('z')
cu = (z + 1 / z) / 2
su = (z - 1 / z) / (2 * sp.I)


def ddu(expr):
    return sp.expand(sp.I * z * sp.diff(expr, z))


def build(e, b, l2, sw, cw):
    """returns (A/?, pieces) in terms of given e, b, l^2 expressions"""
    r = 1 - e * cu
    zt = sw * (cu - e) + b * cw * su
    rho2 = sp.expand(r**2 - l2 * zt**2)
    return r, zt, rho2


if chart == 'C1':
    eps, cth, sth, b, ib, sw, cw = sp.symbols('eps cth sth b ib sw cw', real=True)
    ATOMS = [eps, cth, sth, b, ib, sw, cw]
    e = eps * cth; l2 = eps * sth
    r, zt, rho2 = build(e, b, l2, sw, cw)
    A_over = sp.expand(2 * cth * su * r)                                   # A / eps
    dzt2 = ddu(sp.expand(zt**2))
    FA = sp.expand(A_over**2 * rho2)                                      # Phi / eps^2
    FB = sp.expand((A_over - sth * dzt2)**2)                              # Psi / eps^2
    # d(atom)/d(box var) as polynomials in atoms; box vars x0 = eps, x1 = th, x2 = w
    DER = {eps: [1, 0, 0], cth: [0, -sth, 0], sth: [0, cth, 0],
           b: [-eps * cth**2 * ib, eps**2 * cth * sth * ib, 0], ib: [eps * cth**2 * ib**3, -eps**2 * cth * sth * ib**3, 0],
           sw: [0, 0, cw], cw: [0, 0, -sw]}
elif chart == 'C1p':
    eps, l, b, ib, sw, cw = sp.symbols('eps l b ib sw cw', real=True)
    ATOMS = [eps, l, b, ib, sw, cw]
    e = eps; l2 = l**2
    r, zt, rho2 = build(e, b, l2, sw, cw)
    A_over = sp.expand(2 * su * r)                                        # A / e
    dzt2 = ddu(sp.expand(zt**2))
    PhiE = sp.expand(A_over**2 * rho2)                                    # Phi / e^2
    Psi = sp.expand((eps * A_over - l2 * dzt2)**2)
    FA = PhiE
    FB = sp.expand(eps**2 * PhiE + Psi)
    DER = {eps: [1, 0, 0], l: [0, 1, 0], b: [-eps * ib, 0, 0], ib: [eps * ib**3, 0, 0],
           sw: [0, 0, cw], cw: [0, 0, -sw]}
elif chart == 'C1pp':
    epp, l, b, ib, sw, cw = sp.symbols('epp l b ib sw cw', real=True)
    ATOMS = [epp, l, b, ib, sw, cw]
    e = epp * l**2; l2 = l**2
    r, zt, rho2 = build(e, b, l2, sw, cw)
    A_over = sp.expand(2 * su * r)                                        # A / e
    dzt2 = ddu(sp.expand(zt**2))
    PhiE = sp.expand(A_over**2 * rho2)                                    # Phi / e^2
    PsiL = sp.expand((epp * A_over - dzt2)**2)                            # Psi / l^4
    FA = PhiE
    FB = sp.expand(epp**2 * PhiE + PsiL)
    # e = epp l^2, b = sqrt(1-e^2): db/depp = -e l^2 / b = -epp l^4 ib ; db/dl = -e * 2 epp l / b = -2 epp^2 l^3 ib
    DER = {epp: [1, 0, 0], l: [0, 1, 0], b: [-epp * l**4 * ib, -2 * epp**2 * l**3 * ib, 0],
           ib: [epp * l**4 * ib**3, 2 * epp**2 * l**3 * ib**3, 0], sw: [0, 0, cw], cw: [0, 0, -sw]}
else:
    raise SystemExit("unknown chart")


def real_form(F):
    P = sp.Poly(sp.expand(F * z**6), z)
    assert P.degree() <= 12, P.degree()
    for k in range(13):
        pass
    c = {k - 6: sp.expand(P.coeff_monomial(z**k)) for k in range(13)}
    a = [sp.expand(sp.re(c[0]))] + [sp.expand(2 * sp.re(c[k])) for k in range(1, 7)]
    s = [sp.Integer(0)] + [sp.expand(-2 * sp.im(c[k])) for k in range(1, 7)]
    assert sp.expand(sp.im(c[0])) == 0
    return a, s


polys = {}
for fi, F in enumerate((FA, FB)):
    a, s = real_form(F)
    for part, arr in ((0, a), (1, s)):
        for k in range(7):
            ex = arr[k]
            polys[(0, fi, part, k)] = ex
            for v in range(3):
                d = sum(sp.diff(ex, at) * DER[at][v] for at in ATOMS)
                polys[(1 + v, fi, part, k)] = sp.expand(d)

names = [str(a) for a in ATOMS]
maxdeg = {n: 0 for n in names}
for ex in polys.values():
    if ex == 0:
        continue
    for mon in sp.Poly(ex, *ATOMS).monoms():
        for n, d in zip(names, mon):
            maxdeg[n] = max(maxdeg[n], d)


def emit(fn, use_abs):
    L = ['@njit(cache=True)', f'def {fn}(a, out):']
    for i, n in enumerate(names):
        L.append(f'    {n}0 = 1.0')
        if maxdeg[n] >= 1:
            L.append(f'    {n}1 = a[{i}]')
        for d in range(2, maxdeg[n] + 1):
            L.append(f'    {n}{d} = {n}{d-1} * {n}1')
    for (dv, fi, part, k), ex in polys.items():
        idx = f'{dv}, {fi}, {part}, {k}'
        if ex == 0:
            L.append(f'    out[{idx}] = 0.0')
            continue
        terms = []
        for mon, coef in sp.Poly(ex, *ATOMS).terms():
            cval = float(abs(coef)) if use_abs else float(coef)
            mono = ' * '.join(f'{n}{d}' for n, d in zip(names, mon) if d > 0) or '1.0'
            terms.append(f'{cval!r} * {mono}')
        L.append(f'    out[{idx}] = ' + ' + '.join(terms))
    L.append('')
    return L


src = [f'"""Generated by gen_chart.py, chart {chart}. atoms = {names}."""', 'from numba import njit', '',
       f'ATOMS = {names!r}', '']
src += emit('coef_val', False)
src += emit('coef_abs', True)
open(f'nbcoef_{chart}.py', 'w').write('\n'.join(src) + '\n')
print(f'nbcoef_{chart}.py written; atoms {names}; max degrees {maxdeg}')

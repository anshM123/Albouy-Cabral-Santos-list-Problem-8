"""C3 core charts (e -> 1, mu -> 0, aphelion direction -> circle axis).  See C3_CORE.md.

R-normalised corner function (circle radius 1, pericentre distance mu, k = 2e(1+e)):
  P(f) = mu^2 k^2 sin^2 f sigma^2 - D^2 B^2,  zeta = Pz cos f + Qz sin f, sigma^2 = 1 - zeta^2,
  B = (sigma^2)' + 2e (sin f - Qz zeta),  D = 1 + e cos f.
Anomaly variable th = tan((f - pi)/2):  Q(th) = P (1 + th^2)^6, a polynomial of degree <= 12 (f = 0 <-> th = oo).
Parameters near C3: Qz = l cos w >= 0, kappa = 1 - Pz >= 0, eps = 1 - e >= 0, mu >= 0; weights (1, 2, 2, 3).
Weighted radius delta = max(Qz, sqrt(kappa), sqrt(eps), mu^(1/3)).

Level-1 parameter cone charts (nu = delta in [0, 0.37]; ratios in [0, R]); root variable y = th / nu in [-Y0, Y0]:
  Ca: Qz = nu,        kappa = nu^2 K, eps = nu^2 E, mu = nu^3 m      box (nu, K, E, m)
  CK: kappa = nu^2,   Qz = nu a,      eps = nu^2 E, mu = nu^3 m      box (nu, a, E, m)
  CE: eps = nu^2,     Qz = nu a,      kappa = nu^2 K, mu = nu^3 m    box (nu, a, K, m)
  Cm: mu = nu^3,      Qz = nu a,      kappa = nu^2 K, eps = nu^2 E   box (nu, a, K, E)
  G = Q(nu y) / nu^10  (exact; only even powers of nu occur).
Intermediate cone charts (th dominant): th = +-Y0 nu, (Qz, kappa, eps, mu) = (nu a, nu^2 K, nu^2 E, nu^3 m), ratios in
  [0, 1]; G = Q(+-Y0 nu)/nu^10; certificate: no zero (root variable m in [0, 1]); box (nu, a, K, E).
Level 2 (nested point y = E = m = 0 inside Ca and CK): weights (1, 2, 2) for (y, E, m), scales (Y2, E2, m2):
  *E: E = E2 eta^2, m = m2 eta^2 mp, y = Y2 eta Y     box (nu, K|a, eta, mp), root Y in [-R, R]
  *m: m = m2 eta^2, E = E2 eta^2 Ep, y = Y2 eta Y     box (nu, K|a, eta, Ep), root Y in [-R, R]
  *Yp/*Ym: y = +-Y2 eta, E = E2 eta^2 Ep, m = m2 eta^2 mp   box (nu, K|a, eta, Ep), root variable mp in [0, R]: no zero
  G2 = G / eta^6 (exact).
Writes nbcoef_c3_<chart>.py: coef_val(x, out) with out[0, k] = coefficient of root^k, out[1 + v, k] = d/dx_v;
coef_abs(x, out) the absolute-value polynomials.  Self-test against the direct trigonometric formula.
usage: python gen_c3.py CHART [CHART ...]"""
import sys
import math
import random
import sympy as sp

Y0 = 3.0
YF = 3.0        # face charts (e -> 1 and mu -> 0 at general Pz): root range |y| <= YF
Y2, E2, M2 = 0.25, 0.01, 0.01
th, Pz, Qz, e, mu = sp.symbols('th Pz Qz e mu', real=True)


def corner_Q(centre='aphelion'):
    """Q(th) = P(f)(1 + th^2)^6 with f = f_c + 2 atan(th): f_c = pi (aphelion, default), 3 pi/2 ('right') or pi/2
    ('left')"""
    cf, sf = sp.symbols('cf sf', real=True)
    zeta = Pz * cf + Qz * sf
    dzeta = -Pz * sf + Qz * cf
    B = -2 * zeta * dzeta + 2 * e * (sf - Qz * zeta)
    D = 1 + e * cf
    k = 2 * e * (1 + e)
    P = sp.expand(mu**2 * k**2 * sf**2 * (1 - zeta**2) - D**2 * B**2)
    den = 1 + th**2
    if centre == 'right':       # f = 3 pi/2 + x: cos f = sin x, sin f = -cos x
        sb = {cf: 2 * th / den, sf: -(1 - th**2) / den}
    elif centre == 'left':      # f = pi/2 + x: cos f = -sin x, sin f = cos x
        sb = {cf: -2 * th / den, sf: (1 - th**2) / den}
    else:                       # f = pi + x
        sb = {cf: -(1 - th**2) / den, sf: -2 * th / den}
    return sp.expand(sp.cancel(P.subs(sb) * den**6))


nu, y, K, E, m, a, eta, mp, Ep, Yv = sp.symbols('nu y K E m a eta mp Ep Y', real=True)
sF = sp.symbols('s', real=True)


def level1(chart):
    """returns (G as polynomial expr in y and box vars, box vars, physical map)"""
    Qx = corner_Q()
    if chart == 'Ca':
        sub = {Qz: nu, Pz: 1 - nu**2 * K, e: 1 - nu**2 * E, mu: nu**3 * m}; box = [nu, K, E, m]
    elif chart == 'CK':
        sub = {Qz: nu * a, Pz: 1 - nu**2, e: 1 - nu**2 * E, mu: nu**3 * m}; box = [nu, a, E, m]
    elif chart == 'CE':
        sub = {Qz: nu * a, Pz: 1 - nu**2 * K, e: 1 - nu**2, mu: nu**3 * m}; box = [nu, a, K, m]
    elif chart == 'Cm':
        sub = {Qz: nu * a, Pz: 1 - nu**2 * K, e: 1 - nu**2 * E, mu: nu**3}; box = [nu, a, K, E]
    else:
        raise ValueError(chart)
    sub[th] = nu * y
    W = sp.Poly(sp.expand(Qx.subs(sub)), nu)
    low = min(mon[0] for mon in W.monoms())
    assert low == 10, (chart, low)
    G = sp.expand(W.as_expr() / nu**10)
    return G, box


def face_G(sub, root):
    """face blow-up (weights th:1, eps:2, mu:2): G = Q / s^6, exact"""
    W = sp.Poly(sp.expand(corner_Q().subs(sub)), sF)
    assert min(mon[0] for mon in W.monoms()) == 6
    return sp.expand(W.as_expr() / sF**6)


def chart_def(chart):
    # Face charts at e = 1, mu = 0 for general Pz (C3_CORE.md section 6): eps = 1 - e and mu have weight 2, th weight 1.
    #   FE: eps = s^2, mu = s^2 m       box (s, m, Pz, Qz), root y = th / s        (eps dominant)
    #   FM: mu = s^2, eps = s^2 E       box (s, E, Pz, Qz), root y = th / s        (mu dominant)
    #   FIP/FIM: th = +-YF s, eps = s^2 E, mu = s^2 m   box (s, E, Pz, Qz), root variable m: no zero (intermediate)
    if chart == 'FE':
        return dict(G=face_G({th: sF * y, e: 1 - sF**2, mu: sF**2 * m}, y), box=[sF, m, Pz, Qz], root=y)
    if chart == 'FM':
        return dict(G=face_G({th: sF * y, e: 1 - sF**2 * E, mu: sF**2}, y), box=[sF, E, Pz, Qz], root=y)
    if chart in ('FIP', 'FIM'):
        sg = 1 if chart == 'FIP' else -1
        return dict(G=face_G({th: sg * YF * sF, e: 1 - sF**2 * E, mu: sF**2 * m}, m), box=[sF, E, Pz, Qz], root=m)
    if chart in ('Ca', 'CK', 'CE', 'Cm'):
        G, box = level1(chart)
        return dict(G=G, box=box, root=y)
    if chart == 'KBQ':
        # band chart: physical variables, no rescaling (used on 0.1 <= |f - pi| <= 0.9)
        return dict(G=corner_Q(), box=[e, Pz, Qz, mu], root=th)
    if chart == 'KBQe':
        # band chart in eps = 1 - e (no cancellation at e -> 1)
        epsv = sp.symbols('eps', real=True)
        return dict(G=sp.expand(corner_Q().subs(e, 1 - epsv)), box=[epsv, Pz, Qz, mu], root=th)
    if chart in ('KBRe', 'KBLe'):
        # side charts in eps = 1 - e, centred at f = 3 pi/2 (R) or pi/2 (L): G(r) = P(f)(1 + r^2)^6, f = f_c + 2 atan r
        # (kq_run --quad, C3_CORE.md section 7)
        epsv = sp.symbols('eps', real=True)
        Gs = sp.expand(corner_Q('right' if chart == 'KBRe' else 'left').subs(e, 1 - epsv))
        return dict(G=Gs, box=[epsv, Pz, Qz, mu], root=th)
    if chart == 'KBPe':
        # perihelion chart in eps = 1 - e: root tp = tan(f/2) = -1/th, G(tp) = tp^12 Q(-1/tp) = P(f)(1 + tp^2)^6, i.e. the
        # coefficients of KBQe reversed with signs (-1)^k (kq_run --peri, C3_CORE.md section 7)
        epsv = sp.symbols('eps', real=True)
        Pq = sp.Poly(sp.expand(corner_Q().subs(e, 1 - epsv)), th)
        assert Pq.degree() <= 12
        G = sp.expand(sum(cf * (-1)**k * th**(12 - k) for (k,), cf in Pq.terms()))
        return dict(G=G, box=[epsv, Pz, Qz, mu], root=th)
    if chart in ('PhiP', 'PhiM'):
        sg = 1 if chart == 'PhiP' else -1
        Qx = corner_Q()
        sub = {th: sg * Y0 * nu, Qz: nu * a, Pz: 1 - nu**2 * K, e: 1 - nu**2 * E, mu: nu**3 * m}
        W = sp.Poly(sp.expand(Qx.subs(sub)), nu)
        assert min(mon[0] for mon in W.monoms()) == 10
        G = sp.expand(W.as_expr() / nu**10)
        return dict(G=G, box=[nu, a, K, E], root=m)
    # level 2: parent chart Ca (ratio K) or CK (ratio a)
    parent = 'Ca' if chart.startswith('L2a') else 'CK'
    G1, box1 = level1(parent)
    rvar = K if parent == 'Ca' else a
    kind = chart[3:]
    if kind == 'E':
        sub = {E: E2 * eta**2, m: M2 * eta**2 * mp, y: Y2 * eta * Yv}; box = [nu, rvar, eta, mp]; root = Yv
    elif kind == 'm':
        sub = {m: M2 * eta**2, E: E2 * eta**2 * Ep, y: Y2 * eta * Yv}; box = [nu, rvar, eta, Ep]; root = Yv
    elif kind in ('Yp', 'Ym'):
        sg = 1 if kind == 'Yp' else -1
        sub = {y: sg * Y2 * eta, E: E2 * eta**2 * Ep, m: M2 * eta**2 * mp}; box = [nu, rvar, eta, Ep]; root = mp
    else:
        raise ValueError(chart)
    W = sp.Poly(sp.expand(G1.subs(sub)), eta)
    low = min(mon[0] for mon in W.monoms())
    assert low == 6, (chart, low)
    G = sp.expand(W.as_expr() / eta**6)
    return dict(G=G, box=box, root=root)


def phys(chart, x, r):
    """map (box point x, root value r) -> (e, Pz, Qz, mu, th, norm) with Q(th) = norm * G."""
    if chart == 'FE':
        return 1 - x[0]**2, x[2], x[3], x[0]**2 * x[1], x[0] * r, x[0]**6
    if chart == 'FM':
        return 1 - x[0]**2 * x[1], x[2], x[3], x[0]**2, x[0] * r, x[0]**6
    if chart in ('FIP', 'FIM'):
        sg = 1 if chart == 'FIP' else -1
        return 1 - x[0]**2 * x[1], x[2], x[3], x[0]**2 * r, sg * YF * x[0], x[0]**6
    if chart in ('Ca', 'CK', 'CE', 'Cm'):
        n_ = x[0]
        if chart == 'Ca':
            qz, kap, ep, mm = n_, n_**2 * x[1], n_**2 * x[2], n_**3 * x[3]
        elif chart == 'CK':
            qz, kap, ep, mm = n_ * x[1], n_**2, n_**2 * x[2], n_**3 * x[3]
        elif chart == 'CE':
            qz, kap, ep, mm = n_ * x[1], n_**2 * x[2], n_**2, n_**3 * x[3]
        else:
            qz, kap, ep, mm = n_ * x[1], n_**2 * x[2], n_**2 * x[3], n_**3
        return 1 - ep, 1 - kap, qz, mm, n_ * r, n_**10
    if chart == 'KBQ':
        return x[0], x[1], x[2], x[3], r, 1.0
    if chart == 'KBQe':
        return 1 - x[0], x[1], x[2], x[3], r, 1.0
    if chart == 'KBPe':
        return 1 - x[0], x[1], x[2], x[3], -1.0 / r, 1.0 / r**12
    if chart == 'KBRe':          # th = tan((f - pi)/2) = (1 + r)/(1 - r), Q(th) = (2/(1 - r)^2)^6 G(r)
        return 1 - x[0], x[1], x[2], x[3], (1 + r) / (1 - r), 64.0 / (1 - r)**12
    if chart == 'KBLe':          # th = (r - 1)/(1 + r), Q(th) = (2/(1 + r)^2)^6 G(r)
        return 1 - x[0], x[1], x[2], x[3], (r - 1) / (1 + r), 64.0 / (1 + r)**12
    if chart in ('PhiP', 'PhiM'):
        sg = 1 if chart == 'PhiP' else -1
        n_, a_, K_, E_ = x
        return 1 - n_**2 * E_, 1 - n_**2 * K_, n_ * a_, n_**3 * r, sg * Y0 * n_, n_**10
    n_, rv, et, w4 = x
    kind = chart[3:]
    if kind == 'E':
        E_, m_, y_ = E2 * et**2, M2 * et**2 * w4, Y2 * et * r
    elif kind == 'm':
        E_, m_, y_ = E2 * et**2 * w4, M2 * et**2, Y2 * et * r
    else:
        sg = 1 if kind == 'Yp' else -1
        E_, m_, y_ = E2 * et**2 * w4, M2 * et**2 * r, sg * Y2 * et
    if chart.startswith('L2a'):
        qz, kap = n_, n_**2 * rv
    else:
        qz, kap = n_ * rv, n_**2
    return 1 - n_**2 * E_, 1 - kap, qz, n_**3 * m_, n_ * y_, n_**10 * et**6


def P_direct(ee, pz, qz, mm, f):
    zeta = pz * math.cos(f) + qz * math.sin(f)
    dzeta = -pz * math.sin(f) + qz * math.cos(f)
    B = -2 * zeta * dzeta + 2 * ee * (math.sin(f) - qz * zeta)
    D = 1 + ee * math.cos(f)
    k = 2 * ee * (1 + ee)
    return mm**2 * k**2 * math.sin(f)**2 * (1 - zeta**2) - D**2 * B**2


def main(chart):
    D = chart_def(chart)
    G, box, root = D['G'], D['box'], D['root']
    GP = sp.Poly(G, root)
    deg = GP.degree()
    coeffs = [sp.expand(GP.coeff_monomial(root**k)) for k in range(deg + 1)]
    for c in coeffs:
        assert not (c.free_symbols - set(box)), (chart, c.free_symbols - set(box))
    polys = {}
    for k, c in enumerate(coeffs):
        polys[(0, k)] = c
        for v in range(4):
            polys[(1 + v, k)] = sp.expand(sp.diff(c, box[v]))
    names = [str(s) for s in box]
    maxdeg = {n: 0 for n in names}
    nterms = 0
    for ex in polys.values():
        if ex == 0:
            continue
        Pp = sp.Poly(ex, *box)
        nterms = max(nterms, len(Pp.terms()))
        for mon in Pp.monoms():
            for n, d in zip(names, mon):
                maxdeg[n] = max(maxdeg[n], d)

    def emit(fn, use_abs):
        L = ['@njit(cache=True)', f'def {fn}(x, out):']
        for i, n in enumerate(names):
            if maxdeg[n] >= 1:
                L.append(f'    {n}1 = x[{i}]')
            for d in range(2, maxdeg[n] + 1):
                L.append(f'    {n}{d} = {n}{d-1} * {n}1')
        for (dv, k), ex in polys.items():
            if ex == 0:
                L.append(f'    out[{dv}, {k}] = 0.0')
                continue
            terms = []
            for mon, coef in sp.Poly(ex, *box).terms():
                cval = float(abs(coef)) if use_abs else float(coef)
                mono = ' * '.join(f'{n}{d}' for n, d in zip(names, mon) if d > 0) or '1.0'
                terms.append(f'{cval!r} * {mono}')
            L.append(f'    out[{dv}, {k}] = ' + ' + '.join(terms))
        L.append('')
        return L

    src = [f'"""Generated by gen_c3.py, chart {chart}: G = sum_k c_k(x) {root}^k, box vars {names}, degree {deg};'
           f' max terms per polynomial {nterms}."""', 'from numba import njit', '', f'DEG = {deg}',
           f'BOX = {names!r}', f'ROOT = {str(root)!r}', f'MAXTERMS = {nterms}', '']
    src += emit('coef_val', False)
    src += emit('coef_abs', True)
    open(f'nbcoef_c3_{chart}.py', 'w').write('\n'.join(src) + '\n')
    # self-test with mpmath-free double check at moderate parameters (norm not tiny)
    fG = sp.lambdify([root] + box, G, 'math')
    random.seed(2)
    worst = 0.0
    for _ in range(60):
        x = [random.uniform(0.05, 0.37)] + [random.uniform(0, 1.05) for _ in range(3)]
        if chart == 'KBQe':
            x = [random.uniform(0, 0.01), random.uniform(0.4, 0.99), random.uniform(0.1, 0.9), random.uniform(0, 0.05)]
            r = random.uniform(-0.5, 0.5)
        elif chart in ('KBRe', 'KBLe'):
            x = [random.uniform(0, 0.2), random.uniform(0.0, 0.99), random.uniform(0.0, 0.9), random.uniform(0, 0.5)]
            r = random.uniform(-0.6, 0.6)
        elif chart == 'KBPe':
            x = [random.uniform(0, 0.2), random.uniform(0.0, 0.99), random.uniform(0.0, 0.9), random.uniform(0, 0.5)]
            r = random.choice([-1.0, 1.0]) * random.uniform(0.05, 0.8)
        elif chart == 'KBQ':
            x = [random.uniform(0.99, 1.0), random.uniform(0.4, 0.99), random.uniform(0.1, 0.9), random.uniform(0, 0.05)]
            r = random.uniform(-0.5, 0.5)
        elif chart in ('Ca', 'CK', 'CE', 'Cm'):
            r = random.uniform(-Y0, Y0)
        elif chart in ('FE', 'FM', 'FIP', 'FIM'):
            pz_ = random.uniform(0, 0.95)
            x = [random.uniform(0.1, 0.3), random.uniform(0, 1.0), pz_, random.uniform(0, math.sqrt(1 - pz_ * pz_))]
            r = random.uniform(-YF, YF) if chart in ('FE', 'FM') else random.uniform(0, 1)
        elif chart in ('PhiP', 'PhiM'):
            r = random.uniform(0, 1)
        elif chart[3:] in ('E', 'm'):
            x[2] = random.uniform(0.05, 1.0); r = random.uniform(-1.05, 1.05)
        else:
            x[2] = random.uniform(0.05, 1.0); r = random.uniform(0, 1.05)
        ee, pz, qz, mm, thv, norm = phys(chart, x, r)
        f = math.pi + 2 * math.atan(thv)
        ref = P_direct(ee, pz, qz, mm, f) * (1 + thv * thv)**6 / norm
        val = fG(r, *x)
        worst = max(worst, abs(val - ref) / (1 + abs(ref)))
    print(f'nbcoef_c3_{chart}.py: root {root}, degree {deg}, box {names}, max degrees {maxdeg}, max terms {nterms};'
          f' self-test max rel. error {worst:.2e}', flush=True)
    assert worst < 1e-8, "self-test failed"


if __name__ == '__main__':
    for ch in sys.argv[1:]:
        main(ch)

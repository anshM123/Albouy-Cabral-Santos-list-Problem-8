"""Blow-up charts for the near-coincident corner (D1a) and its degenerate direction (D1b).

Physical model (a-chart, semi-major axis 1):  r = 1 - e cos u,  Z = p (cos u - e) + q b sin u   (p = l sin w, q = l cos w),
rho^2 = r^2 - Z^2,  A = (r^2)' = 2 e sin u r,  Phi = A^2 rho^2,  Psi = ((rho^2)')^2,  f = (1 - t) Phi - t Psi,  t = s^2/(1+s^2).
Write tau = t - 1/2, so f = (Phi - Psi)/2 - tau (Phi + Psi).

Corner D1a: e -> 0, theta := l^2/e -> 0, tau -> 0 (all w).  In the coordinates (e, theta, tau) the function f/e^2 is
  f/e^2 = 2 sin u [theta sin 2(u+w) - 4 tau sin u - e sin 2u] + (second order),
homogeneous of degree one, so the corner is resolved by the cone charts
  B1   (e dominant):      p = e ph, q = e qh (theta = e (ph^2 + qh^2)), tau = e tauh;   F = f / e^3        [affine in tauh]
  B2   (theta dominant):  e = th eh, l^2 = e th, tau = th tauh;                         F = f / (e^2 th)   [affine in tauh]
  B3p/B3m (tau = +-eta):  e = eta eh, theta = eta thh;                                   F = f / (e^2 eta)  [4 mean-value vars]
The leading term of B1 vanishes identically only at ph = 0, qh = 1, tauh = 0, e = 0 (direction D1b: l = e, w = 0, s = 1).
Level-2 cone charts around that point, in (e, ph, a = qh - 1, tauh), with F = F_B1 / nu (nu = dominant variable):
  L2E  e = nu, ph = nu P, a = nu A, tauh = nu T          [affine in T]
  L2P  ph = nu, e = nu E, a = nu A, tauh = nu T          [affine in T]
  L2Ap/L2Am  a = +-nu, e = nu E, ph = nu P, tauh = nu T  [affine in T]
  L2Tp/L2Tm  tauh = +-nu, e = nu E, ph = nu P, a = nu A  [4 mean-value vars]
For the affine charts the box variable x3 in [0, 1] parametrises the affine variable V = V0 + (V1 - V0) x3 and
  F = (1 - x3) F_A - x3 F_B  with  F_A = F|_{V = V0},  F_B = -F|_{V = V1}.
Writes nbcoef_<chart>.py with coef_val(a, out) / coef_abs(a, out); out[0, Fi, part, k] = coefficient of
cos(k u) (part 0) / sin(k u) (part 1); out[1 + v, Fi, part, k] = derivative w.r.t. box variable v (chain rule).

usage: python gen_bu.py CHART [CHART ...]
"""
import sys
import random
import math
import sympy as sp

z = sp.symbols('z')
cu = (z + 1 / z) / 2
su = (z - 1 / z) / (2 * sp.I)


def ddu(expr):
    return sp.expand(sp.I * z * sp.diff(expr, z))


# ---------- range parameters of the affine variables (must match the run scripts / coverage lemma) ----------
T_B1 = 2.0        # tauh in [-T_B1, T_B1]  (cone condition needs 1.95)
T_B2 = 0.6825     # tauh in [-T_B2, T_B2]
K_E = 1.05        # level-2 E chart: P in [0, K], A in [-K, K], T in [-K, K] (scaled ratios, see SC2)
SD2 = 5e-3       # D2b corner cube scale (all four coordinates)
SC2 = (0.01, 0.25, 0.25, 0.5)   # level-2 scale vector for (e, ph, a, tauh); cones are taken in scaled coordinates


def B1_GG(e, ph, qh, b):
    """G1 = (Phi - Psi)/e^3, G2 = (Phi + Psi)/e^2 for p = e ph, q = e qh (exact polynomial identities)."""
    r = 1 - e * cu
    At = 2 * su * r                                   # A / e
    zh = ph * (cu - e) + qh * b * su                  # Z / e
    zh2 = sp.expand(zh**2)
    dzh2 = ddu(zh2)
    G1 = sp.expand(At**2 * (-2 * cu + e * cu**2 - e * zh2) + 2 * At * dzh2 - e * dzh2**2)
    G2 = sp.expand(At**2 * (r**2 - e**2 * zh2) + (At - e * dzh2)**2)
    return G1, G2


def reduce_b(expr, b, e):
    """reduce modulo b^2 = 1 - e^2 (so that b appears at most linearly)."""
    return sp.expand(sp.rem(sp.Poly(sp.expand(expr), b), sp.Poly(b**2 - (1 - e**2), b)).as_expr())


def chart_def(chart):
    """returns dict(atoms, DER, F (list of 1 or 2 sympy exprs), nvar, phys(x) -> (e, p, q, t, norm), affine)"""
    if chart == 'B1':
        e, ph, qh, b, ib = sp.symbols('e ph qh b ib', real=True)
        G1, G2 = B1_GG(e, ph, qh, b)
        V0, V1 = -T_B1, T_B1
        FA = sp.expand(G1 / 2 - V0 * G2); FB = sp.expand(-(G1 / 2 - V1 * G2))
        DER = {e: [1, 0, 0], ph: [0, 1, 0], qh: [0, 0, 1], b: [-e * ib, 0, 0], ib: [e * ib**3, 0, 0]}

        def phys(x):
            ee, p_, q_, x3 = x
            tauh = V0 + (V1 - V0) * x3
            return ee, ee * p_, ee * q_, 0.5 + ee * tauh, ee**3
        return dict(atoms=[e, ph, qh, b, ib], DER=DER, F=[FA, FB], phys=phys, affine=True,
                    atom_fn='B1')
    if chart == 'B2':
        th, eh, b, ib, sw, cw = sp.symbols('th eh b ib sw cw', real=True)
        e = th * eh
        r = 1 - e * cu
        At = 2 * su * r
        zt = sw * (cu - e) + b * cw * su                  # Z / l
        zt2 = sp.expand(zt**2); dzt2 = ddu(zt2)
        G1 = sp.expand(At**2 * eh * (-2 * cu + e * cu**2 - th * zt2) + 2 * At * dzt2 - th * dzt2**2)
        G2 = sp.expand(At**2 * (r**2 - e * th * zt2) + (At - th * dzt2)**2)
        V0, V1 = -T_B2, T_B2
        FA = sp.expand(G1 / 2 - V0 * G2); FB = sp.expand(-(G1 / 2 - V1 * G2))
        DER = {th: [1, 0, 0], eh: [0, 1, 0], b: [-th * eh**2 * ib, -th**2 * eh * ib, 0],
               ib: [th * eh**2 * ib**3, th**2 * eh * ib**3, 0], sw: [0, 0, cw], cw: [0, 0, -sw]}

        def phys(x):
            t_, e_h, w, x3 = x
            ee = t_ * e_h; L = ee * t_; l = math.sqrt(L)
            tauh = V0 + (V1 - V0) * x3
            return ee, l * math.sin(w), l * math.cos(w), 0.5 + t_ * tauh, ee**2 * t_
        return dict(atoms=[th, eh, b, ib, sw, cw], DER=DER, F=[FA, FB], phys=phys, affine=True, atom_fn='B2')
    if chart in ('B3p', 'B3m'):
        sg = 1 if chart == 'B3p' else -1
        eta, eh, thh, b, ib, sw, cw = sp.symbols('eta eh thh b ib sw cw', real=True)
        e = eta * eh; th = eta * thh
        r = 1 - e * cu
        At = 2 * su * r
        zt = sw * (cu - e) + b * cw * su
        zt2 = sp.expand(zt**2); dzt2 = ddu(zt2)
        G1 = sp.expand(At**2 * eh * (-2 * cu + e * cu**2 - th * zt2) + 2 * thh * At * dzt2 - eta * thh**2 * dzt2**2)
        G2 = sp.expand(At**2 * (r**2 - e * th * zt2) + (At - th * dzt2)**2)
        F = sp.expand(G1 / 2 - sg * G2)
        DER = {eta: [1, 0, 0, 0], eh: [0, 1, 0, 0], thh: [0, 0, 1, 0],
               b: [-eta * eh**2 * ib, -eta**2 * eh * ib, 0, 0], ib: [eta * eh**2 * ib**3, eta**2 * eh * ib**3, 0, 0],
               sw: [0, 0, 0, cw], cw: [0, 0, 0, -sw]}

        def phys(x):
            et, e_h, t_h, w = x
            ee = et * e_h; th_ = et * t_h; L = ee * th_; l = math.sqrt(L)
            return ee, l * math.sin(w), l * math.cos(w), 0.5 + sg * et, ee**2 * et
        return dict(atoms=[eta, eh, thh, b, ib, sw, cw], DER=DER, F=[F], phys=phys, affine=False, atom_fn=chart)
    if chart.startswith('L2'):
        # F_B1 = G1/2 - tauh G2 with b reduced (b at most linear; the b-linear part carries a factor ph)
        e0, ph0, qh0, b = sp.symbols('e0 ph0 qh0 b', real=True)
        G1, G2 = B1_GG(e0, ph0, qh0, b)
        G1 = reduce_b(G1, b, e0); G2 = reduce_b(G2, b, e0)
        nu, E, P, A, T, ib = sp.symbols('nu E P A T ib', real=True)
        s0, s1, s2, s3 = SC2
        # scaled level-2 variables: e = s0 X0, ph = s1 X1, a = s2 X2, tauh = s3 X3; the dominant scaled variable is nu
        if chart == 'L2E':
            X = (nu, P, A, T); dom = 0
        elif chart == 'L2P':
            X = (E, nu, A, T); dom = 1
        elif chart in ('L2Ap', 'L2Am'):
            sgn = 1 if chart == 'L2Ap' else -1
            X = (E, P, sgn * nu, T); dom = 2
        else:
            sgn = 1 if chart == 'L2Tp' else -1
            X = (E, P, A, sgn * nu); dom = 3
        # actual values: scaled ratio variables are multiplied by nu (except the dominant one)
        Xv = [nu if i == dom else X[i] * nu for i in range(4)]
        if dom in (2, 3):
            Xv[dom] = X[dom]
        ev = s0 * Xv[0]; phv = s1 * Xv[1]; qhv = 1 + s2 * Xv[2]; tauv = s3 * Xv[3]
        sub = {e0: ev, ph0: phv, qh0: qhv}
        g1 = sp.expand(G1.subs(sub)); g2 = sp.expand(G2.subs(sub))
        Fb1 = sp.expand(g1 / 2 - tauv * g2)
        c0 = sp.expand(Fb1.subs(nu, 0))
        assert c0 == 0, chart + ': F_B1 not divisible by nu'
        Fn = sp.expand(Fb1 / nu)
        assert not any(sp.Poly(sp.expand(Fn * z**6), nu, z).as_expr().has(1 / nu) for _ in [0]), chart
        if chart in ('L2Tp', 'L2Tm'):
            # 4 mean-value variables: nu, E, P, A
            atoms = [nu, E, P, A, b, ib]
            DER = {nu: [1, 0, 0, 0], E: [0, 1, 0, 0], P: [0, 0, 1, 0], A: [0, 0, 0, 1]}
            # e = s0 E nu
            DER[b] = [-(s0 * E * nu) * ib * s0 * E, -(s0 * E * nu) * ib * s0 * nu, 0, 0]
            DER[ib] = [(s0 * E * nu) * ib**3 * s0 * E, (s0 * E * nu) * ib**3 * s0 * nu, 0, 0]
            F = [Fn]

            def phys(x, sgn=sgn):
                n_, E_, P_, A_ = x
                ee = s0 * E_ * n_; ph_ = s1 * P_ * n_; qh_ = 1 + s2 * A_ * n_; tauh = s3 * sgn * n_
                return ee, ee * ph_, ee * qh_, 0.5 + ee * tauh, ee**3 * n_
            return dict(atoms=atoms, DER=DER, F=F, phys=phys, affine=False, atom_fn=chart)
        # affine charts: variable T in [-K, K] -> x3
        V0, V1 = -K_E, K_E
        FA = sp.expand(Fn.subs(T, V0)); FB = sp.expand(-Fn.subs(T, V1))
        if chart == 'L2E':
            atoms = [nu, P, A, b, ib]
            DER = {nu: [1, 0, 0], P: [0, 1, 0], A: [0, 0, 1], b: [-(s0 * nu) * ib * s0, 0, 0],
                   ib: [(s0 * nu) * ib**3 * s0, 0, 0]}

            def phys(x):
                n_, P_, A_, x3 = x
                T_ = V0 + (V1 - V0) * x3
                ee = s0 * n_; ph_ = s1 * P_ * n_; qh_ = 1 + s2 * A_ * n_; tauh = s3 * T_ * n_
                return ee, ee * ph_, ee * qh_, 0.5 + ee * tauh, ee**3 * n_
        elif chart == 'L2P':
            atoms = [nu, E, A, b, ib]
            DER = {nu: [1, 0, 0], E: [0, 1, 0], A: [0, 0, 1],
                   b: [-(s0 * E * nu) * ib * s0 * E, -(s0 * E * nu) * ib * s0 * nu, 0],
                   ib: [(s0 * E * nu) * ib**3 * s0 * E, (s0 * E * nu) * ib**3 * s0 * nu, 0]}

            def phys(x):
                n_, E_, A_, x3 = x
                T_ = V0 + (V1 - V0) * x3
                ee = s0 * E_ * n_; ph_ = s1 * n_; qh_ = 1 + s2 * A_ * n_; tauh = s3 * T_ * n_
                return ee, ee * ph_, ee * qh_, 0.5 + ee * tauh, ee**3 * n_
        else:
            atoms = [nu, E, P, b, ib]
            DER = {nu: [1, 0, 0], E: [0, 1, 0], P: [0, 0, 1],
                   b: [-(s0 * E * nu) * ib * s0 * E, -(s0 * E * nu) * ib * s0 * nu, 0],
                   ib: [(s0 * E * nu) * ib**3 * s0 * E, (s0 * E * nu) * ib**3 * s0 * nu, 0]}

            def phys(x, sgn=sgn):
                n_, E_, P_, x3 = x
                T_ = V0 + (V1 - V0) * x3
                ee = s0 * E_ * n_; ph_ = s1 * P_ * n_; qh_ = 1 + s2 * sgn * n_; tauh = s3 * T_ * n_
                return ee, ee * ph_, ee * qh_, 0.5 + ee * tauh, ee**3 * n_
        return dict(atoms=atoms, DER=DER, F=[FA, FB], phys=phys, affine=True, atom_fn=chart)
    if chart.startswith('D2'):
        # Corner D2b inside C1p (e -> 0, s ~ e, l -> 1, w -> 0).  C1p function: F_C1p = (1-tau) Phi/e^2 - tau (Phi + Psi)
        # (t = e^2 tau/(1-tau)).  Coordinates (e, p, kap = 1 - q, tp = tau - 1/2), all weight 1, scale SD2 each.
        e0, p0, k0, b = sp.symbols('e0 p0 k0 b', real=True)
        r = 1 - e0 * cu
        At = 2 * su * r
        q0 = 1 - k0
        Zx = p0 * (cu - e0) + q0 * b * su
        Z2 = reduce_b(sp.expand(Zx**2), b, e0)
        dZ2 = ddu(Z2)
        rho2 = sp.expand(r**2 - Z2)
        G0 = reduce_b(sp.expand(At**2 * rho2 - e0**2 * At**2 * rho2 - (e0 * At - dZ2)**2), b, e0)   # Phi/e^2 - Phi - Psi
        G1 = reduce_b(sp.expand(At**2 * rho2 + e0**2 * At**2 * rho2 + (e0 * At - dZ2)**2), b, e0)   # Phi/e^2 + Phi + Psi
        nu, E, P, K, T, ib = sp.symbols('nu E P K T ib', real=True)
        sd = SD2
        if chart in ('D2E', 'D2E2'):
            X = (nu, P * nu, K * nu, T * nu); atoms = [nu, P, K, b, ib]; dom = 'E'
        elif chart == 'D2P':
            X = (E * nu, nu, K * nu, T * nu); atoms = [nu, E, K, b, ib]; dom = 'P'
        elif chart == 'D2K':
            X = (E * nu, P * nu, nu, T * nu); atoms = [nu, E, P, b, ib]; dom = 'K'
        else:
            sgn = 1 if chart == 'D2Tp' else -1
            X = (E * nu, P * nu, K * nu, sgn * nu); atoms = [nu, E, P, K, b, ib]; dom = 'T'
        ev, pv, kv, tv = sd * X[0], sd * X[1], sd * X[2], sd * X[3]
        sub = {e0: ev, p0: pv, k0: kv}
        Fb = sp.expand(sp.expand(G0.subs(sub)) / 2 - tv * sp.expand(G1.subs(sub)))
        assert sp.expand(Fb.subs(nu, 0)) == 0, chart + ': not divisible by nu'
        Fn = sp.expand(Fb / nu)
        Ex = sd * (nu if dom == 'E' else E * nu)               # e as an expression in the atoms
        dE = [sp.diff(Ex, v) for v in ([nu, P, K] if dom == 'E' else [nu, E, K] if dom == 'P' else [nu, E, P] if dom == 'K' else [nu, E, P, K])]
        DER = {a_: [1 if a_ == v else 0 for v in atoms[:len(dE)]] for a_ in atoms[:len(dE)]}
        DER[b] = [sp.expand(-Ex * ib * d) for d in dE]
        DER[ib] = [sp.expand(Ex * ib**3 * d) for d in dE]

        def phys(x, dom=dom, chart=chart):
            if dom == 'T':
                n_, E_, P_, K_ = x; sg = 1 if chart == 'D2Tp' else -1
                ee, pp, kk, tp = sd * E_ * n_, sd * P_ * n_, sd * K_ * n_, sd * sg * n_
            else:
                n_, a1, a2, x3 = x
                T_ = (1.04 + 0.24 * x3) if chart == 'D2E2' else (-1.05 + 2.1 * x3)
                if dom == 'E':
                    ee, pp, kk = sd * n_, sd * a1 * n_, sd * a2 * n_
                elif dom == 'P':
                    ee, pp, kk = sd * a1 * n_, sd * n_, sd * a2 * n_
                else:
                    ee, pp, kk = sd * a1 * n_, sd * a2 * n_, sd * n_
                tp = sd * T_ * n_
            qq = 1 - kk; tau = 0.5 + tp
            t = ee * ee * tau / (1 - tau)
            # F * nu = (1 - tau) f / e^2  with f = (1 - t) Phi - t Psi
            return ee, pp, qq, t, ee * ee * n_ / (1 - tau)
        if dom == 'T':
            return dict(atoms=atoms, DER=DER, F=[Fn], phys=phys, affine=False, atom_fn=chart, d2=True)
        # D2E2: D2E extended to T = tau'/e in [1.04, 1.28] (covers the part of D2Tp's degenerate box with e/tau' < 1/1.05)
        V0, V1 = (sp.Rational(104, 100), sp.Rational(128, 100)) if chart == 'D2E2' else (-1.05, 1.05)
        FA = sp.expand(Fn.subs(T, V0)); FB = sp.expand(-Fn.subs(T, V1))
        return dict(atoms=atoms, DER=DER, F=[FA, FB], phys=phys, affine=True, atom_fn=chart, d2=True)
    raise SystemExit('unknown chart ' + chart)


def real_form(F):
    P = sp.Poly(sp.expand(F * z**6), z)
    assert P.degree() <= 12, P.degree()
    c = {k - 6: sp.expand(P.coeff_monomial(z**k)) for k in range(13)}
    a = [sp.expand(sp.re(c[0]))] + [sp.expand(2 * sp.re(c[k])) for k in range(1, 7)]
    s = [sp.Integer(0)] + [sp.expand(-2 * sp.im(c[k])) for k in range(1, 7)]
    assert sp.expand(sp.im(c[0])) == 0
    return a, s


def f_phys(e, p, q, t, u):
    b = math.sqrt(1 - e * e)
    r = 1 - e * math.cos(u)
    Z = p * (math.cos(u) - e) + q * b * math.sin(u)
    Zp = -p * math.sin(u) + q * b * math.cos(u)
    A = 2 * e * math.sin(u) * r
    rho2 = r * r - Z * Z
    d = A - 2 * Z * Zp
    return (1 - t) * A * A * rho2 - t * d * d


def atom_values(chart, x):
    if chart.startswith('D2'):
        if chart in ('D2Tp', 'D2Tm'):
            n_, E_, P_, K_ = x; ee = SD2 * E_ * n_; b = math.sqrt(1 - ee * ee); return [n_, E_, P_, K_, b, 1 / b]
        n_, a1, a2, _ = x
        ee = SD2 * n_ if chart in ('D2E', 'D2E2') else SD2 * a1 * n_
        b = math.sqrt(1 - ee * ee)
        return [n_, a1, a2, b, 1 / b]
    """numeric atom values for the self-test (must match atoms_box in make_engines_bu.py)."""
    if chart == 'B1':
        ee, ph, qh, _ = x; b = math.sqrt(1 - ee * ee); return [ee, ph, qh, b, 1 / b]
    if chart == 'B2':
        th, eh, w, _ = x; ee = th * eh; b = math.sqrt(1 - ee * ee); return [th, eh, b, 1 / b, math.sin(w), math.cos(w)]
    if chart in ('B3p', 'B3m'):
        et, eh, thh, w = x; ee = et * eh; b = math.sqrt(1 - ee * ee)
        return [et, eh, thh, b, 1 / b, math.sin(w), math.cos(w)]
    s0 = SC2[0]
    if chart == 'L2E':
        n_, P_, A_, _ = x; ee = s0 * n_; b = math.sqrt(1 - ee * ee); return [n_, P_, A_, b, 1 / b]
    if chart == 'L2P':
        n_, E_, A_, _ = x; ee = s0 * E_ * n_; b = math.sqrt(1 - ee * ee); return [n_, E_, A_, b, 1 / b]
    if chart in ('L2Ap', 'L2Am'):
        n_, E_, P_, _ = x; ee = s0 * E_ * n_; b = math.sqrt(1 - ee * ee); return [n_, E_, P_, b, 1 / b]
    n_, E_, P_, A_ = x; ee = s0 * E_ * n_; b = math.sqrt(1 - ee * ee); return [n_, E_, P_, A_, b, 1 / b]


def sample_box(chart):
    if chart == 'B1':
        return [random.uniform(0.001, 0.07), random.uniform(0, 1.75), random.uniform(0, 1.75), random.random()]
    if chart == 'B2':
        return [random.uniform(0.001, 0.2), random.uniform(0, 0.37), random.uniform(0, math.pi / 2), random.random()]
    if chart in ('B3p', 'B3m'):
        return [random.uniform(0.001, 0.13), random.uniform(0, 0.57), random.uniform(0, 1.62), random.uniform(0, math.pi / 2)]
    return [random.uniform(0.001, 1.0), random.uniform(0, 1.05), random.uniform(-1.05, 1.05), random.random()] \
        if chart in ('L2E',) else [random.uniform(0.001, 1.0), random.uniform(0, 1.05), random.uniform(-1.05, 1.05),
                                   random.random()]


def main(chart):
    D = chart_def(chart)
    atoms = D['atoms']; DER = D['DER']; Fs = D['F']
    nv = len(next(iter(DER.values())))
    polys = {}
    for fi, F in enumerate(Fs):
        a, s = real_form(F)
        for part, arr in ((0, a), (1, s)):
            for k in range(7):
                ex = arr[k]
                assert not (ex.free_symbols - set(atoms)), (chart, ex.free_symbols - set(atoms))
                polys[(0, fi, part, k)] = ex
                for v in range(nv):
                    dd = sum(sp.diff(ex, at) * DER[at][v] for at in atoms if at in DER)
                    polys[(1 + v, fi, part, k)] = sp.expand(dd)
    # every atom that occurs must have a DER entry
    used = set().union(*[ex.free_symbols for ex in polys.values() if ex != 0])
    for at in used:
        assert at in DER, (chart, 'no derivative rule for', at)
    names = [str(a) for a in atoms]
    maxdeg = {n: 0 for n in names}
    nterms = 0
    for ex in polys.values():
        if ex == 0:
            continue
        P = sp.Poly(ex, *atoms)
        nterms = max(nterms, len(P.terms()))
        for mon in P.monoms():
            for n, d in zip(names, mon):
                maxdeg[n] = max(maxdeg[n], d)

    def emit(fn, use_abs):
        L = ['@njit(cache=True)', f'def {fn}(a, out):']
        for i, n in enumerate(names):
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
            for mon, coef in sp.Poly(ex, *atoms).terms():
                cval = float(abs(coef)) if use_abs else float(coef)
                mono = ' * '.join(f'{n}{d}' for n, d in zip(names, mon) if d > 0) or '1.0'
                terms.append(f'{cval!r} * {mono}')
            L.append(f'    out[{idx}] = ' + ' + '.join(terms))
        L.append('')
        return L

    src = [f'"""Generated by gen_bu.py, chart {chart}. atoms = {names}; nF = {len(Fs)}; nvar = {nv};'
           f' max terms per polynomial = {nterms}."""', 'from numba import njit', '',
           f'ATOMS = {names!r}', f'NF = {len(Fs)}', f'NVAR = {nv}', f'MAXTERMS = {nterms}', '']
    src += emit('coef_val', False)
    src += emit('coef_abs', True)
    open(f'nbcoef_{chart}.py', 'w').write('\n'.join(src) + '\n')
    # ---- self-test: chart function * normalisation == physical f, at random points ----
    fl = [sp.lambdify([z] + atoms, F, 'math') for F in Fs]
    random.seed(1)
    worst = 0.0
    for _ in range(40):
        x = sample_box(chart)
        ee, p, q, t, norm = D['phys'](x)
        av = atom_values(chart, x)
        for _ in range(3):
            u = random.uniform(0, 2 * math.pi)
            zz = complex(math.cos(u), math.sin(u))
            if D['affine']:
                x3 = x[3]
                val = (1 - x3) * fl[0](zz, *av) - x3 * fl[1](zz, *av)
            else:
                val = fl[0](zz, *av)
            ref = f_phys(ee, p, q, t, u) / norm
            scale = 1.0 + abs(ref)
            worst = max(worst, abs(complex(val).real - ref) / scale, abs(complex(val).imag) / scale)
    print(f'nbcoef_{chart}.py written; atoms {names}; max degrees {maxdeg}; max terms {nterms}; '
          f'self-test max rel. error {worst:.2e}', flush=True)
    assert worst < 1e-9, 'self-test failed'


if __name__ == '__main__':
    for ch in sys.argv[1:]:
        main(ch)

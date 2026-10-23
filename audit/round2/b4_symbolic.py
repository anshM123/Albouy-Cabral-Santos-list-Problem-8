"""AUDIT-2 b4: symbolic checks for gen_bu.py (independent derivation with real-trig dicts, as in AUDIT 1 a1).

1. B1 identities: Phi - Psi == e^3 G1, Phi + Psi == e^2 G2 for p = e ph, q = e qh (G1, G2 as written in gen_bu).
2. b-reduction: G1, G2 mod b^2 = 1 - e^2 are R0 + b R1 with ph | R1 (so nu | F_B1 needs only R0 at the D1b point).
3. Leading term F0 = (G1/2 - tauh G2)|_{e=0, b=1} = 2 sin u [(qh^2 - ph^2 - 1) sin 2u + 2 ph qh cos 2u - 4 tauh sin u];
   it vanishes identically only at ph = 0, qh = 1, tauh = 0 (for ph, qh >= 0).
4. CT expansion: f/e^2 = 2 sin u [theta sin 2(u+w) - 4 tau sin u - e sin 2u] + O(2) (first order of e*F0).
5. Lemma E0 (e = 0): Pol = -S l^4 sin^2 2(u+w); critical points of the circle-circle distance all nondegenerate for
   0 < l <= 1 (Hessian determinants on the 8 points).
6. B2 and B3 identities: Phi - Psi == e^2 theta G1, Phi + Psi == e^2 G2 (B2) and e^2 eta G1 (B3), with l^2 = e theta.
"""
import sympy as sp

e, ph, qh, b, th, eh, eta, thh, sw, cw, tauh, S, l, w = sp.symbols('e ph qh b th eh eta thh sw cw tauh S l w', real=True)
u = sp.symbols('u', real=True)
cu, su = sp.cos(u), sp.sin(u)


def D(x):
    return sp.diff(x, u)


def trig_zero(expr):
    """exact test that a trig polynomial in u is identically zero (rewrite in exp and expand)"""
    ex = sp.expand(sp.expand_trig(expr).rewrite(sp.exp))
    return sp.simplify(ex) == 0


# ---------- 1. B1 ----------
r = 1 - e * cu
Z = e * (ph * (cu - e) + qh * b * su)
rho2 = r ** 2 - Z ** 2
A = D(r ** 2)
Phi = A ** 2 * rho2
Psi = D(rho2) ** 2
At = 2 * su * r
zh = ph * (cu - e) + qh * b * su
zh2 = zh ** 2
dzh2 = D(zh2)
G1 = At ** 2 * (-2 * cu + e * cu ** 2 - e * zh2) + 2 * At * dzh2 - e * dzh2 ** 2
G2 = At ** 2 * (r ** 2 - e ** 2 * zh2) + (At - e * dzh2) ** 2
print("1. B1: Phi - Psi == e^3 G1:", sp.expand(Phi - Psi - e ** 3 * G1) == 0,
      "; Phi + Psi == e^2 G2:", sp.expand(Phi + Psi - e ** 2 * G2) == 0)

# ---------- 2. b-reduction ----------
for name, G in (("G1", G1), ("G2", G2)):
    P = sp.Poly(sp.expand(G), b)
    R = sp.rem(P, sp.Poly(b ** 2 - (1 - e ** 2), b))
    R0 = R.coeff_monomial(1); R1 = R.coeff_monomial(b)
    print(f"2. {name} mod b^2=1-e^2: b-linear part divisible by ph:", sp.expand(R1.subs(ph, 0)) == 0,
          "; R0 at D1b point (e=0, ph=0, qh=1) == 0:", sp.expand(R0.subs({e: 0, ph: 0, qh: 1})) == 0 if name == "G1" else "n/a")
FB1 = G1 / 2 - tauh * G2
F0 = sp.expand(FB1.subs({e: 0, b: 1}))
F0_claim = 2 * su * ((qh ** 2 - ph ** 2 - 1) * sp.sin(2 * u) + 2 * ph * qh * sp.cos(2 * u) - 4 * tauh * su)
print("3. leading term F0 formula:", trig_zero(F0 - F0_claim))
# vanishing: coefficients of the independent functions sin u sin 2u, sin u cos 2u, sin^2 u
sol = sp.solve([qh ** 2 - ph ** 2 - 1, ph * qh, tauh], [ph, qh, tauh], dict=True)
print("3. F0 == 0 identically iff:", [s_ for s_ in sol if all(v.is_nonnegative is not False for k, v in s_.items() if k in (ph, qh))])

# ---------- 4. CT expansion ----------
# f/e^2 at first order in (e, theta, tau) equals e*F0 with tauh = tau/e, ph = l sin w/e, qh = l cos w/e, theta = l^2/e
thv, tau = sp.symbols('theta tau', real=True)
first = sp.expand(e * F0.subs({tauh: tau / e, ph: sp.sqrt(thv * e) * sp.sin(w) / e, qh: sp.sqrt(thv * e) * sp.cos(w) / e}))
claim = 2 * su * (thv * sp.sin(2 * (u + w)) - 4 * tau * su - e * sp.sin(2 * u))
print("4. CT first-order expansion f/e^2 = 2 sin u [theta sin2(u+w) - 4 tau sin u - e sin 2u]:", trig_zero(first - claim))

# ---------- 5. Lemma E0 ----------
P0 = (l * sp.sin(w)); Q0 = (l * sp.cos(w))
Z0 = P0 * cu + Q0 * su                      # e = 0, b = 1
rho2_0 = 1 - Z0 ** 2
Pol0 = 0 - S * D(rho2_0) ** 2               # A = 0 at e = 0
print("5. E0: Pol == -S l^4 sin^2 2(u+w):", trig_zero(Pol0 + S * l ** 4 * sp.sin(2 * (u + w)) ** 2))
# circle-circle: d^2 = 1 + s^2 - 2 s (X . C), X(u) = (cos u) N + (sin u) M with node line N; f = x.y = cos u cos v + c sin u sin v
v_, c = sp.symbols('v c', real=True)
fxy = sp.cos(u) * sp.cos(v_) + c * sp.sin(u) * sp.sin(v_)
H = sp.hessian(fxy, (u, v_))
pts = [(0, 0), (0, sp.pi), (sp.pi, 0), (sp.pi, sp.pi), (sp.pi / 2, sp.pi / 2), (sp.pi / 2, 3 * sp.pi / 2),
       (3 * sp.pi / 2, sp.pi / 2), (3 * sp.pi / 2, 3 * sp.pi / 2)]
grads_ok = all(sp.simplify(sp.diff(fxy, u).subs({u: a, v_: bb})) == 0 and sp.simplify(sp.diff(fxy, v_).subs({u: a, v_: bb})) == 0
               for a, bb in pts)
dets = [sp.factor(sp.simplify(H.subs({u: a, v_: bb}).det())) for a, bb in pts]
print("5. E0 circle-circle (c = cos I): 8 critical points:", grads_ok, "; Hessian determinants:", sorted(set(dets), key=str),
      "(nonzero for |c| < 1, i.e. 0 < l <= 1)")

# ---------- 6. B2 / B3 ----------
for chart in ("B2", "B3"):
    if chart == "B2":
        ee = th * eh; L2 = ee * th; thv_ = th
    else:
        ee = eta * eh; thv_ = eta * thh; L2 = ee * thv_
    r = 1 - ee * cu
    zt = sw * (cu - ee) + b * cw * su
    rho2 = r ** 2 - L2 * zt ** 2
    A = D(r ** 2)
    Phi = A ** 2 * rho2; Psi = D(rho2) ** 2
    At = 2 * su * r; zt2 = zt ** 2; dzt2 = D(zt2)
    if chart == "B2":
        G1c = At ** 2 * eh * (-2 * cu + ee * cu ** 2 - th * zt2) + 2 * At * dzt2 - th * dzt2 ** 2
        nrm = ee ** 2 * th
    else:
        G1c = At ** 2 * eh * (-2 * cu + ee * cu ** 2 - thv_ * zt2) + 2 * thh * At * dzt2 - eta * thh ** 2 * dzt2 ** 2
        nrm = ee ** 2 * eta
    G2c = At ** 2 * (r ** 2 - ee * thv_ * zt2) + (At - thv_ * dzt2) ** 2
    # use sw^2 + cw^2 = 1 where needed: compare after substituting sw = sin(w), cw = cos(w)
    sub = {sw: sp.sin(w), cw: sp.cos(w)}
    ok1 = trig_zero((Phi - Psi - nrm * G1c).subs(sub))
    ok2 = trig_zero((Phi + Psi - ee ** 2 * G2c).subs(sub))
    print(f"6. {chart}: Phi - Psi == norm G1: {ok1}; Phi + Psi == e^2 G2: {ok2}")

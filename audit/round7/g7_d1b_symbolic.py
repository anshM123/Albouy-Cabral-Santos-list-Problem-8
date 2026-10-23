"""AUDIT-7 g7: the D1b claim of PROOF.md §4 / paper §5: at the D1b point (e, ph, qh, tauh) = (0, 0, 1, 0) of the B1
chart (p = e ph, q = e qh, t = 1/2 + e tauh, chart function F_B1 = f / e^3 = G1/2 - tauh G2), the first-order part is
  2 sin u [2a sin 2u + 2p cos 2u - 4 tau sin u + e r1(u)],   a = qh - 1, p = ph, tau = tauh,
with r1 containing sin 3u, and the four functions are linearly independent (so the leading term never vanishes
identically).  Independent derivation from the a-chart closed forms (b = sqrt(1 - e^2), f = (1-t) Phi - t Psi)."""
import sympy as sp

u, e, ph, qh, tauh, a = sp.symbols('u e ph qh tauh a', real=True)
su, cu = sp.sin(u), sp.cos(u)
b = sp.sqrt(1 - e ** 2)
r = 1 - e * cu
Z = e * (ph * (cu - e) + qh * b * su)
rho2 = r ** 2 - Z ** 2
A = sp.diff(r ** 2, u)
Phi = A ** 2 * rho2; Psi = sp.diff(rho2, u) ** 2
t = sp.Rational(1, 2) + e * tauh
f = (1 - t) * Phi - t * Psi
F = sp.expand(f / e ** 3)                              # chart function of B1 (exact division, checked below)
print('f/e^3 polynomial in e (no negative powers):', all(m.as_powers_dict().get(e, 0) >= 0
                                                        for m in sp.Add.make_args(sp.expand(sp.series(F, e, 0, 2).removeO()))))
pt = {e: 0, ph: 0, qh: 1, tauh: 0}
F_at = sp.simplify(sp.limit(F.subs({ph: 0, qh: 1, tauh: 0}), e, 0))
print('F at the D1b point:', F_at)
lin = {}
for v in (qh, ph, tauh):
    lin[v] = sp.simplify(sp.limit(sp.diff(F, v).subs({ph: 0, qh: 1, tauh: 0}), e, 0))
# derivative in e at the point: F(e) with ph = 0, qh = 1, tauh = 0, expanded in e
Fe = sp.series(F.subs({ph: 0, qh: 1, tauh: 0}), e, 0, 2).removeO()
de = sp.simplify(Fe.coeff(e, 1))
print('dF/da   =', sp.simplify(sp.expand_trig(lin[qh] - 2 * su * 2 * sp.sin(2 * u))) == 0, '(claim 2 sin u * 2 sin 2u)')
print('dF/dp   =', sp.simplify(sp.expand_trig(lin[ph] - 2 * su * 2 * sp.cos(2 * u))) == 0, '(claim 2 sin u * 2 cos 2u)')
print('dF/dtau =', sp.simplify(sp.expand_trig(lin[tauh] + 2 * su * 4 * su)) == 0, '(claim -2 sin u * 4 sin u)')
r1 = sp.simplify(de / (2 * su))
r1t = sp.fu(sp.expand_trig(r1)) if r1 != 0 else 0
print('dF/de / (2 sin u) = r1(u) =', sp.simplify(sp.expand(sp.expand_trig(r1).rewrite(sp.exp))).rewrite(sp.cos))
# Fourier content of r1 and independence of the four functions
xs = [sp.Rational(k, 7) + sp.Rational(1, 11) for k in range(1, 9)]
funcs = [2 * su * 2 * sp.sin(2 * u), 2 * su * 2 * sp.cos(2 * u), -8 * su ** 2, de]
M = sp.Matrix([[sp.N(fn.subs(u, x), 40) for fn in funcs] for x in xs])
print('rank of the four first-order functions at 8 points:', M.rank(iszerofunc=lambda z: abs(z) < 1e-30))
c3 = sp.integrate(sp.expand_trig(r1) * sp.sin(3 * u), (u, 0, 2 * sp.pi)) / sp.pi
print('sin 3u coefficient of r1:', sp.simplify(c3))
coef = {}
for k in range(0, 7):
    ck = sp.simplify(sp.integrate(sp.expand_trig(r1) * sp.cos(k * u), (u, 0, 2 * sp.pi)) / (sp.pi * (2 if k == 0 else 1)))
    sk = sp.simplify(sp.integrate(sp.expand_trig(r1) * sp.sin(k * u), (u, 0, 2 * sp.pi)) / sp.pi)
    if ck != 0: coef['cos%du' % k] = ck
    if sk != 0: coef['sin%du' % k] = sk
print('r1(u) Fourier coefficients:', coef)

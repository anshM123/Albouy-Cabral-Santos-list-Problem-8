"""AUDIT a9: symbolic checks of CORNER_E1.md, the axis-aligned factorisation, Pol != 0, and the symmetry reductions.

1. f-chart normalisation of CORNER_E1 (R = 1, mu = q/R): F_-+ = 2 mu e (1+e) sin f sigma -+ D B is a positive
   multiple of g_-+' = (r^2)' -+ 2 rho' (factor D^3 sigma / (mu (1+e))), and F_- F_+ = mu^2 (PhiF - kappa^2 PsiF).
2. B(pi) = -2 P_z Q_z (1-e),  B'(pi) = 2 P_z^2 - 2 e + 2 (e-1) Q_z^2,  D(pi) = 1 - e.
3. Axis-aligned case P_z = 1, Q_z = 0:  B = 2 sin f (cos f + e)  (= 2 sin f D only at e = 1);
   PhiF - kappa^2 PsiF = 4 sin^2 f H1 H2,  H1 = (1+e) e sin f - kappa D (cos f + e),  H2(f) = -H1(-f);
   P_f = (1-t)(PhiF - kappa^2 PsiF).  At e = 1, mu = 0 (kappa = inf): P ~ -4 sin^2 f D^4.
4. AM-GM threshold in Lemma N2: 2|P_z Q_z| eps |phi| <= 2 c^2 eps + 3 c^2 phi^2 whenever eps <= 24 c^4.
5. Pol is not identically 0 for e in (0,1): c_0 = 0 forces q = 0, p = e, and then c_1 = e^5 (1-e^2)/4 != 0.
6. Evolute: for the ellipse with semi-axes a (along the apse line) and b, the evolute's cusps on the apse line are
   at distance a e^2 from the centre; the focus is at distance a e > a e^2.
"""
import sympy as sp

e, p, q, f, mu, kap, S, phi, eps, c2 = sp.symbols('e p q f mu kappa S phi epsilon c2', real=True)
zeta = p * sp.cos(f) + q * sp.sin(f)
sig2 = 1 - zeta ** 2
sig = sp.sqrt(sig2)
D = 1 + e * sp.cos(f)
B = sp.diff(sig2, f) + 2 * e * (sp.sin(f) - q * zeta)

# 1.
r = mu * (1 + e) / D
rho = r * sig
for sgn in (-1, 1):
    gprime = sp.diff(r ** 2, f) + sgn * 2 * sp.diff(rho, f)      # sgn=-1: g_-' ; +1: g_+'
    F = 2 * mu * e * (1 + e) * sp.sin(f) * sig + sgn * D * B
    ratio = sp.simplify(gprime * D ** 3 * sig / (mu * (1 + e)) - F)
    print(f"1. F_{'-' if sgn < 0 else '+'} == g' * D^3 sigma/(mu(1+e)):", sp.simplify(sp.expand(ratio)) == 0)
Fm = 2 * mu * e * (1 + e) * sp.sin(f) * sig - D * B
Fp = 2 * mu * e * (1 + e) * sp.sin(f) * sig + D * B
PhiF = 4 * (1 + e) ** 2 * e ** 2 * sp.sin(f) ** 2 * sig2
PsiF = D ** 2 * B ** 2
print("1. F_- F_+ == mu^2 (PhiF - kappa^2 PsiF) with kappa = 1/mu:",
      sp.simplify(sp.expand(Fm * Fp - mu ** 2 * (PhiF - PsiF / mu ** 2))) == 0)
# 2.
print("2. B(pi) == -2 p q (1-e):", sp.simplify(B.subs(f, sp.pi) + 2 * p * q * (1 - e)) == 0)
Bp = sp.diff(B, f)
print("2. B'(pi) == 2p^2 - 2e + 2(e-1) q^2:", sp.simplify(Bp.subs(f, sp.pi) - (2 * p ** 2 - 2 * e + 2 * (e - 1) * q ** 2)) == 0)
print("2. D(pi) == 1-e:", sp.simplify(D.subs(f, sp.pi) - (1 - e)) == 0)
# 3.
Bax = sp.simplify(B.subs({p: 1, q: 0}))
print("3. B at P_z=1,Q_z=0:", sp.factor(sp.expand_trig(Bax)), "; equals 2 sin f D ?", sp.simplify(Bax - 2 * sp.sin(f) * D) == 0,
      "; at e=1:", sp.simplify((Bax - 2 * sp.sin(f) * D).subs(e, 1)) == 0)
H1 = (1 + e) * e * sp.sin(f) - kap * D * (sp.cos(f) + e)
H2 = (1 + e) * e * sp.sin(f) + kap * D * (sp.cos(f) + e)
lhs = (PhiF - kap ** 2 * PsiF).subs({p: 1, q: 0})
print("3. PhiF - kappa^2 PsiF == 4 sin^2 f H1 H2 at P_z=1,Q_z=0:", sp.simplify(sp.expand(lhs - 4 * sp.sin(f) ** 2 * H1 * H2)) == 0)
print("3. H2(f) == -H1(-f):", sp.simplify(H2 + H1.subs(f, -f)) == 0)
t = sp.symbols('t', positive=True)
Pf = (1 - t) * PhiF - t * PsiF
print("3. P_f == (1-t)(PhiF - kappa^2 PsiF) when kappa^2 = t/(1-t):",
      sp.simplify(sp.expand(Pf - (1 - t) * (PhiF - t / (1 - t) * PsiF))) == 0)
Pm0 = (Fm * Fp).subs({p: 1, q: 0, e: 1, mu: 0})
print("3. at e=1, mu=0, P_z=1: F_-F_+ == -4 sin^2 f D^4:", sp.simplify(Pm0 + 4 * sp.sin(f) ** 2 * (1 + sp.cos(f)) ** 4) == 0)
# 4. AM-GM: minimum over phi of 2 c2 eps + 3 c2 phi^2 - 2 k eps |phi| with k = |P_z Q_z| <= 1/2
k = sp.symbols('k', positive=True)
x = sp.symbols('x', positive=True)
expr = 2 * c2 * eps + 3 * c2 * x ** 2 - 2 * k * eps * x
xm = sp.solve(sp.diff(expr, x), x)[0]
mn = sp.simplify(expr.subs(x, xm))
print("4. min over |phi| of (2c^2 eps + 3c^2 phi^2 - 2k eps |phi|) =", mn, "; >= 0 iff eps <= 6 c^4/k^2 (=24 c^4 at k=1/2):",
      sp.simplify(sp.solve(sp.Eq(mn, 0), eps)[0]))
# 5. Pol != 0
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'code'))
import coeffs_ce
bb = sp.symbols('b', positive=True)
re_, im_ = coeffs_ce.coeffs(e, bb, p, q, S)
c0 = sp.nsimplify(re_[0]) + sp.I * sp.nsimplify(im_[0])
print("5. c_0 =", sp.factor(sp.nsimplify(re_[0])), "+ i", sp.factor(sp.nsimplify(im_[0])))
c1q0 = sp.factor(sp.nsimplify(re_[1]).subs({q: 0, p: e}))
print("5. Re c_1 at q = 0, p = e:", c1q0, "; Im c_1 there:", sp.nsimplify(im_[1]).subs({q: 0, p: e}))
# 6. evolute cusp: centre of curvature at the apse (a, 0): x = a - b^2/a = a e^2 (with b^2 = a^2 (1-e^2))
a = sp.symbols('a', positive=True)
bsq = a ** 2 * (1 - e ** 2)
print("6. evolute cusp distance from centre (apse line) =", sp.simplify(a - bsq / a), "; focus distance a e;",
      "difference a e - a e^2 =", sp.factor(a * e - (a - bsq / a)))

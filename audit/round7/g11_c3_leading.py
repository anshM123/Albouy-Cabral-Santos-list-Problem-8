"""AUDIT-7 g11: the C3 leading form printed in C3_CORE.md §2 and paper §5, from MY OWN physics:
Q(tau) = P(f) (1 + tau^2)^6, f = pi + 2 atan(tau), P = mu^2 k^2 sin^2 f sigma^2 - D^2 B^2; substitute
tau = nu y, Qz = nu a, Pz = 1 - nu^2 K, e = 1 - nu^2 E, mu = nu^3 m.  Claims:
  Q = nu^10 [P0(2y) + nu^2 P1 + ...] (only even powers of nu),
  P0(x) = 16 m^2 x^2 (lam + (x - a)^2) - (E + x^2/2)^2 Bh(x)^2,  lam = 2K - a^2,
  Bh(x) = 2E(x - a) - 2 lam x - 2 a^2 x + 3 a x^2 - x^3 = 2E(x - a) - 4Kx + 3ax^2 - x^3   (paper form);
  at E = m = 0: P0 = -(x^6/4)(x^2 - 3ax + 4K)^2;
  level 2 (weights (1,2,2) for (y, E, m)): leading part y^2 [4m^2 - 2K(E + 2y^2)^2] up to a positive factor."""
import sympy as sp

tau, nu, y, a, K, E, m, x = sp.symbols('tau nu y a K E m x', real=True)
cf = -(1 - tau ** 2) / (1 + tau ** 2); sf = -2 * tau / (1 + tau ** 2)
dcf, dsf = -sf, cf                                   # d/df of cos f and sin f
Pz, Qz, e, mu = 1 - nu ** 2 * K, nu * a, 1 - nu ** 2 * E, nu ** 3 * m
zeta = Pz * cf + Qz * sf
dzeta = Pz * dcf + Qz * dsf
B = -2 * zeta * dzeta + 2 * e * (sf - Qz * zeta)
D = 1 + e * cf
k = 2 * e * (1 + e)
P = mu ** 2 * k ** 2 * sf ** 2 * (1 - zeta ** 2) - D ** 2 * B ** 2
Q = sp.together(P * (1 + tau ** 2) ** 6)
Q = sp.expand(sp.cancel(Q))
print('Q is a polynomial in tau:', sp.Poly(Q, tau).degree(), '(degree in tau)')
Qn = sp.expand(Q.subs(tau, nu * y))
poly = sp.Poly(Qn, nu)
degs = sorted(set(mon[0] for mon in poly.monoms()))
print('powers of nu present:', degs[:6], '...', degs[-3:], '; all even:', all(d % 2 == 0 for d in degs), '; lowest:', degs[0])
lead = sp.expand(poly.coeff_monomial(nu ** degs[0]))
lam = 2 * K - a ** 2
Bh = 2 * E * (x - a) - 4 * K * x + 3 * a * x ** 2 - x ** 3
P0 = 16 * m ** 2 * x ** 2 * (lam + (x - a) ** 2) - (E + x ** 2 / 2) ** 2 * Bh ** 2
print('nu^10 coefficient == P0(2y):', sp.expand(lead - P0.subs(x, 2 * y)) == 0)
print('Bh forms agree:', sp.expand(Bh - (2 * E * (x - a) - 2 * lam * x - 2 * a ** 2 * x + 3 * a * x ** 2 - x ** 3)) == 0)
print('E = m = 0 form:', sp.expand(P0.subs({E: 0, m: 0}) + x ** 6 / 4 * (x ** 2 - 3 * a * x + 4 * K) ** 2) == 0)
eta = sp.symbols('eta', positive=True)
L2 = sp.expand(lead.subs({y: eta * y, E: eta ** 2 * E, m: eta ** 2 * m}))
p2 = sp.Poly(L2, eta)
low = min(mn[0] for mn in p2.monoms())
lead2 = sp.factor(p2.coeff_monomial(eta ** low))
print('level-2 weighted order:', low, '; leading part:', lead2,
      '; == 32 K y^2 [4m^2 - 2K(E + 2y^2)^2]:', sp.expand(lead2 - 32 * K * y ** 2 * (4 * m ** 2 - 2 * K * (E + 2 * y ** 2) ** 2)) == 0)

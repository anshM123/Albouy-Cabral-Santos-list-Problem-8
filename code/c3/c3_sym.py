"""C3 core: tan-half-angle polynomial of the R-normalised corner function and its weighted expansion.
P(f) = mu^2 k^2 sin^2 f sigma^2 - D^2 B^2  (circle radius 1, pericentre distance mu, k = 2e(1+e)),
zeta = Pz cos f + Qz sin f, sigma^2 = 1 - zeta^2, B = (sigma^2)' + 2e (sin f - Qz zeta), D = 1 + e cos f.
f = pi + phi, th = tan(phi/2):  Q(th) = P (1 + th^2)^6  is a polynomial of degree <= 12 in th."""
import sympy as sp
th, Pz, Qz, e, mu = sp.symbols('th Pz Qz e mu', real=True)
cf, sf = sp.symbols('cf sf', real=True)
zeta = Pz * cf + Qz * sf
dzeta = -Pz * sf + Qz * cf
s2 = 1 - zeta**2
B = -2 * zeta * dzeta + 2 * e * (sf - Qz * zeta)
D = 1 + e * cf
k = 2 * e * (1 + e)
P = sp.expand(mu**2 * k**2 * sf**2 * s2 - D**2 * B**2)
den = 1 + th**2
Q = sp.expand(sp.cancel(P.subs({cf: -(1 - th**2) / den, sf: -2 * th / den}) * den**6))
QP = sp.Poly(Q, th)
print('degree in th:', QP.degree())
# weighted substitution
nu, y, a, K, E, m = sp.symbols('nu y a K E m', real=True)
W = sp.expand(Q.subs({th: nu * y, Qz: nu * a, Pz: 1 - nu**2 * K, e: 1 - nu**2 * E, mu: nu**3 * m}))
WP = sp.Poly(W, nu)
degs = sorted(set(mon[0] for mon in WP.monoms()))
print('nu-degrees present:', degs[:6], '... max', degs[-1])
lead = sp.expand(WP.coeff_monomial(nu**degs[0]))
x = sp.symbols('x', real=True)
# expected leading polynomial (C3_CORE.md §2) in x = phi/delta with lam = 2K - a^2, psi = a
lam = 2 * K - a**2; psi = a
Bh = 2 * E * (x - psi) - 2 * lam * x - 2 * psi**2 * x + 3 * psi * x**2 - x**3
P0 = 16 * m**2 * x**2 * (lam + (x - psi)**2) - (E + x**2 / 2)**2 * Bh**2
cmp = sp.expand(lead - sp.expand(P0.subs(x, 2 * y)) * sp.Rational(1, 1))
ratio = sp.simplify(lead / sp.expand(P0.subs(x, 2 * y)))
print('leading / P0(2y) =', ratio)

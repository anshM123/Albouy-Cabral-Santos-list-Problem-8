"""AUDIT a3b: (i) the ambiguous W0/AX cases are exact double roots at u = 0, pi (50-digit check);
(ii) families near the axis (l -> 1, rho small but nonzero) and near the maximum N = 10."""
import numpy as np, mpmath as mp
from a3_reduction import compare, crit_2d, pol_real_roots, run_family
mp.mp.dps = 50
rng = np.random.default_rng(4242)

def pol_and_deriv(e, p, q, S, u0):
    e, p, q, S, u0 = [mp.mpf(x) for x in (e, p, q, S, u0)]
    b = mp.sqrt(1 - e*e)
    def pol(u):
        r = 1 - e*mp.cos(u); Z = p*(mp.cos(u)-e) + b*q*mp.sin(u); Zp = -p*mp.sin(u) + b*q*mp.cos(u)
        A = 2*e*mp.sin(u)*r; rho2 = r*r - Z*Z; d = A - 2*Z*Zp
        return A*A*rho2 - S*d*d
    return pol(u0), mp.diff(pol, u0), mp.diff(pol, u0, 2)

# (i) W0 and AX: exact double root at u = 0 and u = pi
worst = 0
for it in range(40):
    e = rng.uniform(0.05, 0.95); I = rng.uniform(0, np.pi); w = rng.choice([0.0, np.pi/2]); s = np.exp(rng.uniform(np.log(.1), np.log(10)))
    if it % 2: I, w = np.pi/2, np.pi/2   # AX
    p, q = np.sin(w)*np.sin(I), np.cos(w)*np.sin(I)
    for u0 in (0.0, np.pi):
        v0, v1, v2 = pol_and_deriv(e, p, q, s*s, u0)
        worst = max(worst, float(abs(v0)), float(abs(v1)))
print("(i) max |Pol|, |Pol'| at u = 0, pi over 40 W0/AX configs (50 digits):", worst)

# (ii) near-axis family: l = 1 - delta
lu = lambda a, b: np.exp(rng.uniform(np.log(a), np.log(b)))
run_family("L1-", lambda: (rng.uniform(0.05, 0.95), np.pi/2 - lu(1e-4, 0.1), rng.uniform(0, 2*np.pi), lu(0.1, 10)), 100)

# (iii) near-maximum family: random (e, I, w), choose s maximising the Pol count on a log grid
def best_s(e, I, w):
    p, q = np.sin(w)*np.sin(I), np.cos(w)*np.sin(I)
    best = (-1, None)
    for s in np.exp(np.linspace(np.log(0.05), np.log(20), 160)):
        roots, amb = pol_real_roots(e, p, q, s*s)
        if amb == 0 and len(roots) > best[0]:
            best = (len(roots), s)
    return best
T = []
while len(T) < 40:
    e = rng.uniform(0.3, 0.97); I = rng.uniform(0.4, np.pi - 0.4); w = rng.uniform(0, 2*np.pi)
    n, s = best_s(e, I, w)
    if n == 10:
        T.append((e, I, w, s))
it = iter(T)
run_family("T(N=10)", lambda: next(it), len(T))

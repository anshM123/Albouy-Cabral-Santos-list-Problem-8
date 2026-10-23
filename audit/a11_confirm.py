"""AUDIT a11: high-precision confirmation of the two counterexamples to LOG/STATEMENT numerical claims.
(i)  delta = 0.554 < 0.836 at e = 0.95: 2-D Newton count = 10 and the off-circle pair from 60-digit polyroots.
(ii) parabola-circle relative gap -> sqrt(2) < sqrt(3): exact factorisation at l = 1, w = 0, and a nearby
     interior point (l < 1, w > 0) with 8 real roots and gap < 1.5."""
import numpy as np, mpmath as mp
from a3_reduction import crit_2d
mp.mp.dps = 60
def Q_coeffs(e, l, w, S):
    e, l, w, S = [mp.mpf(x) for x in (e, l, w, S)]
    p, q, b = l*mp.sin(w), l*mp.cos(w), mp.sqrt(1-e*e)
    N = 64
    vals = []
    for n in range(N):
        u = 2*mp.pi*n/N
        r = 1 - e*mp.cos(u); Z = p*(mp.cos(u)-e) + b*q*mp.sin(u); Zp = -p*mp.sin(u) + b*q*mp.cos(u)
        A = 2*e*mp.sin(u)*r; rho2 = r*r - Z*Z; d = A - 2*Z*Zp
        vals.append(A*A*rho2 - S*d*d)
    c = [mp.fsum(vals[n]*mp.exp(-1j*k*2*mp.pi*n/N) for n in range(N))/N for k in range(-6, 7)]
    return c
for (e, l, w, t) in [(0.95000, 0.999999, 0.001157, 0.008854), (0.94998, 0.999996, 0.000364, 0.009032)]:
    S = t/(1-t)
    c = Q_coeffs(e, l, w, S)
    rts = mp.polyroots(c[::-1], maxsteps=400, extraprec=200)
    lz = sorted([abs(mp.log(abs(z))) for z in rts])
    n_on = sum(1 for x in lz if x < mp.mpf('1e-20'))
    cps = crit_2d(e, np.arcsin(l), w, np.sqrt(S), n=72)
    print(f"(i) e={e} l={l} w={w} t={t}: roots on circle {n_on}, smallest off-circle |log|z|| = {mp.nstr(lz[n_on], 6)}, 2-D count {len(cps)}")
# (ii)
from numpy.polynomial import polynomial as Pn
def parab(k, l, w):
    pz, qz = l*np.sin(w), l*np.cos(w)
    X2 = np.array([1.0, 0, 2.0, 0, 1.0]); Z = np.array([pz, 2*qz, -pz])
    rho2 = Pn.polysub(X2, Pn.polymul(Z, Z)); A = Pn.polyder(X2); dr = Pn.polyder(rho2)
    return Pn.polysub(Pn.polymul(Pn.polymul(A, A), rho2), k*k*Pn.polymul(dr, dr))
for (k, l, w) in [(1.0, 1.0, 0.0), (1.05, 1.0, 0.0), (1.05, 0.999, 0.01), (1.1, 0.995, 0.02), (2.0, 1.0, 0.0)]:
    c = [mp.mpf(x) for x in parab(k, l, w)]
    while c[-1] == 0: c.pop()
    try:
        rts = mp.polyroots(c[::-1], maxsteps=400, extraprec=300)
    except mp.libmp.libhyper.NoConvergence:
        print(f"(ii) k={k} l={l} w={w}: polyroots no convergence (exact multiple roots); exact factorisation 16 w^2 (1-w^2)^2 ((1+w^2)^2 - k^2) gives gap sqrt(1+k) = {np.sqrt(1+k):.5f}")
        continue
    rel = sorted([(abs(mp.im(z))/(1+abs(mp.re(z))), z) for z in rts], key=lambda a: a[0])
    nreal = sum(1 for g, z in rel if g < mp.mpf('1e-25'))
    cpx = [g for g, z in rel if g >= mp.mpf('1e-25')]
    print(f"(ii) k={k} l={l} w={w}: real roots {nreal}, min relative gap of complex pair {mp.nstr(min(cpx), 6) if cpx else None}")

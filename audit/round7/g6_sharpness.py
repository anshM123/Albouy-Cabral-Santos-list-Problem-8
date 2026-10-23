"""AUDIT-7 g6: the sharpness example of PROOF.md / paper Theorem 1: e = 0.673, sin I = 0.992, omega = 0.338,
R/a = 0.7772 (I taken in (0, pi/2]; the count does not depend on the choice I or pi - I by symmetry).
Three independent counts:
  (1) Pol(u) = A^2 rho^2 - s^2 ((rho^2)')^2 from MY closed forms at 80 digits: zeros on |z| = 1 (with the root
      separation and the distance of the remaining complex zeros from the unit circle);
  (2) sign changes of the branch functions h_-+ = A rho -+ s (rho^2)' (= rho g_-+') on a fine grid, refined;
  (3) 2-D damped Newton on grad d^2(u, v) from a 56 x 56 grid on the torus (audit a3 crit_2d, 3-D vectors only),
      plus the Hessian signature of each critical point (minima / saddles / maxima), and a 2-D residual at 80 digits.
"""
import os, sys
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
sys.dont_write_bytecode = True
from a3_reduction import crit_2d
mp.mp.dps = 80
e, l, w, s = mp.mpf('0.673'), mp.mpf('0.992'), mp.mpf('0.338'), mp.mpf('0.7772')
I = mp.asin(l)
b = mp.sqrt(1 - e * e); p, q = l * mp.sin(w), l * mp.cos(w)


def parts(u):
    r = 1 - e * mp.cos(u)
    Z = p * (mp.cos(u) - e) + q * b * mp.sin(u)
    Zp = -p * mp.sin(u) + q * b * mp.cos(u)
    A = 2 * e * mp.sin(u) * r
    rho2 = r * r - Z * Z
    d = A - 2 * Z * Zp
    return A, rho2, d


def pol(u):
    A, rho2, d = parts(u)
    return A * A * rho2 - s * s * d * d


N = 64
vals = [pol(2 * mp.pi * j / N) for j in range(N)]
c = [mp.fsum(vals[j] * mp.expj(-2 * mp.pi * j * k / N) for j in range(N)) / N for k in range(-6, 7)]
z = mp.polyroots(c[::-1], maxsteps=800, extraprec=800)
on = sorted([z0 for z0 in z if abs(mp.log(abs(z0))) < mp.mpf(10) ** -40], key=lambda t: float(mp.arg(t) % (2 * mp.pi)))
off = [z0 for z0 in z if abs(mp.log(abs(z0))) >= mp.mpf(10) ** -40]
us = [mp.arg(z0) % (2 * mp.pi) for z0 in on]
sep = min(abs(mp.arg(on[i] / on[(i + 1) % len(on)])) for i in range(len(on)))
print('(1) Pol: zeros on the unit circle', len(on), '; min angular separation', mp.nstr(sep, 5),
      '; other zeros |log|z|| =', [mp.nstr(abs(mp.log(abs(z0))), 5) for z0 in off])
print('    real roots u =', [mp.nstr(u, 8) for u in us])

# (2) branch sign changes
def h(u, sg):
    A, rho2, d = parts(u)
    return A * mp.sqrt(rho2) - sg * s * d       # sg = +1: h_- = A rho - s (rho^2)'; sg = -1: h_+


for sg, name in ((1, 'g_-'), (-1, 'g_+')):
    M = 4000
    xs = [2 * mp.pi * (j + mp.mpf(1) / 7) / M for j in range(M)]
    hv = [h(x, sg) for x in xs]
    ch = sum(1 for j in range(M) if mp.sign(hv[j]) != mp.sign(hv[(j + 1) % M]))
    print(f'(2) sign changes of {name}\' : {ch}')

# (3) 2-D Newton on the torus (float, from audit a3) and Hessian signature
cps = crit_2d(float(e), float(I), float(w), float(s))
print('(3) 2-D Newton critical points:', len(cps))
P = [mp.cos(w), mp.sin(w) * mp.cos(I), mp.sin(w) * mp.sin(I)]
Q = [-mp.sin(w), mp.cos(w) * mp.cos(I), mp.cos(w) * mp.sin(I)]
kinds = {}
maxres = 0
for (u0, v0) in cps:
    u, v = mp.mpf(u0), mp.mpf(v0)
    for _ in range(40):                     # Newton refinement at 80 digits
        X = [(mp.cos(u) - e) * P[i] + b * mp.sin(u) * Q[i] for i in range(3)]
        Xp = [-mp.sin(u) * P[i] + b * mp.cos(u) * Q[i] for i in range(3)]
        Xpp = [-mp.cos(u) * P[i] - b * mp.sin(u) * Q[i] for i in range(3)]
        C = [s * mp.cos(v), s * mp.sin(v), 0]; Cp = [-s * mp.sin(v), s * mp.cos(v), 0]; Cpp = [-C[0], -C[1], 0]
        D = [X[i] - C[i] for i in range(3)]
        gu = 2 * sum(D[i] * Xp[i] for i in range(3)); gv = -2 * sum(D[i] * Cp[i] for i in range(3))
        huu = 2 * (sum(Xp[i] ** 2 for i in range(3)) + sum(D[i] * Xpp[i] for i in range(3)))
        hvv = 2 * (sum(Cp[i] ** 2 for i in range(3)) - sum(D[i] * Cpp[i] for i in range(3)))
        huv = -2 * sum(Xp[i] * Cp[i] for i in range(3))
        det = huu * hvv - huv * huv
        du = (hvv * gu - huv * gv) / det; dv = (huu * gv - huv * gu) / det
        u -= du; v -= dv
    res = abs(gu) + abs(gv)
    maxres = max(maxres, float(res))
    k = 'saddle' if det < 0 else ('min' if huu > 0 else 'max')
    kinds[k] = kinds.get(k, 0) + 1
    dmin = min(abs(mp.arg(mp.expj(u - uu))) for uu in us)
    print(f'    u = {mp.nstr(u % (2 * mp.pi), 10)}, v = {mp.nstr(v % (2 * mp.pi), 10)}, {k}, |grad| = {mp.nstr(res, 3)},'
          f' distance to the nearest Pol root {mp.nstr(dmin, 3)}')
print('    kinds', kinds, '; max |grad| after refinement', maxres,
      '; Euler characteristic check (min - saddle + max = 0 on the torus):', kinds.get('min', 0) - kinds.get('saddle', 0) + kinds.get('max', 0))

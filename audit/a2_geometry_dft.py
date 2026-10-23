"""AUDIT a2: check the generated coefficients against the raw 3-D geometry (no formula for Z or rho^2 reused).

For random (e, I, omega) the ellipse X(u) = (cos u - e) P + beta sin u Q is built from the rotation vectors,
rho^2 = X_x^2 + X_y^2 and r^2 = |X|^2 are computed from the vector, derivatives are taken with mpmath.diff
(numerical differentiation at 40 digits), and Phi = ((r^2)')^2 rho^2, Psi = ((rho^2)')^2 are sampled at
N = 32 equally spaced u and transformed by a DFT.  Compared with coeffs3_ce.allc(e, beta, P_z, Q_z).

f-chart: X(f) = r (cos f P + sin f Q), r = (1+e)/(1+e cos f) (pericentre distance 1), derivatives in f, and
  D^8/(1+e)^4 ((r^2)')^2 rho^2  = PhiF,     D^8/(1+e)^4 ((rho^2)')^2 = PsiF   (claimed)
compared with coeffsF.allcF(e, P_z, Q_z).
Also checks the a-chart / f-chart root correspondence: real roots of Pol(u) mapped by the anomaly map
tan(f/2) = sqrt((1+e)/(1-e)) tan(u/2) coincide with the real roots of P_f at kappa = s/(1-e).
"""
import os
import sys
import numpy as np
import mpmath as mp

CODE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'code')
sys.path.insert(0, CODE)
from coeffs3_ce import allc
from coeffsF import allcF

mp.mp.dps = 40
rng = np.random.default_rng(20261022)


def vecs(I, w):
    P = [mp.cos(w), mp.sin(w) * mp.cos(I), mp.sin(w) * mp.sin(I)]
    Q = [-mp.sin(w), mp.cos(w) * mp.cos(I), mp.cos(w) * mp.sin(I)]
    return P, Q


def dft(vals):
    N = len(vals)
    a = []; s = []
    for k in range(7):
        ck = mp.fsum(vals[n] * mp.cos(2 * mp.pi * k * n / N) for n in range(N)) / N
        sk = mp.fsum(vals[n] * mp.sin(2 * mp.pi * k * n / N) for n in range(N)) / N
        a.append(ck if k == 0 else 2 * ck)
        s.append(mp.mpf(0) if k == 0 else 2 * sk)
    return a, s


def check_a(e, I, w):
    e = mp.mpf(e); I = mp.mpf(I); w = mp.mpf(w)
    beta = mp.sqrt(1 - e * e)
    P, Q = vecs(I, w)
    X = lambda u: [(mp.cos(u) - e) * P[i] + beta * mp.sin(u) * Q[i] for i in range(3)]
    r2 = lambda u: sum(x * x for x in X(u))
    rho2 = lambda u: X(u)[0] ** 2 + X(u)[1] ** 2
    N = 32
    us = [2 * mp.pi * n / N for n in range(N)]
    phi = []; psi = []
    for u in us:
        A = mp.diff(r2, u); dr = mp.diff(rho2, u)
        phi.append(A * A * rho2(u)); psi.append(dr * dr)
        # side identity: |X| = 1 - e cos u
        assert abs(mp.sqrt(r2(u)) - (1 - e * mp.cos(u))) < mp.mpf(10) ** -30
    want = [dft(phi), dft(psi)]
    V, _ = allc(float(e), float(beta), float(P[2]), float(Q[2]))
    err = 0.0; mag = 0.0
    for f in range(2):
        for part in range(2):
            for k in range(7):
                err = max(err, abs(float(want[f][part][k]) - float(V[f][part][k])))
                mag = max(mag, abs(float(want[f][part][k])))
    return err, mag


def check_f(e, I, w):
    e = mp.mpf(e); I = mp.mpf(I); w = mp.mpf(w)
    P, Q = vecs(I, w)
    rr = lambda f: (1 + e) / (1 + e * mp.cos(f))
    X = lambda f: [rr(f) * (mp.cos(f) * P[i] + mp.sin(f) * Q[i]) for i in range(3)]
    r2 = lambda f: sum(x * x for x in X(f))
    rho2 = lambda f: X(f)[0] ** 2 + X(f)[1] ** 2
    N = 32
    fs = [2 * mp.pi * n / N for n in range(N)]
    phi = []; psi = []
    for f in fs:
        D = 1 + e * mp.cos(f)
        fac = D ** 8 / (1 + e) ** 4
        A = mp.diff(r2, f); dr = mp.diff(rho2, f)
        phi.append(fac * A * A * rho2(f)); psi.append(fac * dr * dr)
    want = [dft(phi), dft(psi)]
    V, _ = allcF(float(e), float(P[2]), float(Q[2]))
    err = 0.0; mag = 0.0
    for f in range(2):
        for part in range(2):
            for k in range(7):
                err = max(err, abs(float(want[f][part][k]) - float(V[f][part][k])))
                mag = max(mag, abs(float(want[f][part][k])))
    return err, mag


print("== a-chart: generated coefficients vs DFT of the raw geometry ==")
worst = 0.0
for it in range(12):
    e, I, w = rng.uniform(0.01, 0.99), rng.uniform(0, np.pi), rng.uniform(0, 2 * np.pi)
    err, mag = check_a(e, I, w)
    worst = max(worst, err / mag)
    print(f"  e={e:.4f} I={I:.4f} w={w:.4f}: max abs err {err:.2e} (max |coef| {mag:.2e})")
print("  worst relative error:", f"{worst:.2e}")

print("== f-chart: generated coefficients vs DFT of D^8/(1+e)^4 x (geometry) ==")
worst = 0.0
for it in range(12):
    e, I, w = rng.uniform(0.01, 0.99), rng.uniform(0, np.pi), rng.uniform(0, 2 * np.pi)
    err, mag = check_f(e, I, w)
    worst = max(worst, err / mag)
    print(f"  e={e:.4f} I={I:.4f} w={w:.4f}: max abs err {err:.2e} (max |coef| {mag:.2e})")
print("  worst relative error:", f"{worst:.2e}")


# ---------- root correspondence a-chart <-> f-chart ----------
def trig_vals(a, s, x):
    return a[0] + sum(a[k] * np.cos(k * x) + s[k] * np.sin(k * x) for k in range(1, 7))


def roots_on_circle(a, s, tol=1e-7):
    """roots of the trig polynomial on the circle via the degree-12 companion polynomial"""
    c = np.zeros(13, complex)
    c[6] = a[0]
    for k in range(1, 7):
        c[6 + k] = (a[k] - 1j * s[k]) / 2
        c[6 - k] = (a[k] + 1j * s[k]) / 2
    top = 12
    while abs(c[top]) < 1e-14 * np.abs(c).max():
        top -= 1
    low = 0
    while abs(c[low]) < 1e-14 * np.abs(c).max():
        low += 1
    z = np.roots(c[low:top + 1][::-1])
    on = np.abs(np.log(np.abs(z))) < tol
    return np.sort(np.mod(np.angle(z[on]), 2 * np.pi))


print("== root correspondence: Pol(u) (a-chart) vs P_f (f-chart, kappa = s/(1-e)) ==")
agree = 0; tot = 0; bad = []
for it in range(400):
    e = rng.uniform(0.02, 0.98); l = rng.uniform(0, 1); w = rng.uniform(0, np.pi / 2)
    s = np.exp(rng.uniform(np.log(0.05), np.log(20)))
    p, q = l * np.sin(w), l * np.cos(w)
    V, _ = allc(e, np.sqrt(1 - e * e), p, q)
    S = s * s
    a = [V[0][0][k] - S * V[1][0][k] for k in range(7)]; sn = [V[0][1][k] - S * V[1][1][k] for k in range(7)]
    ru = roots_on_circle(a, sn)
    kap = s / (1 - e)
    VF, _ = allcF(e, p, q)
    K2 = kap * kap
    aF = [VF[0][0][k] - K2 * VF[1][0][k] for k in range(7)]; sF = [VF[0][1][k] - K2 * VF[1][1][k] for k in range(7)]
    rf = roots_on_circle(aF, sF)
    mapped = np.sort(np.mod(2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(ru / 2)), 2 * np.pi))
    tot += 1
    if len(mapped) == len(rf) and (len(rf) == 0 or np.max(np.abs(np.angle(np.exp(1j * (mapped - rf))))) < 1e-5):
        agree += 1
    else:
        bad.append((e, l, w, s, len(ru), len(rf)))
print(f"  agreement on {agree}/{tot} random configurations (count and location)")
for b in bad[:8]:
    print("   mismatch e=%.4f l=%.4f w=%.4f s=%.4f: #u-roots %d, #f-roots %d" % b)

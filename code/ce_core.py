"""Circle-ellipse Keplerian distance: critical points of d^2 on the torus.

Setting (common focus at the origin, ellipse semi-major axis a = 1):
  circle  C(v) = s (cos v, sin v, 0),                     s = R / a > 0
  ellipse X(u) = (cos u - e) P + beta sin u Q,            beta = sqrt(1 - e^2)
  P = (cos w, sin w cos I, sin w sin I),  Q = (-sin w, cos w cos I, cos w sin I)
  d^2(u, v) = |X(u)|^2 + s^2 - 2 s (X_x cos v + X_y sin v).

For fixed u the critical v are the nearest / farthest points of the circle, so the
critical points of d^2 are the critical points of g_-(u) and g_+(u),
  g_-+(u) = |X|^2 + s^2 -+ 2 s rho(u),   rho^2 = X_x^2 + X_y^2,
i.e. the solutions of k(u) = +s (branch g_-) and k(u) = -s (branch g_+), where
  k(u) = A rho / (rho^2)',   A = (|X|^2)' = 2 e sin u (1 - e cos u).
Squaring gives the degree-6 trigonometric polynomial
  Pol(u) = A^2 rho^2 - s^2 ((rho^2)')^2,
so there are at most 12 critical points.  Conjecture (Albouy-Cabral-Santos 2012,
Problem 8, second question): at most 10.
"""
import numpy as np


def rho2_parts(u, e, sig, w):
    """rho^2 and its u-derivative; sig = sin^2 I."""
    beta = np.sqrt(1.0 - e * e)
    r = 1.0 - e * np.cos(u)                      # |X| = r
    zb = (np.cos(u) - e) * np.sin(w) + beta * np.sin(u) * np.cos(w)   # X_z / sin I
    zbp = -np.sin(u) * np.sin(w) + beta * np.cos(u) * np.cos(w)
    rho2 = r * r - sig * zb * zb
    drho2 = 2.0 * e * np.sin(u) * r - 2.0 * sig * zb * zbp
    return rho2, drho2, r


def k_func(u, e, sig, w):
    rho2, drho2, r = rho2_parts(u, e, sig, w)
    A = 2.0 * e * np.sin(u) * r
    rho = np.sqrt(np.maximum(rho2, 0.0))
    return A * rho / drho2


def pol_coeffs(e, sig, w, s, n=64):
    """Coefficients c_{-6..6} of Pol(u) = sum c_k e^{iku} (exact up to rounding: Pol is a trig poly of degree 6)."""
    u = 2 * np.pi * np.arange(n) / n
    rho2, drho2, r = rho2_parts(u, e, sig, w)
    A = 2.0 * e * np.sin(u) * r
    vals = A * A * rho2 - s * s * drho2 * drho2
    c = np.fft.fft(vals) / n
    ks = np.arange(-6, 7)
    return np.array([c[k % n] for k in ks])


def unit_roots(e, sig, w, s, tol=1e-6):
    """Roots of z^6 Pol on/off the unit circle. Returns (angles of unit roots, min | log|z| | of the others)."""
    c = pol_coeffs(e, sig, w, s)
    # polynomial in z of degree 12: sum_{k=-6}^{6} c_k z^{k+6}; numpy wants highest degree first
    roots = np.roots(c[::-1])
    lz = np.log(np.abs(roots))
    on = np.abs(lz) < tol
    ang = np.sort(np.mod(np.angle(roots[on]), 2 * np.pi))
    off = np.abs(lz[~on])
    return ang, (off.min() if off.size else np.inf), roots


def count_branch(e, sig, w, s, m=200000):
    """Count critical points of g_- and g_+ by sign changes of their derivatives on a fine grid."""
    u = 2 * np.pi * (np.arange(m) + 0.5) / m
    rho2, drho2, r = rho2_parts(u, e, sig, w)
    A = 2.0 * e * np.sin(u) * r
    rho = np.sqrt(np.maximum(rho2, 1e-300))
    gm = A - s * drho2 / rho       # g_-'
    gp = A + s * drho2 / rho       # g_+'
    def sc(f):
        return int(np.sum(np.sign(f) != np.sign(np.roll(f, -1))))
    return sc(gm), sc(gp)


def nmax_over_s(e, sig, w, m=40000):
    """For fixed (e, sig, w): max over s>0 of #{u : |k(u)| = s}, and the maximizing s."""
    u = 2 * np.pi * (np.arange(m) + 0.5) / m
    kk = np.abs(k_func(u, e, sig, w))
    # candidate levels: midpoints between sorted local extremum values of |k|
    d = np.diff(np.concatenate([kk, kk[:1]]))
    ext = np.where(np.sign(d) != np.sign(np.roll(d, 1)))[0]
    vals = np.unique(np.concatenate([kk[ext], [0.0]]))
    vals = vals[np.isfinite(vals)]
    best, bs = 0, None
    levels = (vals[:-1] + vals[1:]) / 2 if vals.size > 1 else vals
    for lv in levels:
        if lv <= 0:
            continue
        f = kk - lv
        cnt = int(np.sum(np.sign(f) != np.sign(np.roll(f, -1))))
        if cnt > best:
            best, bs = cnt, lv
    return best, bs


def brute_force_2d(e, sig, w, s, grid=60, seed=0):
    """Independent check: Newton on grad d^2 over the torus from a grid of starts."""
    I = np.arcsin(np.sqrt(sig))
    beta = np.sqrt(1 - e * e)
    P = np.array([np.cos(w), np.sin(w) * np.cos(I), np.sin(w) * np.sin(I)])
    Q = np.array([-np.sin(w), np.cos(w) * np.cos(I), np.cos(w) * np.sin(I)])

    def X(u):
        return (np.cos(u) - e) * P + beta * np.sin(u) * Q

    def Xp(u):
        return -np.sin(u) * P + beta * np.cos(u) * Q

    def Xpp(u):
        return -np.cos(u) * P - beta * np.sin(u) * Q

    def C(v):
        return s * np.array([np.cos(v), np.sin(v), 0.0])

    def Cp(v):
        return s * np.array([-np.sin(v), np.cos(v), 0.0])

    def grad_hess(u, v):
        D = X(u) - C(v)
        gu = 2 * D @ Xp(u)
        gv = -2 * D @ Cp(v)
        huu = 2 * (Xp(u) @ Xp(u) + D @ Xpp(u))
        hvv = 2 * (Cp(v) @ Cp(v) + D @ C(v))
        huv = -2 * Xp(u) @ Cp(v)
        return np.array([gu, gv]), np.array([[huu, huv], [huv, hvv]])

    found = []
    for i in range(grid):
        for j in range(grid):
            z = np.array([2 * np.pi * (i + 0.5) / grid, 2 * np.pi * (j + 0.5) / grid])
            for _ in range(60):
                g, H = grad_hess(*z)
                try:
                    step = np.linalg.solve(H, g)
                except np.linalg.LinAlgError:
                    break
                z = z - step
                if np.linalg.norm(step) < 1e-13:
                    break
            g, H = grad_hess(*z)
            if np.linalg.norm(g) < 1e-10:
                z = np.mod(z, 2 * np.pi)
                if not any(np.linalg.norm(np.angle(np.exp(1j * (z - f)))) < 1e-6 for f in found):
                    found.append(z)
    return found

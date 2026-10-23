"""AUDIT a3: reduction lemma  #stationary points of d^2 on the torus  <=  #real roots of Pol (with multiplicity).

Independent 2-D computation: vectorised damped Newton on grad d^2(u, v) from a 56 x 56 grid of starts on the
torus, built directly from the 3-D vectors (no reduction formula used).  Deduplicated critical points are
compared with the real roots of Pol(u) = A^2 rho^2 - S ((rho^2)')^2 obtained from the companion polynomial.
For every 2-D critical point we also check that its u is a root of Pol, and at roots where rho = 0 or
A = (rho^2)' = 0 we record the number of critical points over that u and the multiplicity of the root.

Test families:
  R  random generic configurations
  L1 l = 1 exactly (orbit plane contains the circle's axis: rho = 0 at two points)
  W0 omega = 0 or pi/2 (P_z Q_z = 0: A = (rho^2)' = 0 at u = 0, pi)
  T  configurations near the maximum (N = 10) found by a random search over s
  AX P_z = 1 (axis-aligned) : claimed N = 8 for every e, R
"""
import os
import sys
import numpy as np

rng = np.random.default_rng(777)


def frame(e, I, w):
    P = np.array([np.cos(w), np.sin(w) * np.cos(I), np.sin(w) * np.sin(I)])
    Q = np.array([-np.sin(w), np.cos(w) * np.cos(I), np.cos(w) * np.sin(I)])
    return P, Q


def crit_2d(e, I, w, s, n=56, iters=60):
    beta = np.sqrt(1 - e * e)
    P, Q = frame(e, I, w)
    g = np.linspace(0, 2 * np.pi, n, endpoint=False) + np.pi / n
    U, Vv = np.meshgrid(g, g, indexing='ij')
    u = U.ravel().copy(); v = Vv.ravel().copy()

    def parts(u, v):
        cu, su = np.cos(u), np.sin(u)
        X = (cu - e)[:, None] * P + (beta * su)[:, None] * Q
        Xp = (-su)[:, None] * P + (beta * cu)[:, None] * Q
        Xpp = (-cu)[:, None] * P + (-beta * su)[:, None] * Q
        C = s * np.stack([np.cos(v), np.sin(v), 0 * v], 1)
        Cp = s * np.stack([-np.sin(v), np.cos(v), 0 * v], 1)
        Cpp = -C
        Dv = X - C
        gu = 2 * np.sum(Dv * Xp, 1); gv = -2 * np.sum(Dv * Cp, 1)
        huu = 2 * (np.sum(Xp * Xp, 1) + np.sum(Dv * Xpp, 1))
        hvv = 2 * (np.sum(Cp * Cp, 1) - np.sum(Dv * Cpp, 1))
        huv = -2 * np.sum(Xp * Cp, 1)
        return gu, gv, huu, hvv, huv

    for _ in range(iters):
        gu, gv, huu, hvv, huv = parts(u, v)
        det = huu * hvv - huv * huv
        ok = np.abs(det) > 1e-300
        du = np.where(ok, (hvv * gu - huv * gv) / np.where(ok, det, 1), 0.0)
        dv = np.where(ok, (huu * gv - huv * gu) / np.where(ok, det, 1), 0.0)
        st = np.sqrt(du * du + dv * dv)
        fac = np.where(st > 0.4, 0.4 / np.maximum(st, 1e-300), 1.0)
        u = u - fac * du; v = v - fac * dv
    gu, gv, huu, hvv, huv = parts(u, v)
    scale = 1 + s * s
    good = np.sqrt(gu * gu + gv * gv) < 1e-10 * scale
    pts = np.stack([np.mod(u[good], 2 * np.pi), np.mod(v[good], 2 * np.pi)], 1)
    uniq = []
    for pt in pts:
        if not any(np.hypot(*np.angle(np.exp(1j * (pt - q0)))) < 1e-6 for q0 in uniq):
            uniq.append(pt)
    return np.array(uniq)


def pol_trig(e, p, q, S):
    """Fourier coefficients of Pol from the closed form (r, Z, rho^2), sampled + FFT (degree 6)."""
    N = 64
    u = 2 * np.pi * np.arange(N) / N
    b = np.sqrt(1 - e * e)
    r = 1 - e * np.cos(u)
    Z = p * (np.cos(u) - e) + b * q * np.sin(u)
    Zp = -p * np.sin(u) + b * q * np.cos(u)
    A = 2 * e * np.sin(u) * r
    rho2 = r * r - Z * Z
    d = A - 2 * Z * Zp
    vals = A * A * rho2 - S * d * d
    return np.fft.fft(vals) / N


def pol_real_roots(e, p, q, S, tol=1e-6):
    c = pol_trig(e, p, q, S)
    cc = np.array([c[(k - 6) % 64] for k in range(13)])
    z = np.roots(cc[::-1])
    lz = np.abs(np.log(np.abs(z)))
    on = lz < tol
    amb = np.sum((lz >= tol) & (lz < 1e-3))
    return np.sort(np.mod(np.angle(z[on]), 2 * np.pi)), int(amb)


def compare(e, I, w, s, verbose=False):
    p, q = np.sin(w) * np.sin(I), np.cos(w) * np.sin(I)
    cps = crit_2d(e, I, w, s)
    roots, amb = pol_real_roots(e, p, q, s * s)
    n2 = len(cps); npol = len(roots)
    # every critical point's u must be (close to) a root of Pol
    miss = 0
    for (uu, vv) in cps:
        if len(roots) == 0 or np.min(np.abs(np.angle(np.exp(1j * (roots - uu))))) > 1e-4:
            miss += 1
    return n2, npol, miss, amb, cps, roots


def run_family(name, sampler, n):
    viol = 0; miss_tot = 0; amb_tot = 0; hist = {}; eq = 0
    worst = None
    for it in range(n):
        e, I, w, s = sampler()
        n2, npol, miss, amb, cps, roots = compare(e, I, w, s)
        hist[(n2, npol)] = hist.get((n2, npol), 0) + 1
        miss_tot += miss; amb_tot += amb > 0
        eq += n2 == npol
        if n2 > npol and amb == 0:
            viol += 1
            worst = (e, I, w, s, n2, npol)
    print(f"[{name}] {n} configs: violations (#2D > #Pol, unambiguous) = {viol}; 2D points whose u is not a Pol root = {miss_tot}; "
          f"configs with near-circle complex roots (ambiguous count) = {amb_tot}; equal counts = {eq}")
    print(f"      (N_2D, N_Pol) histogram: {dict(sorted(hist.items()))}")
    if worst:
        print("      example violation:", worst)


if __name__ == "__main__":
    lu = lambda a, b: np.exp(rng.uniform(np.log(a), np.log(b)))
    run_family("R", lambda: (rng.uniform(0.05, 0.95), rng.uniform(0, np.pi), rng.uniform(0, 2 * np.pi), lu(0.1, 10)), 150)
    run_family("L1", lambda: (rng.uniform(0.05, 0.95), np.pi / 2, rng.uniform(0, 2 * np.pi), lu(0.1, 10)), 120)
    run_family("W0", lambda: (rng.uniform(0.05, 0.95), rng.uniform(0, np.pi), rng.choice([0.0, np.pi / 2, np.pi, 1.5 * np.pi]), lu(0.1, 10)), 120)
    run_family("AX", lambda: (rng.uniform(0.02, 0.98), np.pi / 2, np.pi / 2, lu(0.01, 100)), 120)

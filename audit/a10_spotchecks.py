"""AUDIT a10: numerical spot checks of the claims in STATEMENT.md and the working notes (independent code).

 S1 symmetry/reduction of the orientation: N_2D invariant under Omega (rotation about the circle axis),
    I -> pi - I, omega -> -omega, omega -> omega + pi.
 S2 axis-aligned case: H1 = (1+e) e sin f - kappa D (cos f + e) has exactly 2 simple zeros (e in (0,1), kappa > 0).
 S3 parabola-circle problem: max real roots of the degree-10 polynomial, and min relative gap of the complex pair.
 S4 branch counts: max #crit(g_-) and #crit(g_+), joint histogram.
 S5 corner probes C1, C1', C2, C2b, C3 (max root count, strict and loose tolerance).
"""
import sys
import numpy as np
from numpy.polynomial import polynomial as Pn
from scipy.optimize import minimize
sys.path.insert(0, '.')
from a3_reduction import crit_2d
from nbport import phipsi_coeffs

rng = np.random.default_rng(31337)
lu = lambda a, b: np.exp(rng.uniform(np.log(a), np.log(b)))
which = sys.argv[1:] or ['S1', 'S2', 'S3', 'S4', 'S5']


# ---------- S1 ----------
def crit_general(e, Om, I, w, s, n=56):
    """2-D Newton count with a general Omega: rotate the circle by -Omega instead (equivalent)."""
    return len(crit_2d(e, I, w, s, n=n)) if Om == 0 else len(crit_rot(e, Om, I, w, s, n))


def crit_rot(e, Om, I, w, s, n):
    # ellipse with full (Omega, I, omega) frame; circle in the xy-plane
    cO, sO = np.cos(Om), np.sin(Om)
    P = np.array([cO * np.cos(w) - sO * np.sin(w) * np.cos(I), sO * np.cos(w) + cO * np.sin(w) * np.cos(I), np.sin(w) * np.sin(I)])
    Q = np.array([-cO * np.sin(w) - sO * np.cos(w) * np.cos(I), -sO * np.sin(w) + cO * np.cos(w) * np.cos(I), np.cos(w) * np.sin(I)])
    beta = np.sqrt(1 - e * e)
    g = np.linspace(0, 2 * np.pi, n, endpoint=False) + np.pi / n
    U, V = np.meshgrid(g, g, indexing='ij'); u = U.ravel().copy(); v = V.ravel().copy()

    def parts(u, v):
        cu, su = np.cos(u), np.sin(u)
        X = (cu - e)[:, None] * P + (beta * su)[:, None] * Q
        Xp = (-su)[:, None] * P + (beta * cu)[:, None] * Q
        Xpp = (-cu)[:, None] * P + (-beta * su)[:, None] * Q
        C = s * np.stack([np.cos(v), np.sin(v), 0 * v], 1); Cp = s * np.stack([-np.sin(v), np.cos(v), 0 * v], 1)
        Dv = X - C
        return (2 * np.sum(Dv * Xp, 1), -2 * np.sum(Dv * Cp, 1), 2 * (np.sum(Xp * Xp, 1) + np.sum(Dv * Xpp, 1)),
                2 * (np.sum(Cp * Cp, 1) + np.sum(Dv * C, 1)), -2 * np.sum(Xp * Cp, 1))

    for _ in range(60):
        gu, gv, huu, hvv, huv = parts(u, v)
        det = huu * hvv - huv * huv
        ok = np.abs(det) > 1e-300
        du = np.where(ok, (hvv * gu - huv * gv) / np.where(ok, det, 1), 0); dv = np.where(ok, (huu * gv - huv * gu) / np.where(ok, det, 1), 0)
        st = np.hypot(du, dv); fac = np.where(st > 0.4, 0.4 / np.maximum(st, 1e-300), 1.0)
        u -= fac * du; v -= fac * dv
    gu, gv, *_ = parts(u, v)
    good = np.hypot(gu, gv) < 1e-10 * (1 + s * s)
    pts = np.stack([np.mod(u[good], 2 * np.pi), np.mod(v[good], 2 * np.pi)], 1)
    uniq = []
    for pt in pts:
        if not any(np.hypot(*np.angle(np.exp(1j * (pt - q0)))) < 1e-6 for q0 in uniq):
            uniq.append(pt)
    return uniq


if 'S1' in which:
    bad = 0
    for it in range(40):
        e = rng.uniform(0.05, 0.95); I = rng.uniform(0.05, np.pi - 0.05); w = rng.uniform(0, 2 * np.pi); s = lu(0.2, 5)
        base = len(crit_rot(e, 0.0, I, w, s, 56))
        others = [len(crit_rot(e, rng.uniform(0, 2 * np.pi), I, w, s, 56)), len(crit_rot(e, 0.0, np.pi - I, w, s, 56)),
                  len(crit_rot(e, 0.0, I, -w, s, 56)), len(crit_rot(e, 0.0, I, w + np.pi, s, 56))]
        if any(o != base for o in others):
            bad += 1
            print("   S1 mismatch", e, I, w, s, base, others)
    print(f"S1 orientation symmetries (Omega, I->pi-I, w->-w, w->w+pi): {40 - bad}/40 configurations invariant")

# ---------- S2 ----------
if 'S2' in which:
    f = np.linspace(-np.pi, np.pi, 200001)[:-1]
    hist = {}; mind = np.inf
    for it in range(20000):
        e = rng.uniform(1e-3, 1 - 1e-3) if it % 2 else 1 - lu(1e-6, 0.5); kap = lu(1e-4, 1e4)
        H1 = (1 + e) * e * np.sin(f) - kap * (1 + e * np.cos(f)) * (np.cos(f) + e)
        n = int(np.sum(np.sign(H1) != np.sign(np.roll(H1, -1))))
        hist[n] = hist.get(n, 0) + 1
    print("S2 zeros of H1 over 20000 (e, kappa) incl. e -> 1:", dict(sorted(hist.items())))


# ---------- S3 ----------
def parab_poly(k, l, w):
    pz, qz = l * np.sin(w), l * np.cos(w)
    X2 = np.array([1.0, 0, 2.0, 0, 1.0])                 # (1 + w^2)^2
    Z = np.array([pz, 2 * qz, -pz])
    rho2 = Pn.polysub(X2, Pn.polymul(Z, Z))
    A = Pn.polyder(X2); dr = Pn.polyder(rho2)
    return Pn.polysub(Pn.polymul(Pn.polymul(A, A), rho2), k * k * Pn.polymul(dr, dr))


def parab_count(k, l, w):
    c = np.trim_zeros(parab_poly(k, l, w), 'b')
    r = np.roots(c[::-1])
    rel = np.abs(r.imag) / (1 + np.abs(r.real))
    real = rel < 1e-7
    return int(real.sum()), (np.min(rel[~real]) if (~real).any() else np.inf)


if 'S3' in which:
    hist = {}; cands = []
    for it in range(30000):
        k = lu(1e-3, 1e3); l = rng.uniform(0, 1); w = rng.uniform(0, np.pi / 2)
        n, g = parab_count(k, l, w)
        hist[n] = hist.get(n, 0) + 1
        if n >= 8:
            cands.append((g, k, l, w))
    cands.sort()
    best = []
    for g, k, l, w in cands[:60]:
        fun = lambda x: (lambda n_g: n_g[1] if n_g[0] >= 8 else 10.0)(parab_count(np.exp(x[0]), np.clip(x[1], 0, 1), x[2]))
        r = minimize(fun, [np.log(k), l, w], method='Nelder-Mead', options=dict(xatol=1e-10, fatol=1e-12, maxiter=3000))
        best.append((r.fun, np.exp(r.x[0]), np.clip(r.x[1], 0, 1), r.x[2]))
    best.sort()
    print("S3 parabola-circle real-root histogram (30000 random):", dict(sorted(hist.items())))
    print("S3 smallest relative gap of the complex pair (>= 8 real):", ["gap=%.4f k=%.4f l=%.4f w=%.4f" % b for b in best[:3]])


# ---------- S4 ----------
def branch_counts(e, l, w, s, m=200000):
    u = 2 * np.pi * (np.arange(m) + 0.5) / m
    b = np.sqrt(1 - e * e); p, q = l * np.sin(w), l * np.cos(w)
    r = 1 - e * np.cos(u); Z = p * (np.cos(u) - e) + b * q * np.sin(u); Zp = -p * np.sin(u) + b * q * np.cos(u)
    A = 2 * e * np.sin(u) * r; rho2 = r * r - Z * Z; d = A - 2 * Z * Zp
    rho = np.sqrt(np.maximum(rho2, 0))
    Fm = A * rho - s * d; Fp = A * rho + s * d          # zeros of g_-' and g_+' (rho > 0)
    sc = lambda F: int(np.sum(np.sign(F) != np.sign(np.roll(F, -1))))
    return sc(Fm), sc(Fp)


if 'S4' in which:
    hist = {}; mm = (0, None); mp_ = (0, None)
    for it in range(12000):
        e = 1 - lu(1e-4, 0.5) if it % 3 == 0 else rng.uniform(0.0, 0.999)
        l = rng.uniform(0, 1); w = rng.uniform(0, np.pi / 2); s = lu(1e-3, 1e2)
        nm, npl = branch_counts(e, l, w, s, m=40000)
        hist[(nm, npl)] = hist.get((nm, npl), 0) + 1
        if nm > mm[0]: mm = (nm, (e, l, w, s))
        if npl > mp_[0]: mp_ = (npl, (e, l, w, s))
    # targeted: maximise each branch count over s for random (e, l, w)
    for it in range(1500):
        e = rng.uniform(0.2, 0.99); l = rng.uniform(0.3, 1); w = rng.uniform(0, np.pi / 2)
        for s in np.exp(np.linspace(np.log(0.02), np.log(20), 60)):
            nm, npl = branch_counts(e, l, w, s, m=20000)
            if nm > mm[0]: mm = (nm, (e, l, w, s))
            if npl > mp_[0]: mp_ = (npl, (e, l, w, s))
    print("S4 max #crit(g_-) =", mm, "; max #crit(g_+) =", mp_)
    print("S4 joint histogram (random part):", dict(sorted(hist.items())))


# ---------- S5 ----------
def count_q(e, l, w, S):
    t = S / (1 + S)
    (aP, sP), (aS, sS) = phipsi_coeffs(0, np.arcsin(e), l, w)
    a = (1 - t) * aP - t * aS; s = (1 - t) * sP - t * sS
    c = np.zeros(13, complex); c[6] = a[0]
    for k in range(1, 7):
        c[6 + k] = (a[k] - 1j * s[k]) / 2; c[6 - k] = (a[k] + 1j * s[k]) / 2
    mx = np.abs(c).max(); top = 12
    while top > 0 and abs(c[top]) < 1e-14 * mx: top -= 1
    low = 0
    while low < top and abs(c[low]) < 1e-14 * mx: low += 1
    z = np.roots(c[low:top + 1][::-1]); lz = np.abs(np.log(np.abs(z)))
    return int(np.sum(lz < 1e-7)), int(np.sum(lz < 1e-3))


if 'S5' in which:
    probes = {
        "C1  (e,l small)": lambda: (lu(1e-4, 0.05), lu(1e-4, 0.2), rng.uniform(0, np.pi / 2), lu(1e-3, 1e3)),
        "C1' (e small, s small)": lambda: (lu(1e-4, 0.05), rng.uniform(0, 1), rng.uniform(0, np.pi / 2), lu(1e-8, 1e-2)),
        "C2  (e->1, s small)": lambda: (1 - lu(1e-6, 1e-2), rng.uniform(0, 1), rng.uniform(0, np.pi / 2), lu(1e-10, 1e-2)),
        "C2b (e->1, s any)": lambda: (1 - lu(1e-6, 1e-2), rng.uniform(0, 1), rng.uniform(0, np.pi / 2), lu(1e-4, 1e4)),
        "C3  (e->1, P_z->1)": lambda: (1 - lu(1e-6, 1e-2), 1 - lu(1e-6, 1e-2), np.pi / 2 - lu(1e-6, 1e-2), lu(1e-4, 1e4)),
    }
    for name, smp in probes.items():
        hs = {}; hl = {}
        for _ in range(8000):
            e, l, w, S = smp()
            a, b = count_q(e, l, w, S)
            hs[a] = hs.get(a, 0) + 1; hl[b] = hl.get(b, 0) + 1
        print(f"S5 {name}: strict {dict(sorted(hs.items()))}  loose(|log|z||<1e-3) {dict(sorted(hl.items()))}")

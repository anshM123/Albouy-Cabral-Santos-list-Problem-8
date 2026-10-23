"""AUDIT a16: empirical soundness of the new engines (run with the .venv python, NUMBA_CACHE_DIR in audit/).

Part 1  nbeng_C1 / nbeng_C1p / nbeng_C1pp .certify_box: random boxes in each chart domain; when certified with bound
        B, the real roots of the chart function (1-x3) F_A - x3 F_B (independent evaluator of a15, eps > 0 samples)
        are counted at sampled points (near-circle roots counted as real: upward bias) and must be <= B.
Part 2  nb_corner.certify_corner_box: random corner boxes (e in [0.98,1], l in [0,0.9], w in [0,pi/2],
        mu in [0, 0.05], phi0 = 0.9, phi_in = 0.12).  When ok, at sampled points with e < 1, mu > 0:
          #zeros F_- on |phi| <= phi0                    <= nm
          #zeros F_+ on |phi| <= phi_in                  <= 1       (Lemma N2)
          #zeros F_+ on phi_in <= |phi| <= phi0          <= po
          #real zeros of P = F_- F_+ on |phi| >= phi0    <= aw
          total stationary points (2-D Newton, R = 1)   <= tot
        zeros counted as sign changes on fine grids (F_-, F_+ are smooth there since sigma > 0).
usage: python a16_engines.py N1 N2 SEED
"""
import os
import sys
import math
import time
import itertools
import importlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, 'numba_cache'))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(HERE, '..', 'code'))
sys.path.insert(0, HERE)

N1 = int(sys.argv[1]) if len(sys.argv) > 1 else 60
N2 = int(sys.argv[2]) if len(sys.argv) > 2 else 40
rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 16)

# ---------------- independent evaluator (copied from a15) ----------------
def phipsi(e, l, w, N=32):
    u = 2 * np.pi * np.arange(N) / N
    b = np.sqrt(1 - e * e)
    p, q = l * np.sin(w), l * np.cos(w)
    r = 1 - e * np.cos(u); Z = p * (np.cos(u) - e) + b * q * np.sin(u); Zp = -p * np.sin(u) + b * q * np.cos(u)
    A = 2 * e * np.sin(u) * r; rho2 = r * r - Z * Z; d = A - 2 * Z * Zp
    return A * A * rho2, d * d, u


def proj(F, u):
    N = len(u); a = np.zeros(7); s = np.zeros(7)
    a[0] = F.sum() / N
    for k in range(1, 7):
        a[k] = 2 * (F * np.cos(k * u)).sum() / N; s[k] = 2 * (F * np.sin(k * u)).sum() / N
    return a, s


def chart_coeffs(chart, x0, x1, x2, x3):
    if chart == 'C1':
        e = x0 * math.cos(x1); l = math.sqrt(x0 * math.sin(x1))
        Phi, Psi, u = phipsi(e, l, x2); FA, FB = Phi / x0 ** 2, Psi / x0 ** 2
    elif chart == 'C1p':
        e = x0; Phi, Psi, u = phipsi(e, x1, x2); FA, FB = Phi / e ** 2, x0 ** 2 * Phi / e ** 2 + Psi
    else:
        e = x0 * x1 ** 2; Phi, Psi, u = phipsi(e, x1, x2); FA, FB = Phi / e ** 2, x0 ** 2 * Phi / e ** 2 + Psi / x1 ** 4
    (aA, sA), (aB, sB) = proj(FA, u), proj(FB, u)
    return (1 - x3) * aA - x3 * aB, (1 - x3) * sA - x3 * sB


def count_trig(a, s):
    c = np.zeros(13, complex); c[6] = a[0]
    for k in range(1, 7):
        c[6 + k] = (a[k] - 1j * s[k]) / 2; c[6 - k] = (a[k] + 1j * s[k]) / 2
    mx = np.abs(c).max()
    top = 12
    while top > 0 and abs(c[top]) < 1e-15 * mx: top -= 1
    low = 0
    while low < top and abs(c[low]) < 1e-15 * mx: low += 1
    z = np.roots(c[low:top + 1][::-1])
    return int(np.sum(np.abs(np.log(np.abs(z))) < 1e-3))


# ---------------- Part 1 ----------------
DOMS = {'C1': [(0.0, 0.26), (0.0, math.atan(2.0)), (0.0, math.pi / 2), (0.0, 1.0)],
        'C1p': [(0.0, 0.2), (0.4, 1.0), (0.0, math.pi / 2), (0.0, 1.0)],
        'C1pp': [(0.0, 0.5), (0.0, 0.52), (0.0, math.pi / 2), (0.0, 1.0)]}
t0 = time.time()
for chart in ('C1', 'C1p', 'C1pp'):
    eng = importlib.import_module('nbeng_' + chart)
    st = dict(boxes=0, cert=0, viol=0); hist = {}; ex = None
    for it in range(N1):
        dom = DOMS[chart]
        c = [rng.uniform(*dom[i]) for i in range(4)]
        if it % 3 == 0:
            c[0] = dom[0][0] + 10 ** rng.uniform(-4, -1) * (dom[0][1] - dom[0][0])     # near the blow-up face
        r = np.array([dom[i][1] - dom[i][0] for i in range(4)]) * 0.01 * np.exp(rng.uniform(-3, 0, 4))
        bx = []
        for i in range(4):
            bx += [max(c[i] - r[i], dom[i][0]), min(c[i] + r[i], dom[i][1])]
        ok, B = eng.certify_box(*bx, 6, 64, 8)
        st['boxes'] += 1
        if not ok:
            continue
        st['cert'] += 1; hist[int(B)] = hist.get(int(B), 0) + 1
        pts = list(itertools.product(*[bx[2 * i:2 * i + 2] for i in range(4)]))
        pts += [tuple(rng.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(4)) for _ in range(50)]
        for P_ in pts:
            if P_[0] <= 1e-7 or (chart == 'C1pp' and P_[1] <= 1e-3):
                continue                                  # evaluator divides by eps (or e = epp l^2)
            a, s = chart_coeffs(chart, *P_)
            n = count_trig(a, s)
            if n > B:
                st['viol'] += 1; ex = (bx, P_, n, int(B))
    print(f"[{chart}] {st}, certified bounds {dict(sorted(hist.items()))}" + (f"  EXAMPLE {ex}" if ex else ""), flush=True)
print(f"part 1: {time.time()-t0:.0f}s", flush=True)

# ---------------- Part 2 ----------------
import nb_corner
PHI0, PHI_IN = 0.9, 0.12


def Fpm(f, e, l, w, mu):
    p, q = l * np.sin(w), l * np.cos(w)
    zeta = p * np.cos(f) + q * np.sin(f); zp = -p * np.sin(f) + q * np.cos(f)
    sig = np.sqrt(np.maximum(1 - zeta ** 2, 0)); D = 1 + e * np.cos(f)
    B = -2 * zeta * zp + 2 * e * (np.sin(f) - q * zeta)
    T = mu * 2 * e * (1 + e) * np.sin(f) * sig
    return T - D * B, T + D * B


def nzero(v):
    return int(np.sum(np.sign(v[:-1]) != np.sign(v[1:])))


st = dict(boxes=0, cert=0, viol_nm=0, viol_n2=0, viol_po=0, viol_aw=0, viol_tot=0)
hist = {}; exs = []
t0 = time.time()
fn = np.linspace(math.pi - PHI0, math.pi + PHI0, 400001)
fi = np.linspace(math.pi - PHI_IN, math.pi + PHI_IN, 200001)
fo1 = np.linspace(math.pi + PHI_IN, math.pi + PHI0, 100001); fo2 = np.linspace(math.pi - PHI0, math.pi - PHI_IN, 100001)
fa = np.linspace(math.pi + PHI0, 3 * math.pi - PHI0, 400001)
for it in range(N2):
    c = [rng.uniform(0.98, 1.0), rng.uniform(0, 0.9), rng.uniform(0, math.pi / 2), rng.uniform(0, 0.05)]
    r = np.array([0.02, 0.9, math.pi / 2, 0.05]) * 0.02 * np.exp(rng.uniform(-2, 0, 4))
    lims = [(0.98, 1.0), (0.0, 0.9), (0.0, math.pi / 2), (0.0, 0.05)]
    bx = []
    for i in range(4):
        bx += [max(c[i] - r[i], lims[i][0]), min(c[i] + r[i], lims[i][1])]
    ok, tot, nm, po, aw = nb_corner.certify_corner_box(*bx, PHI0, PHI_IN)
    st['boxes'] += 1
    if not ok:
        continue
    st['cert'] += 1; hist[int(tot)] = hist.get(int(tot), 0) + 1
    for _ in range(12):
        e = rng.uniform(bx[0], min(bx[1], 1 - 1e-9)); l = rng.uniform(bx[2], bx[3]); w = rng.uniform(bx[4], bx[5])
        mu = rng.uniform(max(bx[6], 1e-7), bx[7])
        Fm, Fp = Fpm(fn, e, l, w, mu)
        if nzero(Fm) > nm: st['viol_nm'] += 1; exs.append(('nm', bx, (e, l, w, mu), nzero(Fm), nm))
        Fm_i, Fp_i = Fpm(fi, e, l, w, mu)
        if nzero(Fp_i) > 1: st['viol_n2'] += 1; exs.append(('N2', bx, (e, l, w, mu), nzero(Fp_i)))
        npo = nzero(Fpm(fo1, e, l, w, mu)[1]) + nzero(Fpm(fo2, e, l, w, mu)[1])
        if npo > po: st['viol_po'] += 1; exs.append(('po', bx, (e, l, w, mu), npo, po))
        Am, Ap = Fpm(fa, e, l, w, mu)
        naw = nzero(Am * Ap)
        if naw > aw: st['viol_aw'] += 1; exs.append(('aw', bx, (e, l, w, mu), naw, aw))
        if nzero(Fm) + nzero(Fp_i) + npo + naw > tot: st['viol_tot'] += 1
print(f"[corner] {st}, certified totals {dict(sorted(hist.items()))}, {time.time()-t0:.0f}s")
for x in exs[:6]:
    print("   VIOLATION", x)

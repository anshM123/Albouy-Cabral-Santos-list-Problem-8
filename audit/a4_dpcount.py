"""AUDIT a4: dp_count (nb_cert numba version and certify2_ce python version).

1. Exhaustive comparison with a brute-force minimal cover (same rules: start after the first arc with a certified
   sign of f, cover the other M-1 arcs, in cyclic order, by runs of consecutive arcs on which one f^{(j)} has one
   constant nonzero sign; a run costs j).
2. Soundness on real functions: random trigonometric polynomials of degree 6 with many real roots; signs of f^{(j)}
   on arcs are taken from dense sampling with a safety margin (only certain signs are kept); the true number of
   roots on the circle (with multiplicity, from the companion polynomial) must be <= dp_count.
3. Wrap-around: the result must not depend on a rotation of the arc indexing beyond the choice of start
   (dp_count of a rotated pattern is a valid bound too; record the spread).
Run with the .venv python.
"""
import os
import sys
import itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, 'numba_cache'))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(HERE, '..', 'code'))
import nb_cert
from certify2_ce import dp_count as dp_py

rng = np.random.default_rng(99)


def brute(sg):
    J1, M = sg.shape
    nz = [i for i in range(M) if sg[0, i] != 0]
    if not nz:
        return 999
    st = nz[0]
    order = [(st + 1 + i) % M for i in range(M - 1)]
    best = [0] + [999] * (M - 1)          # best[n] = min cost to cover the first n arcs of `order`
    for n in range(1, M):
        for m in range(n):                # last run = order[m..n-1]
            run = order[m:n]
            costs = [j for j in range(J1) if sg[j, run[0]] != 0 and all(sg[j, a] == sg[j, run[0]] for a in run)]
            if costs:
                best[n] = min(best[n], best[m] + min(costs))
    return best[M - 1] if M > 1 else 0


# 1. exhaustive / random comparison
mism = 0; tests = 0
for M in range(1, 8):
    for J1 in (1, 2, 3):
        for _ in range(4000 if M > 3 else 500):
            sg = rng.integers(-1, 2, size=(J1, M)).astype(np.int64)
            a = nb_cert.dp_count(sg); b = dp_py(sg); c = min(brute(sg), 999)
            a = min(a, 999); b = min(b, 999)
            tests += 1
            if not (a == b == c):
                mism += 1
                if mism <= 5:
                    print("MISMATCH", sg.tolist(), a, b, c)
print(f"1. dp_count vs brute force: {tests} patterns, mismatches {mism}")


# 2. soundness on real trig polynomials
def tp_eval(a, s, u, j):
    ks = np.arange(7, dtype=float)
    ang = np.outer(u, ks) + j * np.pi / 2
    return (np.cos(ang) * ks ** j) @ a + (np.sin(ang) * ks ** j) @ s


def nroots(a, s):
    c = np.zeros(13, complex); c[6] = a[0]
    for k in range(1, 7):
        c[6 + k] = (a[k] - 1j * s[k]) / 2; c[6 - k] = (a[k] + 1j * s[k]) / 2
    z = np.roots(c[::-1])
    return int(np.sum(np.abs(np.log(np.abs(z))) < 1e-5))


viol = 0; tests = 0; tight = 0
for it in range(3000):
    # random polynomial with prescribed real roots: product of (cos u - cos r_i)-type factors + small noise
    nr = rng.integers(1, 7)
    roots = rng.uniform(0, 2 * np.pi, 2 * nr)
    u = 2 * np.pi * np.arange(64) / 64
    vals = np.prod([np.sin((u - r) / 2) for r in roots], axis=0) * rng.choice([-1, 1])
    vals = vals + rng.normal(0, 10 ** rng.uniform(-4, -1), size=(1,)) * np.cos(rng.integers(0, 7) * u + rng.uniform(0, 6))
    cf = np.fft.rfft(vals) / 64
    a = np.zeros(7); s = np.zeros(7); a[0] = cf[0].real; a[1:] = 2 * cf[1:7].real; s[1:] = -2 * cf[1:7].imag
    M = int(rng.choice([16, 32, 64, 128]))
    J1 = 7
    sg = np.zeros((J1, M), dtype=np.int64)
    for i in range(M):
        uu = np.linspace(2 * np.pi * i / M, 2 * np.pi * (i + 1) / M, 200)
        for j in range(J1):
            v = tp_eval(a, s, uu, j)
            sc = np.max(np.abs(v)) + 1e-300
            if np.all(v > 1e-6 * sc) and v.min() > 1e-9:
                sg[j, i] = 1
            elif np.all(v < -1e-6 * sc) and v.max() < -1e-9:
                sg[j, i] = -1
    bound = nb_cert.dp_count(sg)
    nr_true = nroots(a, s)
    tests += 1
    if nr_true > bound:
        viol += 1
        print("VIOLATION", nr_true, bound)
    tight += nr_true == bound
print(f"2. soundness on {tests} random trig polynomials: violations {viol}; bound attained in {tight} cases")

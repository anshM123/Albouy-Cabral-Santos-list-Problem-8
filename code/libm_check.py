"""Empirical support for the libm assumption of FLOAT_RIGOUR.md: numba math.sin / math.cos against 40-digit mpmath on
[-64, 64], including arguments of the engine form k*u + j*pi/2 and points next to multiples of pi/2.
Output (2026-10-23, Windows UCRT via numba): points 174467, max error 0.755 ulp of the exact value."""
import math, random, numpy as np, mpmath
from numba import njit
mpmath.mp.dps = 40
@njit(cache=False)
def sc(x):
    s = np.empty(x.size); c = np.empty(x.size)
    for i in range(x.size):
        s[i] = math.sin(x[i]); c[i] = math.cos(x[i])
    return s, c
random.seed(7)
xs = [random.uniform(-64, 64) for _ in range(150000)]
# arguments of the form fl(fl(k*u) + fl(j*pi/2)) as in the engines, plus points near multiples of pi/2
for k in range(1, 7):
    for j in range(0, 8):
        for _ in range(500):
            u = random.uniform(-3 * math.pi, 3 * math.pi)
            xs.append(k * u + j * math.pi / 2)
for m in range(-40, 41):
    for d in (0.0, 1e-12, -1e-12, 1e-8, -1e-8, 1e-4):
        xs.append(m * math.pi / 2 + d)
x = np.array([v for v in xs if abs(v) <= 64.0])
s, c = sc(x)
worst = 0.0; wabs = 0.0
for i in range(x.size):
    X = mpmath.mpf(float(x[i]))
    for val, ex in ((s[i], mpmath.sin(X)), (c[i], mpmath.cos(X))):
        err = abs(mpmath.mpf(float(val)) - ex)
        ulp = math.ulp(abs(float(ex))) if float(ex) != 0 else 5e-324
        worst = max(worst, float(err / ulp)); wabs = max(wabs, float(err))
print('points', x.size, 'max error in ulps of the exact value', worst, 'max abs error', wabs)

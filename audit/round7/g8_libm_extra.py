"""AUDIT-7 g8: libm functions used by the certificates beyond sin/cos (code/libm_check.py covers sin/cos only):
  cosh, sinh at k*L (winding weights, wind_cert.labels: |k L| <= 6*15 = 90, L in [-15, 0)),
  atan, atan2 (winding evaluation, wind_cert.arg_in_branch), tan (kq_run.T_HI = tan(0.45) + 1e-9), atan (bb_run PHI_A).
numba-compiled math.* against 40-digit mpmath; error in ulps of the exact value (cosh/sinh: of the value at the
float argument; the amplification of the rounding of k*L is a separate, documented term (AUDIT_4 R4-3))."""
import os, sys, math, random
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
from numba import njit
mp.mp.dps = 40


@njit(cache=False)
def f_hyp(x):
    a = np.empty(x.size); b = np.empty(x.size)
    for i in range(x.size):
        a[i] = math.cosh(x[i]); b[i] = math.sinh(x[i])
    return a, b


@njit(cache=False)
def f_at(y, x):
    a = np.empty(x.size); b = np.empty(x.size)
    for i in range(x.size):
        a[i] = math.atan2(y[i], x[i]); b[i] = math.atan(y[i] / x[i])
    return a, b


@njit(cache=False)
def f_misc():
    return math.tan(0.9 / 2), math.atan(3.0 * 0.37)


def ulp_err(val, ex):
    if ex == 0:
        return 0.0 if val == 0 else float('inf')
    return float(abs(mp.mpf(float(val)) - ex) / math.ulp(abs(float(ex))))


random.seed(11)
xs = []
for k in range(1, 7):
    for _ in range(4000):
        L = -random.uniform(1e-3, 15.0)
        xs.append(k * L)
xs += [random.uniform(-90, 90) for _ in range(20000)] + [-1e-3 * k for k in range(1, 7)]
x = np.array(xs)
ch, sh = f_hyp(x)
w_c = max(ulp_err(ch[i], mp.cosh(mp.mpf(float(x[i])))) for i in range(x.size))
w_s = max(ulp_err(sh[i], mp.sinh(mp.mpf(float(x[i])))) for i in range(x.size))
print(f'cosh, sinh on {x.size} arguments in [-90, 90] (incl. k*L, L in [-15, -1e-3]): max error {w_c:.3f} / {w_s:.3f} ulp')
ys = np.array([random.uniform(-1, 1) * 10 ** random.uniform(-12, 12) for _ in range(40000)])
xx = np.array([random.uniform(-1, 1) * 10 ** random.uniform(-12, 12) for _ in range(40000)])
a2, a1 = f_at(ys, xx)
w_a2 = max(ulp_err(a2[i], mp.atan2(mp.mpf(float(ys[i])), mp.mpf(float(xx[i])))) for i in range(ys.size))
w_a1 = max(abs(float(a1[i]) - float(mp.atan(mp.mpf(float(ys[i])) / mp.mpf(float(xx[i]))))) for i in range(ys.size))
print(f'atan2 on {ys.size} random pairs (magnitudes 1e-12..1e12): max error {w_a2:.3f} ulp; '
      f'atan(fl(y/x)) max absolute error {w_a1:.3e} rad (includes the division rounding)')
t, a = f_misc()
print('tan(0.45): error', ulp_err(t, mp.tan(mp.mpf(0.9) / 2)), 'ulp; atan(1.11):', ulp_err(a, mp.atan(mp.mpf(3.0 * 0.37))),
      'ulp (both are used with 1e-9 margins)')

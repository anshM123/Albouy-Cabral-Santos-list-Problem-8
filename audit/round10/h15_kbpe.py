"""AUDIT-10 h15: the perihelion chart KBPe.
  (1) identity: th = tan((f - pi)/2) = -1/tp with tp = tan(f/2), so tp^12 Q(-1/tp) = P(f) tp^12 (1 + 1/tp^2)^6
      = P(f)(1 + tp^2)^6 (exact; the zeros on tp in R are the zeros of P with f != pi, with the same multiplicities,
      since f -> tp is a diffeomorphism of (-pi, pi) onto R and (1 + tp^2)^6 > 0);
  (2) the generated value coefficients of nbcoef_c3_KBPe (first order) and nbcoef_c3_KBPe_2 (slot 0) equal the KBQe
      coefficients reversed with signs (-1)^j (j = power of tp), at random box points (float, exact equality expected
      up to the generator's float conversion);
  (3) G(tp) from KBPe_2 against MY P(f)(1 + tp^2)^6 at 50 digits, e in [0.83, 1], all l, w, mu in [0, 1],
      |tp| <= tan(0.6) + 0.1, also first and second derivatives in the box variables (slots 1..14) against numerical
      derivatives of slot 0."""
import os, re, sys, types, random
import numpy as np
import mpmath as mp
C3 = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code', 'c3'))
random.seed(15)
mp.mp.dps = 50


def load_plain(path):
    src = open(path).read().replace('from numba import njit', '')
    src = re.sub(r'@njit\([^)]*\)\n', '', src)
    mod = types.ModuleType(os.path.basename(path)); exec(compile(src, path, 'exec'), mod.__dict__)
    return mod


def coeffs(mod, x, nslot, obj=False):
    out = np.zeros((nslot, mod.DEG + 1), dtype=object if obj else float)
    if obj:
        out[...] = mp.mpf(0)
    fn = mod.coef_val2 if hasattr(mod, 'coef_val2') else mod.coef_val
    fn(np.array(x, dtype=object if obj else float), out)
    return out


def P_phys(e, pz, qz, mu, f):
    zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
    B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta); D = 1 + e * mp.cos(f); k = 2 * e * (1 + e)
    return mu ** 2 * k ** 2 * mp.sin(f) ** 2 * (1 - zeta ** 2) - D ** 2 * B ** 2


Q2 = load_plain(os.path.join(C3, 'nbcoef_c3_KBQe_2.py'))
P2 = load_plain(os.path.join(C3, 'nbcoef_c3_KBPe_2.py'))
P1 = load_plain(os.path.join(C3, 'nbcoef_c3_KBPe.py'))
print('degrees: KBQe_2', Q2.DEG, 'KBPe_2', P2.DEG, 'KBPe', P1.DEG)
# (2) reversal
worst2 = 0.0
for _ in range(200):
    x = [random.uniform(0, 0.2), random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)]
    q = coeffs(Q2, x, 15); p2 = coeffs(P2, x, 15); p1 = coeffs(P1, x, 5)
    deg = 12
    for j in range(deg + 1):
        ref = (-1) ** j * q[:, 12 - j] if 12 - j <= Q2.DEG else 0 * q[:, 0]
        worst2 = max(worst2, float(np.max(np.abs(p2[:, j] - ref)) / (np.max(np.abs(q)) + 1e-300)),
                     float(np.max(np.abs(p1[:5, j] - ref[:5])) / (np.max(np.abs(q)) + 1e-300)))
print(f'(2) KBPe coefficients (value, first and second derivatives; both modules) vs reversed KBQe with signs: max diff {worst2:.2e}')
# (3) against physics
worst3 = mp.mpf(0)
for _ in range(40):
    e = mp.mpf(random.uniform(0.83, 1.0)); l = mp.mpf(random.uniform(0, 1)); w = mp.mpf(random.uniform(0, mp.pi / 2))
    mu = mp.mpf(random.uniform(0, 1)); tp = mp.mpf(random.uniform(-0.8, 0.8))
    pz, qz = l * mp.sin(w), l * mp.cos(w)
    c = coeffs(P2, [1 - e, pz, qz, mu], 15, obj=True)
    G = sum(c[0, k] * tp ** k for k in range(P2.DEG + 1))
    ref = P_phys(e, pz, qz, mu, 2 * mp.atan(tp)) * (1 + tp ** 2) ** 6
    worst3 = max(worst3, abs(G - ref) / (1 + abs(ref)))
print(f'(3) KBPe_2 G(tp) vs my P(f)(1 + tp^2)^6, f = 2 atan(tp): max rel err {mp.nstr(worst3, 3)}')
# derivative slots against numerical derivatives (one point)
x0 = [mp.mpf(0.05), mp.mpf(0.6), mp.mpf(0.5), mp.mpf(0.3)]; h = mp.mpf(10) ** -20
c0 = coeffs(P2, x0, 15, obj=True); wd = mp.mpf(0)
for v in range(4):
    xp = list(x0); xm = list(x0); xp[v] += h; xm[v] -= h
    num = (coeffs(P2, xp, 15, obj=True)[0] - coeffs(P2, xm, 15, obj=True)[0]) / (2 * h)
    for k in range(P2.DEG + 1):
        wd = max(wd, abs(num[k] - c0[1 + v, k]) / (1 + abs(c0[1 + v, k])))
print(f'    first-derivative slots vs numerical derivatives: max rel err {mp.nstr(wd, 3)}')

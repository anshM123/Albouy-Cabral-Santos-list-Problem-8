"""AUDIT-9 h9: the face charts and the KQ near polynomial against MY OWN physics.
  (1) exactness of the face blow-up: Q(th = s y) with (e, mu) = (1 - s^2, s^2 m) (FE), (1 - s^2 E, s^2) (FM) and at
      th = +-3 s with (1 - s^2 E, s^2 m) (FIP/FIM) is divisible by s^6 (sympy, own Q from P = mu^2 k^2 sin^2 f sigma^2
      - D^2 B^2), and the quotient at s = 0 is not identically zero;
  (2) the generated second-order modules nbcoef_c3_{FE,FM,FIP,FIM}_2 (value coefficients, slot 0) and their first-order
      versions against G = Q / s^6 from my physics at 50 digits, at random box points;
  (3) nbcoef_c3_KBQe_2 against my Q(th) at e in [0.83, 0.96] (engine KQ used by KQ0/KQ0T away from e = 1)."""
import os, re, sys, types, random
import numpy as np
import mpmath as mp
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
C3 = os.path.normpath(os.path.join(HERE, '..', '..', 'code', 'c3'))
random.seed(9)
mp.mp.dps = 50


def P_phys(e, pz, qz, mu, f, M=mp):
    zeta = pz * M.cos(f) + qz * M.sin(f); dz = -pz * M.sin(f) + qz * M.cos(f)
    B = -2 * zeta * dz + 2 * e * (M.sin(f) - qz * zeta); D = 1 + e * M.cos(f); k = 2 * e * (1 + e)
    return mu ** 2 * k ** 2 * M.sin(f) ** 2 * (1 - zeta ** 2) - D ** 2 * B ** 2


def Q_phys(e, pz, qz, mu, th):
    return P_phys(e, pz, qz, mu, mp.pi + 2 * mp.atan(th)) * (1 + th ** 2) ** 6


# (1) symbolic exactness, own Q
th, s, y, m, E, pz, qz = sp.symbols('th s y m E pz qz', real=True)
ee, mu = sp.symbols('e mu', real=True)
cf = -(1 - th ** 2) / (1 + th ** 2); sf = -2 * th / (1 + th ** 2)
zeta = pz * cf + qz * sf; dzeta = -pz * sf + qz * cf
B = -2 * zeta * dzeta + 2 * ee * (sf - qz * zeta); D = 1 + ee * cf; k = 2 * ee * (1 + ee)
Q = sp.expand(sp.cancel((mu ** 2 * k ** 2 * sf ** 2 * (1 - zeta ** 2) - D ** 2 * B ** 2) * (1 + th ** 2) ** 6))
for name, sub in (('FE', {th: s * y, ee: 1 - s ** 2, mu: s ** 2 * m}), ('FM', {th: s * y, ee: 1 - s ** 2 * E, mu: s ** 2}),
                  ('FIP', {th: 3 * s, ee: 1 - s ** 2 * E, mu: s ** 2 * m}), ('FIM', {th: -3 * s, ee: 1 - s ** 2 * E, mu: s ** 2 * m})):
    W = sp.Poly(sp.expand(Q.subs(sub)), s)
    low = min(mo[0] for mo in W.monoms())
    lead = sp.factor(W.coeff_monomial(s ** low))
    print(f'(1) {name}: lowest power of s = {low}; leading coefficient = {lead}', flush=True)


# (2) generated modules vs physics
def load_plain(path):
    src = open(path).read().replace('from numba import njit', '')
    src = re.sub(r'@njit\([^)]*\)\n', '', src)
    mod = types.ModuleType(os.path.basename(path)); exec(compile(src, path, 'exec'), mod.__dict__)
    return mod


def gval(mod, nslot, x, r):
    out = np.zeros((nslot, mod.DEG + 1), dtype=object); out[...] = mp.mpf(0)
    (mod.coef_val2 if hasattr(mod, 'coef_val2') else mod.coef_val)(np.array([mp.mpf(v) for v in x], dtype=object), out)
    return sum(out[0, kk] * mp.mpf(r) ** kk for kk in range(mod.DEG + 1))


for ch in ('FE', 'FM', 'FIP', 'FIM'):
    worst = mp.mpf(0)
    for fn, nslot in ((f'nbcoef_c3_{ch}_2.py', 15), (f'nbcoef_c3_{ch}.py', 5)):
        mod = load_plain(os.path.join(C3, fn))
        for _ in range(12):
            sv = mp.mpf(10) ** mp.mpf(random.uniform(-6, mp.log10(0.02)))
            rat = mp.mpf(random.uniform(0, 1.001)); pzv = mp.mpf(random.uniform(0, 1)); qzv = mp.mpf(random.uniform(0, 1))
            if ch in ('FE', 'FM'):
                r = mp.mpf(random.uniform(-3, 3))
                e_, mu_ = (1 - sv ** 2, sv ** 2 * rat) if ch == 'FE' else (1 - sv ** 2 * rat, sv ** 2)
                th_ = sv * r
            else:
                r = mp.mpf(random.uniform(0, 1.001))
                e_, mu_ = 1 - sv ** 2 * rat, sv ** 2 * r
                th_ = (3 if ch == 'FIP' else -3) * sv
            ref = Q_phys(e_, pzv, qzv, mu_, th_) / sv ** 6
            val = gval(mod, nslot, [sv, rat, pzv, qzv], r)
            worst = max(worst, abs(val - ref) / (1 + abs(ref)))
    print(f'(2) {ch}: generated chart polynomial (both module versions) vs my Q / s^6: max rel err {mp.nstr(worst, 3)}', flush=True)

# (3) KBQe away from e = 1
mod = load_plain(os.path.join(C3, 'nbcoef_c3_KBQe_2.py'))
worst = mp.mpf(0)
for _ in range(30):
    e_ = mp.mpf(random.uniform(0.83, 0.96)); l_ = mp.mpf(random.uniform(0, 1)); w_ = mp.mpf(random.uniform(0, mp.pi / 2))
    mu_ = mp.mpf(random.uniform(0, 1)); th_ = mp.mpf(random.uniform(-0.75, 0.75))
    pzv, qzv = l_ * mp.sin(w_), l_ * mp.cos(w_)
    ref = Q_phys(e_, pzv, qzv, mu_, th_)
    val = gval(mod, 15, [1 - e_, pzv, qzv, mu_], th_)
    worst = max(worst, abs(val - ref) / (1 + abs(ref)))
print(f'(3) KBQe_2 vs my Q(th) at e in [0.83, 0.96], mu in [0, 1], |th| <= 0.75: max rel err {mp.nstr(worst, 3)}', flush=True)

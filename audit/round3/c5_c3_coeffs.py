"""AUDIT-3 c5: generated C3 coefficient files nbcoef_c3_<chart>.py vs my own physics at 100 digits.

For each chart: random box points x and root values r; G(r) = sum_k c_k(x) r^k from coef_val (pure python, mpmath
inputs) must equal Q_phys(tau)/norm (own chart maps, see c2_c3_replay_spot.chart_fn); gradients out[1+v] vs mpmath
central differences of out[0]; abs polynomials checked term by term (textual).
Also: the leading-order claim of c3_sym.py: G(y) at nu -> 0 equals P0hat(2y) (Ca), with an O(nu^2) remainder.
"""
import os, re, sys, types, random
import mpmath as mp
mp.mp.dps = 100
HERE = os.path.dirname(os.path.abspath(__file__))
C3 = os.path.join(HERE, '..', '..', 'code', 'c3')
sys.argv = [sys.argv[0], '0', '0', '0']           # c2 module parses argv at import
from c2_c3_replay_spot import chart_fn, Q_phys
random.seed(77)
CHARTS = ['Ca', 'CK', 'CE', 'Cm', 'PhiP', 'PhiM', 'L2aE', 'L2am', 'L2aYp', 'L2aYm', 'L2KE', 'L2Km', 'L2KYp', 'L2KYm', 'KBQe']


def load_plain(path):
    src = open(path).read().replace('from numba import njit', '')
    src = re.sub(r'@njit\([^)]*\)\n', '', src)
    mod = types.ModuleType(os.path.basename(path)); exec(compile(src, path, 'exec'), mod.__dict__)
    return mod, src


def abs_textual(src):
    bv = src.split('def coef_val')[1].split('def coef_abs')[0]; ba = src.split('def coef_abs')[1]
    lv = dict(re.findall(r'out\[([^\]]+)\] = (.*)', bv)); la = dict(re.findall(r'out\[([^\]]+)\] = (.*)', ba))
    bad = 0; n = 0
    for k in lv:
        tv = lv[k].split(' + '); ta = la[k].split(' + ')
        if len(tv) != len(ta):
            bad += 1; continue
        for a, b in zip(tv, ta):
            n += 1
            ca, ma = a.split(' * ', 1) if ' * ' in a else (a, ''); cb, mb = b.split(' * ', 1) if ' * ' in b else (b, '')
            if ma != mb or abs(float(ca)) != float(cb):
                bad += 1
    return bad, n


def sample_x(ch):
    if ch == 'KBQe':
        return [random.uniform(0, 0.01), random.uniform(0.9, 0.99), random.uniform(0.15, 0.7), random.uniform(0, 0.0505)]
    x = [random.uniform(0.01, 0.37)] + [random.uniform(0, 1.0) for _ in range(3)]
    if ch.startswith('L2'):
        x[2] = random.uniform(0.02, 1.0)
    return x


for ch in CHARTS:
    mod, src = load_plain(os.path.join(C3, f'nbcoef_c3_{ch}.py'))
    deg = mod.DEG
    worst = mp.mpf(0); worst_g = mp.mpf(0)
    for it in range(6):
        x = sample_x(ch)
        xm = [mp.mpf(v) for v in x]
        import numpy as np
        out = np.empty((5, deg + 1), dtype=object); out[...] = mp.mpf(0)
        mod.coef_val(xm, out)
        if ch == 'KBQe':
            fn = lambda r: Q_phys(1 - xm[0], xm[1], xm[2], xm[3], r); lo, hi = -0.5, 0.5
        else:
            fn, (lo, hi), _ = chart_fn(ch, x); lo, hi = float(lo), float(hi)
        vals = []
        for r in [mp.mpf(lo) + (mp.mpf(hi) - mp.mpf(lo)) * mp.mpf(k) / 6 for k in range(7)]:
            g = sum(out[0, k] * r ** k for k in range(deg + 1))
            ref = fn(r)
            vals.append((g, ref))
        sc = max(abs(ref) for _, ref in vals) + mp.mpf(10) ** -80
        worst = max(worst, max(abs(g - ref) for g, ref in vals) / sc)
        if it < 2:
            for v in range(4):
                h = mp.mpf(10) ** -30
                xp = list(xm); xp[v] += h; xmn = list(xm); xmn[v] -= h
                op = np.empty((5, deg + 1), dtype=object); op[...] = mp.mpf(0); mod.coef_val(xp, op)
                om = np.empty((5, deg + 1), dtype=object); om[...] = mp.mpf(0); mod.coef_val(xmn, om)
                for k in range(deg + 1):
                    num = (op[0, k] - om[0, k]) / (2 * h)
                    worst_g = max(worst_g, abs(num - out[1 + v, k]) / (1 + abs(out[1 + v, k])))
    bad, n = abs_textual(src)
    print(f"{ch:6s}: deg {deg}; value vs own physics max rel err {mp.nstr(worst, 3)}; gradient rel err {mp.nstr(worst_g, 3)};"
          f" abs-poly terms {n}, mismatches {bad}", flush=True)

# leading-order claim (c3_sym): G(y; nu) -> P0hat(2y) as nu -> 0 (chart Ca: a = 1, K, E, m ratios)
def P0hat(x, a, K, E, m):
    lam = 2 * K - a ** 2
    Bh = 2 * E * (x - a) - 2 * lam * x - 2 * a ** 2 * x + 3 * a * x ** 2 - x ** 3
    return 16 * m ** 2 * x ** 2 * (lam + (x - a) ** 2) - (E + x ** 2 / 2) ** 2 * Bh ** 2
K, E, m = mp.mpf('0.7'), mp.mpf('0.3'), mp.mpf('0.4')
for y in (mp.mpf('-1.3'), mp.mpf('0.4'), mp.mpf('2.2')):
    errs = []
    for nu in (mp.mpf('1e-2'), mp.mpf('1e-3'), mp.mpf('1e-4')):
        G = Q_phys(1 - nu ** 2 * E, 1 - nu ** 2 * K, nu, nu ** 3 * m, nu * y) / nu ** 10
        errs.append(G - P0hat(2 * y, 1, K, E, m))
    print(f"leading order at y={mp.nstr(y,3)}: G - P0hat(2y) = {[mp.nstr(v, 4) for v in errs]} (nu = 1e-2, 1e-3, 1e-4; ratio ~ nu^2: "
          f"{mp.nstr(errs[0] / errs[1], 4)}, {mp.nstr(errs[1] / errs[2], 4)})")

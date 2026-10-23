"""AUDIT-5 e4: D2b corner.
1. Symbolic: C1p function F = (1-tau) Phi/e^2 - tau (Phi + Psi) at e = 0, p = 0, q = 1 equals 4 s^2 c^2 (1 - 2 tau);
   its first-order part in (e, p, kappa = 1 - q, tau' = tau - 1/2) at the D2b point equals
   -8 tau' s^2 c^2 + 2 kappa (2 s^4 + 4 s^2 c^2) - 4 p s c^3 - 4 e s^2 c^3; the four functions are linearly independent.
2. Generated nbcoef_D2{E,P,K,Tp,Tm}.py vs my own physics at 50 digits (value), gradient consistency, abs-poly terms.
3. Cone coverage of the D2 cube: random points of {e, p, kappa <= SD2, |tau'| <= SD2} lie in some D2 chart domain."""
import os, re, sys, types, random
import sympy as sp
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__)); CODE = os.path.join(HERE, '..', '..', 'code')
u, e, p, k, tp = sp.symbols('u e p kappa taup', real=True)
s, c = sp.sin(u), sp.cos(u)
b = sp.sqrt(1 - e ** 2); q = 1 - k; tau = sp.Rational(1, 2) + tp
r = 1 - e * c; Z = p * (c - e) + q * b * s
rho2 = r ** 2 - Z ** 2; A = sp.diff(r ** 2, u)
PhiE = (2 * s * r) ** 2 * rho2                    # Phi / e^2
Phi = A ** 2 * rho2; Psi = sp.diff(rho2, u) ** 2
F = (1 - tau) * PhiE - tau * (Phi + Psi)


def zero(expr):
    return sp.simplify(sp.expand(sp.expand_trig(expr).rewrite(sp.exp))) == 0


F0 = F.subs({e: 0, p: 0, k: 0})
print("1. F at e=p=kappa=0 equals 4 s^2 c^2 (1 - 2 tau):", zero(F0 - 4 * s ** 2 * c ** 2 * (1 - 2 * tau)))
lin = sum(sp.diff(F, v).subs({e: 0, p: 0, k: 0, tp: 0}) * v for v in (e, p, k, tp))
claim = -8 * tp * s ** 2 * c ** 2 + 2 * k * (2 * s ** 4 + 4 * s ** 2 * c ** 2) - 4 * p * s * c ** 3 - 4 * e * s ** 2 * c ** 3
print("1. linear part equals the claim:", zero(lin - claim))
funcs = [s ** 2 * c ** 2, 2 * s ** 4 + 4 * s ** 2 * c ** 2, s * c ** 3, s ** 2 * c ** 3]
pts = [sp.Rational(j, 7) for j in range(1, 9)]
Mx = sp.Matrix([[f.subs(u, x).evalf(30) for f in funcs] for x in pts])
print("1. rank of the four functions sampled at 8 points:", Mx.rank(iszerofunc=lambda z: abs(z) < 1e-25))

# 2. generated coefficients vs physics at 50 digits
mp.mp.dps = 50
SD = mp.mpf(5e-3)


def load_plain(path):
    src = open(path).read().replace('from numba import njit', '')
    src = re.sub(r'@njit\([^)]*\)\n', '', src)
    mod = types.ModuleType(os.path.basename(path)); exec(compile(src, path, 'exec'), mod.__dict__)
    return mod, src


def fphys(ee, pp, qq, tt, uu):
    bb = mp.sqrt(1 - ee * ee); rr = 1 - ee * mp.cos(uu)
    Zz = pp * (mp.cos(uu) - ee) + qq * bb * mp.sin(uu); Zp = -pp * mp.sin(uu) + qq * bb * mp.cos(uu)
    Aa = 2 * ee * mp.sin(uu) * rr; rh = rr * rr - Zz * Zz; dd = Aa - 2 * Zz * Zp
    return (1 - tt) * Aa * Aa * rh - tt * dd * dd


random.seed(3)
import numpy as np
for ch in ('D2E', 'D2P', 'D2K', 'D2Tp', 'D2Tm'):
    mod, src = load_plain(os.path.join(CODE, f'nbcoef_{ch}.py'))
    names = mod.ATOMS; nvar = mod.NVAR; affine = mod.NF == 2
    worst = mp.mpf(0); worst_g = mp.mpf(0)
    for it in range(6):
        nu = mp.mpf(10) ** mp.mpf(random.uniform(-4, 0.004))
        r1, r2 = mp.mpf(random.uniform(0.01, 1.05)), mp.mpf(random.uniform(0, 1.05))
        if ch in ('D2Tp', 'D2Tm'):
            r3 = mp.mpf(random.uniform(0, 1.05)); sg = 1 if ch == 'D2Tp' else -1
            E_, P_, K_ = r1, r2, r3
            ee, pp, kk, tpv = SD * E_ * nu, SD * P_ * nu, SD * K_ * nu, sg * SD * nu
            x = [nu, E_, P_, K_]; x3 = None
        else:
            x3 = mp.mpf(random.random()); T = -mp.mpf(1.05) + 2 * mp.mpf(1.05) * x3
            if ch == 'D2E':
                ee, pp, kk = SD * nu, SD * r1 * nu, SD * r2 * nu
            elif ch == 'D2P':
                ee, pp, kk = SD * r1 * nu, SD * nu, SD * r2 * nu
            else:
                ee, pp, kk = SD * r1 * nu, SD * r2 * nu, SD * nu
            tpv = SD * T * nu; x = [nu, r1, r2, x3]
        tauv = mp.mpf(1) / 2 + tpv; tt = ee * ee * tauv / (1 - tauv)
        A_ = {'nu': nu, 'b': mp.sqrt(1 - ee * ee)}; A_['ib'] = 1 / A_['b']
        A_.update({'E': ee / (SD * nu)} if ch != 'D2E' else {})
        if ch == 'D2E': A_.update(P=r1, K=r2)
        if ch == 'D2P': A_.update(E=r1, K=r2)
        if ch == 'D2K': A_.update(E=r1, P=r2)
        if ch in ('D2Tp', 'D2Tm'): A_.update(E=E_, P=P_, K=K_)
        out = np.empty((nvar + 1, 2, 2, 7), dtype=object); out[...] = mp.mpf(0)
        mod.coef_val([A_[n] for n in names], out)
        norm = ee * ee * nu / (1 - tauv)
        for jj in range(5):
            uu = 2 * mp.pi * (jj + mp.mpf(1) / 3) / 5
            trig = lambda fi: sum(out[0, fi, 0, kk_] * mp.cos(kk_ * uu) + out[0, fi, 1, kk_] * mp.sin(kk_ * uu) for kk_ in range(7))
            val = (1 - x3) * trig(0) - x3 * trig(1) if affine else trig(0)
            ref = fphys(ee, pp, 1 - kk, tt, uu) / norm
            worst = max(worst, abs(val - ref) / (1 + abs(ref)))
    # abs-poly textual check
    bv = src.split('def coef_val')[1].split('def coef_abs')[0]; ba = src.split('def coef_abs')[1]
    lv = dict(re.findall(r'out\[([^\]]+)\] = (.*)', bv)); la = dict(re.findall(r'out\[([^\]]+)\] = (.*)', ba))
    badt = 0; nt = 0
    for key in lv:
        for a_, b_ in zip(lv[key].split(' + '), la[key].split(' + ')):
            nt += 1
            ca, ma = a_.split(' * ', 1) if ' * ' in a_ else (a_, ''); cb, mb = b_.split(' * ', 1) if ' * ' in b_ else (b_, '')
            badt += (ma != mb) or abs(float(ca)) != float(cb)
    print(f"2. {ch}: value vs physics (50 digits) max rel err {mp.nstr(worst, 3)}; abs-poly terms {nt}, mismatches {badt}")

# 3. cube coverage by the cones (scaled coords X = (e, p, kappa, tau')/SD, nu = max |X_i|, ratios <= 1.05)
viol = 0
for it in range(200000):
    X = [random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1), random.uniform(-1, 1)]
    nu = max(abs(v) for v in X)
    i = max(range(4), key=lambda j: abs(X[j]))
    ratios = [abs(X[j]) / nu for j in range(4) if j != i]
    if not (nu <= 1.01 and all(rr <= 1.05 for rr in ratios)):
        viol += 1
print("3. cube points not in any D2 cone chart domain:", viol, "of 200000")

"""AUDIT-2 b3: generated blow-up coefficient files nbcoef_<B1,B2,B3p,B3m,L2E,L2P,L2Ap,L2Am,L2Tp,L2Tm>.py.

1. VALUES vs PHYSICS at 50 digits: the chart function built from coef_val(atoms(x)) (affine charts:
   (1-x3) F_A - x3 F_B) is compared with f_phys(u) / norm, where f_phys = (1-t) Phi - t Psi is evaluated from the
   a-chart closed forms in mpmath (50 digits) at the physical point given by MY OWN implementation of each chart map.
   Reports max |F_chart - f/norm| / max|f/norm| over 8 u-values x N random box points (also points with tiny nu).
   The chart value itself is evaluated in mpmath from the generated polynomials (float constants are exact dyadics).
2. GRADIENTS: out[1+v] vs mpmath derivative (mp.diff) of out[0] along box variable v through atoms(x).
3. ABS POLYNOMIALS: textual term-by-term check that coef_abs has the same monomials with |coefficient|.
usage: python b3_bu_coeffs.py NPTS
"""
import os
import re
import sys
import types
import random
import mpmath as mp

CODE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code')
mp.mp.dps = 50
NPTS = int(sys.argv[1]) if len(sys.argv) > 1 else 6
random.seed(2323)
CHARTS = ['B1', 'B2', 'B3p', 'B3m', 'L2E', 'L2P', 'L2Ap', 'L2Am', 'L2Tp', 'L2Tm']
S0, S1, S2, S3 = mp.mpf('0.01'), mp.mpf('0.25'), mp.mpf('0.25'), mp.mpf('0.5')


def load_plain(path):
    src = open(path).read()
    src = src.replace('from numba import njit', '')
    src = re.sub(r'@njit\([^)]*\)\n', '', src)
    mod = types.ModuleType(os.path.basename(path))
    exec(compile(src, path, 'exec'), mod.__dict__)
    return mod, src


def chart_map(ch, x):
    """my own map: box point x -> (e, p, q, t, norm, atom dict, affine?)"""
    x = [mp.mpf(v) for v in x]
    A = {}
    if ch == 'B1':
        e, ph, qh, x3 = x; tauh = -2 + 4 * x3
        p, q, t, norm = e * ph, e * qh, mp.mpf(1) / 2 + e * tauh, e ** 3
        A.update(e=e, ph=ph, qh=qh)
    elif ch == 'B2':
        th, eh, w, x3 = x; tauh = -mp.mpf('0.6825') + 2 * mp.mpf('0.6825') * x3
        e = th * eh; l = mp.sqrt(e * th)
        p, q, t, norm = l * mp.sin(w), l * mp.cos(w), mp.mpf(1) / 2 + th * tauh, e ** 2 * th
        A.update(th=th, eh=eh, sw=mp.sin(w), cw=mp.cos(w))
    elif ch in ('B3p', 'B3m'):
        sg = 1 if ch == 'B3p' else -1
        eta, eh, thh, w = x
        e = eta * eh; th = eta * thh; l = mp.sqrt(e * th)
        p, q, t, norm = l * mp.sin(w), l * mp.cos(w), mp.mpf(1) / 2 + sg * eta, e ** 2 * eta
        A.update(eta=eta, eh=eh, thh=thh, sw=mp.sin(w), cw=mp.cos(w))
    else:
        nu = x[0]
        if ch == 'L2E':
            _, P, Aa, x3 = x; T = -mp.mpf('1.05') + mp.mpf('2.1') * x3
            e = S0 * nu; ph = S1 * P * nu; qh = 1 + S2 * Aa * nu; tauh = S3 * T * nu
            A.update(nu=nu, P=P, A=Aa)
        elif ch == 'L2P':
            _, E, Aa, x3 = x; T = -mp.mpf('1.05') + mp.mpf('2.1') * x3
            e = S0 * E * nu; ph = S1 * nu; qh = 1 + S2 * Aa * nu; tauh = S3 * T * nu
            A.update(nu=nu, E=E, A=Aa)
        elif ch in ('L2Ap', 'L2Am'):
            sg = 1 if ch == 'L2Ap' else -1
            _, E, P, x3 = x; T = -mp.mpf('1.05') + mp.mpf('2.1') * x3
            e = S0 * E * nu; ph = S1 * P * nu; qh = 1 + sg * S2 * nu; tauh = S3 * T * nu
            A.update(nu=nu, E=E, P=P)
        else:
            sg = 1 if ch == 'L2Tp' else -1
            _, E, P, Aa = x
            e = S0 * E * nu; ph = S1 * P * nu; qh = 1 + S2 * Aa * nu; tauh = sg * S3 * nu
            A.update(nu=nu, E=E, P=P, A=Aa)
        p, q, t, norm = e * ph, e * qh, mp.mpf(1) / 2 + e * tauh, e ** 3 * nu
    b = mp.sqrt(1 - e * e)
    A.update(b=b, ib=1 / b)
    affine = ch in ('B1', 'B2', 'L2E', 'L2P', 'L2Ap', 'L2Am')
    return e, p, q, t, norm, A, affine


def f_phys(e, p, q, t, u):
    b = mp.sqrt(1 - e * e)
    r = 1 - e * mp.cos(u)
    Z = p * (mp.cos(u) - e) + q * b * mp.sin(u)
    Zp = -p * mp.sin(u) + q * b * mp.cos(u)
    Aq = 2 * e * mp.sin(u) * r
    rho2 = r * r - Z * Z
    d = Aq - 2 * Z * Zp
    return (1 - t) * Aq * Aq * rho2 - t * d * d


def sample(ch):
    if ch == 'B1':
        return [random.uniform(1e-4, 0.07), random.uniform(0, 1.75), random.uniform(0, 1.75), random.random()]
    if ch == 'B2':
        return [random.uniform(1e-4, 0.2), random.uniform(1e-3, 0.3675), random.uniform(0, 1.5707963), random.random()]
    if ch in ('B3p', 'B3m'):
        return [random.uniform(1e-4, 0.13), random.uniform(1e-3, 0.5654), random.uniform(0, 1.6154), random.uniform(0, 1.5707963)]
    nu = 10 ** random.uniform(-6, 0.004)
    if ch in ('L2Tp', 'L2Tm'):
        return [nu, random.uniform(1e-3, 1.05), random.uniform(0, 1.05), random.uniform(-1.05, 1.05)]
    if ch == 'L2E':
        return [nu, random.uniform(0, 1.05), random.uniform(-1.05, 1.05), random.random()]
    if ch == 'L2P':
        return [nu, random.uniform(1e-3, 1.05), random.uniform(-1.05, 1.05), random.random()]
    return [nu, random.uniform(1e-3, 1.05), random.uniform(0, 1.05), random.random()]


def coeff_array(mod, names, A, nvar):
    import numpy as np
    out = np.empty((nvar + 1, 2, 2, 7), dtype=object)
    out[...] = mp.mpf(0)
    mod.coef_val([A[n] for n in names], out)
    return out


def chart_value(out, affine, x3, u, fi_single=0):
    def trig(fi):
        return sum(out[0, fi, 0, k] * mp.cos(k * u) + out[0, fi, 1, k] * mp.sin(k * u) for k in range(7))
    if affine:
        return (1 - x3) * trig(0) - x3 * trig(1)
    return trig(0)


def abs_textual(src):
    """check coef_abs terms == |coef_val terms| line by line"""
    body_v = src.split('def coef_val')[1].split('def coef_abs')[0]
    body_a = src.split('def coef_abs')[1]
    lv = dict(re.findall(r'out\[([^\]]+)\] = (.*)', body_v))
    la = dict(re.findall(r'out\[([^\]]+)\] = (.*)', body_a))
    bad = 0; nterm = 0
    assert lv.keys() == la.keys()
    for k in lv:
        tv = lv[k].split(' + '); ta = la[k].split(' + ')
        if len(tv) != len(ta):
            bad += 1; continue
        for a_, b_ in zip(tv, ta):
            nterm += 1
            ca, ma = a_.split(' * ', 1) if ' * ' in a_ else (a_, '')
            cb, mb = b_.split(' * ', 1) if ' * ' in b_ else (b_, '')
            if ma != mb or abs(float(ca)) != float(cb):
                bad += 1
    return bad, nterm


for ch in CHARTS:
    mod, src = load_plain(os.path.join(CODE, f'nbcoef_{ch}.py'))
    names = mod.ATOMS; nvar = mod.NVAR
    worst_v = mp.mpf(0); worst_g = mp.mpf(0)
    for it in range(NPTS):
        x = sample(ch)
        e, p, q, t, norm, A, affine = chart_map(ch, x)
        out = coeff_array(mod, names, A, nvar)
        x3 = mp.mpf(x[3]) if affine else None
        us = [mp.mpf(2) * mp.pi * (k + mp.mpf(1) / 3) / 8 for k in range(8)]
        ref = [f_phys(e, p, q, t, u) / norm for u in us]
        val = [chart_value(out, affine, x3, u) for u in us]
        sc = max(abs(r) for r in ref) + mp.mpf(10) ** -40
        worst_v = max(worst_v, max(abs(a - b) for a, b in zip(val, ref)) / sc)
        # gradients: d/dx_v of the value coefficients through atoms(x)
        if it < 3:
            for v in range(nvar):
                def coef_of(h, v=v):
                    xx = list(x); xx[v] = mp.mpf(xx[v]) + h
                    _, _, _, _, _, AA, _ = chart_map(ch, xx)
                    return coeff_array(mod, names, AA, nvar)
                hstep = mp.mpf(10) ** -18
                op, om = coef_of(hstep), coef_of(-hstep)
                for fi in range(2 if affine else 1):
                    for part in range(2):
                        for k in range(7):
                            num = (op[0, fi, part, k] - om[0, fi, part, k]) / (2 * hstep)
                            g = out[1 + v, fi, part, k]
                            worst_g = max(worst_g, abs(num - g) / (1 + abs(g)))
    bad, nterm = abs_textual(src)
    print(f"{ch:5s}: value vs physics (50 digits) max rel err {mp.nstr(worst_v, 3)}; gradient max rel err {mp.nstr(worst_g, 3)};"
          f" abs-poly terms checked {nterm}, mismatches {bad}", flush=True)

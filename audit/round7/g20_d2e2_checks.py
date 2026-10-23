"""AUDIT-7 g20 (round-7 addendum B): the new chart D2E2 and the D2Tp coverage rules of check_all.
  1. nbcoef_D2E2 (gen_bu.py): chart value vs MY physics at 50 digits (C1p function (1-tau) Phi / e^2 - tau (Phi + Psi)
     divided by the positive D2 normalisation, as for D2E in AUDIT_5 e4) with T = 1.04 + 0.24 x3; abs-poly terms;
     D2E2 and D2E describe the same function where their T ranges overlap ([1.04, 1.05]).
  2. Rule (ii): every point of every owed D2Tp box that check_all excuses (image in D2E for T <= 1.05, in D2E2 for
     T >= 1.04) - sampled at corners and random points with E > 0 - maps to the same physical point in D2E / D2E2
     coordinates, lies in a certified leaf of the D2E or D2E2 run (found by replaying their trees), and the 100-digit
     count there is <= that leaf's bound (winding leaves: #{|z| < e^L} == W).
  3. Rule (i): every owed D2 box dropped by the e-rule has smallest physical e > the largest e of C1p's owed boxes
     (exact arithmetic).
usage: python g20_d2e2_checks.py NPTS SEED"""
import os, re, sys, types, random, json, math, itertools
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.normpath(os.path.join(HERE, '..', '..', 'code'))
NPTS, SEED = int(sys.argv[1]), int(sys.argv[2])
random.seed(SEED)
mp.mp.dps = 50
SD = mp.mpf(5) / 1000


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


def Tmap(ch, x3):
    return (mp.mpf('1.04') + mp.mpf('0.24') * x3) if ch == 'D2E2' else (-mp.mpf('1.05') + mp.mpf('2.1') * x3)


# ---- 1. generated D2E2 coefficients vs physics ----
for ch in ('D2E2', 'D2E'):
    mod, src = load_plain(os.path.join(CODE, f'nbcoef_{ch}.py'))
    names = mod.ATOMS; nvar = mod.NVAR
    worst = mp.mpf(0)
    for it in range(12):
        nu = mp.mpf(10) ** mp.mpf(random.uniform(-4, math.log10(0.14)))
        r1, r2 = mp.mpf(random.uniform(0, 0.35)), mp.mpf(random.uniform(0.24, 0.68))
        x3 = mp.mpf(random.random()); T = Tmap(ch, x3)
        ee, pp, kk = SD * nu, SD * r1 * nu, SD * r2 * nu
        tauv = mp.mpf(1) / 2 + SD * T * nu; tt = ee * ee * tauv / (1 - tauv)
        A_ = {'nu': nu, 'b': mp.sqrt(1 - ee * ee), 'P': r1, 'K': r2}; A_['ib'] = 1 / A_['b']
        out = np.empty((nvar + 1, 2, 2, 7), dtype=object); out[...] = mp.mpf(0)
        mod.coef_val([A_[n] for n in names], out)
        norm = ee * ee * nu / (1 - tauv)
        for jj in range(5):
            uu = 2 * mp.pi * (jj + mp.mpf(1) / 3) / 5
            trig = lambda fi: sum(out[0, fi, 0, k] * mp.cos(k * uu) + out[0, fi, 1, k] * mp.sin(k * uu) for k in range(7))
            val = (1 - x3) * trig(0) - x3 * trig(1)
            ref = fphys(ee, pp, 1 - kk, tt, uu) / norm
            worst = max(worst, abs(val - ref) / (1 + abs(ref)))
    bv = src.split('def coef_val')[1].split('def coef_abs')[0]; ba = src.split('def coef_abs')[1]
    lv = dict(re.findall(r'out\[([^\]]+)\] = (.*)', bv)); la = dict(re.findall(r'out\[([^\]]+)\] = (.*)', ba))
    badt = nt = 0
    for key in lv:
        for a_, b_ in zip(lv[key].split(' + '), la[key].split(' + ')):
            nt += 1
            ca, ma = a_.split(' * ', 1) if ' * ' in a_ else (a_, ''); cb, mb = b_.split(' * ', 1) if ' * ' in b_ else (b_, '')
            badt += (ma != mb) or abs(float(ca)) != float(cb)
    print(f'1. {ch}: chart value vs physics (50 digits, nu <= 0.14, P <= 0.35, K in [0.24, 0.68]) max rel err '
          f'{mp.nstr(worst, 3)}; abs-poly terms {nt}, mismatches {badt}', flush=True)


# ---- helpers: leaves of a bb_run file ----
def leaves(fn):
    D = json.load(open(fn)); H = D['header']; S = H['scale']; out = []
    if H.get('split_rule'):                   # AUDIT-12: g20 replays the plain bb_run bisection only
        print('   PROBLEM: header split_rule %r not supported by g20 (leaf boxes would be wrong)' % H['split_rule'], flush=True)
    for r in D['results']:
        stack = [tuple(r['box'])]; wl = r.get('wL', []); wi = 0
        for c in r['tree']:
            bx = stack.pop()
            if c == '+':
                w = [(bx[2 * i + 1] - bx[2 * i]) / S[i] for i in range(4)]
                i = max(range(4), key=lambda j: w[j]); mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
                a = list(bx); b = list(bx); a[2 * i + 1] = mid; b[2 * i] = mid
                stack += [tuple(a), tuple(b)]
            elif '1' <= c <= '6':
                out.append((bx, 12 - 2 * int(c), int(c), wl[wi])); wi += 1
            elif 'a' <= c <= 'l' or 'A' <= c <= 'L':
                out.append((bx, ord(c) - (ord('a') if c >= 'a' else ord('A')), None, None))
    return H, D, out


def find(lv, x):
    for bx, b, W, L in lv:
        if all(bx[2 * i] <= x[i] <= bx[2 * i + 1] for i in range(4)):
            return bx, b, W, L
    return None


def roots_phys(ee, pp, qq, tt):
    N = 32
    vals = [fphys(ee, pp, qq, tt, 2 * mp.pi * n / N) for n in range(N)]
    c = [mp.fsum(vals[n] * mp.expj(-k * 2 * mp.pi * n / N) for n in range(N)) / N for k in range(-6, 7)]
    sc = max(abs(v) for v in c); lo = 0
    while abs(c[lo]) < mp.mpf(10) ** -40 * sc:
        lo += 1
    hi = 12
    while abs(c[hi]) < mp.mpf(10) ** -40 * sc:
        hi -= 1
    return list(mp.polyroots(c[lo:hi + 1][::-1], maxsteps=600, extraprec=400)) + [mp.mpf(0)] * lo


# ---- 2. rule (ii) on the owed D2Tp boxes ----
Htp, Dtp, _ = leaves(os.path.join(CODE, 'runs', 'D2Tp.json'))
owed = [tuple(f[0]) for r in Dtp['results'] for f in r['fails']] + [tuple(b) for r in Dtp['results'] for b in r['left']]
C = json.load(open(os.path.join(CODE, 'runs', 'C1p.json')))
c1p_owed = [f[0] for r in C['results'] for f in r['fails']] + [b for r in C['results'] for b in r['left']]
emax = max(mp.mpf(b[1]) for b in c1p_owed)
rule1 = [b for b in owed if SD * mp.mpf(b[0]) * mp.mpf(b[2]) > emax]
rest = [b for b in owed if b not in set(rule1)]
print(f'3. D2Tp owed boxes {len(owed)}; dropped by rule (i) (exact: SD*nu_lo*E_lo > {mp.nstr(emax, 6)}) {len(rule1)}; '
      f'remaining {len(rest)}: {rest}', flush=True)
_, _, LE = leaves(os.path.join(CODE, 'runs', 'D2E.json'))
_, _, LE2 = leaves(os.path.join(CODE, 'runs', 'D2E2.json'))
nck = nmiss = nviol = 0; maxN = 0
for b in rest:
    pts = [tuple(b[2 * i + ((m >> i) & 1)] for i in range(4)) for m in range(16)]
    pts += [tuple(random.uniform(b[2 * i], b[2 * i + 1]) for i in range(4)) for _ in range(NPTS)]
    for (n, E, P, K) in pts:
        if n <= 0 or E <= 0:
            continue
        n, E, P, K = [mp.mpf(v) for v in (n, E, P, K)]
        # D2Tp physical point
        ee, pp, kk, tp = SD * E * n, SD * P * n, SD * K * n, SD * n
        nuE, PE, KE, T = E * n, P / E, K / E, 1 / E
        assert abs(SD * nuE - ee) < mp.mpf(10) ** -45 and abs(SD * PE * nuE - pp) < mp.mpf(10) ** -45 and \
            abs(SD * KE * nuE - kk) < mp.mpf(10) ** -45 and abs(SD * T * nuE - tp) < mp.mpf(10) ** -45
        cands = []
        if T <= mp.mpf('1.05'):
            cands.append(('D2E', LE, (T + mp.mpf('1.05')) / mp.mpf('2.1')))
        if T >= mp.mpf('1.04'):
            cands.append(('D2E2', LE2, (T - mp.mpf('1.04')) / mp.mpf('0.24')))
        for nm, lv, x3 in cands:
            x = [float(nuE), float(PE), float(KE), float(x3)]
            hit = find(lv, x)
            nck += 1
            if hit is None:
                nmiss += 1; print('   NOT IN ANY LEAF', nm, x, flush=True); continue
            bx, bound, W, L = hit
            tau = mp.mpf(1) / 2 + tp; tt = ee * ee * tau / (1 - tau)
            z = roots_phys(ee, pp, 1 - kk, tt)
            Nm = sum(1 for r in z if r != 0 and abs(mp.log(abs(r))) < mp.mpf(10) ** -25)
            maxN = max(maxN, Nm)
            if Nm > min(bound, 10) or (W is not None and sum(1 for r in z if abs(r) < mp.e ** mp.mpf(L)) != W):
                nviol += 1; print('   VIOLATION', nm, x, Nm, bound, W, flush=True)
print(f'2. rule (ii): {nck} (point, chart) checks over the remaining D2Tp boxes: not in a certified leaf {nmiss}; '
      f'100-digit violations {nviol}; max N_mult {maxN}', flush=True)

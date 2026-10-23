"""AUDIT-7 g15 (PRELIMINARY, unrequested): checks of the second-order f-chart enclosure nb_cert2 (code/nb_cert2.py,
code/nbcoef_f2.py), wired into engine KQ (away count and winding) at 08:40 during this round.
  (A) coef_val2 slots 0..3 == nbcoef_f.coef_val slots 0..3 (same polynomials; nbcoef_f was checked against physics
      in AUDIT_1), float evaluation at random points;
  (B) coef_val2 slots 4..9 (second derivatives in e, p, q) == numerical derivatives (mpmath, 50 digits) of slots 1..3;
  (C) coef_abs2: slots 0..9 are the term-wise absolute values of coef_val2 (textual); slots 10..19 equal the
      absolute-coefficient polynomials of the symbolic third derivatives of coef_val2 slot 0 (sympy), at random points;
      term counts and degrees (float error budget of FLOAT_RIGOUR §1);
  (D) the Taylor inequality behind enclose_point2, per coefficient: for random boxes (e in [0.9, 1], all l, w) and random
      points x of the box, |c(x) - cmid - grad.d - 1/2 d'Hd| <= crad + rem (c(x) exact at 50 digits from the plain
      polynomials; d = x - x_c in (e, l, w); cmid, dcen, hcen, rem, crad from nb_cert2.box_coeffs2), and the observed
      ratio |remainder| / rem.
usage: python g15_nbcert2_checks.py NPTS NBOX SEED"""
import os, re, sys, types, random, math
import numpy as np
import mpmath as mp
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.normpath(os.path.join(HERE, '..', '..', 'code'))
sys.path.insert(0, CODE)
NPTS, NBOX, SEED = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
random.seed(SEED)
mp.mp.dps = 50


def load_plain(path):
    src = open(path).read()
    src = src.replace('from numba import njit', '')
    src = re.sub(r'@njit\([^)]*\)\n', '', src)
    mod = types.ModuleType(os.path.basename(path))
    exec(compile(src, path, 'exec'), mod.__dict__)
    return mod, src


F1, src1 = load_plain(os.path.join(CODE, 'nbcoef_f.py'))
F2, src2 = load_plain(os.path.join(CODE, 'nbcoef_f2.py'))


def ev(fn, nslot, e, p, q, obj=False):
    out = np.zeros((nslot, 2, 2, 7), dtype=object if obj else float)
    if obj:
        out[...] = mp.mpf(0)
    fn(e, p, q, out)
    return out


# (A)
worstA = 0.0
for _ in range(NPTS):
    e, p, q = random.uniform(0.0, 1.0), random.uniform(-1, 1), random.uniform(-1, 1)
    a = ev(F1.coef_val, 4, e, p, q); b = ev(F2.coef_val2, 10, e, p, q)
    sc = np.abs(a[:4]).max() + 1e-300
    worstA = max(worstA, float(np.abs(a[:4] - b[:4]).max() / sc))
print(f'(A) coef_val2[0..3] vs nbcoef_f.coef_val[0..3]: max rel diff {worstA:.2e}')

# (B)
worstB = mp.mpf(0)
pairs = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
for _ in range(max(3, NPTS // 10)):
    x0 = [mp.mpf(random.uniform(0.0, 1.0)), mp.mpf(random.uniform(-1, 1)), mp.mpf(random.uniform(-1, 1))]
    b = ev(F2.coef_val2, 10, *x0, obj=True)
    for s, (i, j) in enumerate(pairs):
        h = mp.mpf(10) ** -20
        xp = list(x0); xm = list(x0); xp[j] += h; xm[j] -= h
        dp_ = ev(F2.coef_val2, 10, *xp, obj=True); dm_ = ev(F2.coef_val2, 10, *xm, obj=True)
        for f in range(2):
            for pt in range(2):
                for k in range(7):
                    num = (dp_[1 + i, f, pt, k] - dm_[1 + i, f, pt, k]) / (2 * h)
                    worstB = max(worstB, abs(num - b[4 + s, f, pt, k]) / (1 + abs(b[4 + s, f, pt, k])))
print(f'(B) second derivatives vs 50-digit numerical derivatives: max rel err {mp.nstr(worstB, 3)}')

# (C) textual abs check (slots 0..9) and third derivatives (slots 10..19)
bv = dict(re.findall(r'out\[([^\]]+)\] = (.*)', src2.split('def coef_val2')[1].split('def coef_abs2')[0]))
ba = dict(re.findall(r'out\[([^\]]+)\] = (.*)', src2.split('def coef_abs2')[1]))
bad = 0; nterm = 0; maxterms = 0
for key, ex in bv.items():
    tv = ex.split(' + '); ta = ba[key].split(' + ')
    maxterms = max(maxterms, len(tv))
    if len(tv) != len(ta):
        bad += 1; continue
    for x_, y_ in zip(tv, ta):
        nterm += 1
        cx, mx = x_.split(' * ', 1) if ' * ' in x_ else (x_, '')
        cy, my = y_.split(' * ', 1) if ' * ' in y_ else (y_, '')
        if mx != my or abs(float(cx)) != float(cy):
            bad += 1
E, P, Q = sp.symbols('e p q')
loc = {f'e{d}': E ** d for d in range(12)}; loc.update({f'p{d}': P ** d for d in range(12)}); loc.update({f'q{d}': Q ** d for d in range(12)})
THI = [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 1), (0, 1, 2), (0, 2, 2), (1, 1, 1), (1, 1, 2), (1, 2, 2), (2, 2, 2)]
V = [E, P, Q]
worstC = 0.0; maxdeg = 0
pts = [(random.uniform(0, 1.2), random.uniform(0, 1.2), random.uniform(0, 1.2)) for _ in range(3)]
for f in range(2):
    for pt in range(2):
        for k in range(7):
            ex = sp.sympify(bv[f'0, {f}, {pt}, {k}'].replace('1.0 * ', ''), locals=loc) if bv[f'0, {f}, {pt}, {k}'] != '0.0' else sp.Integer(0)
            for m, (i, j, l) in enumerate(THI):
                d3 = sp.expand(sp.diff(ex, V[i], V[j], V[l]))
                absd3 = sum(abs(c) * E ** a * P ** b_ * Q ** g for (a, b_, g), c in sp.Poly(d3, E, P, Q).terms()) if d3 != 0 else 0
                if d3 != 0:
                    maxdeg = max(maxdeg, sp.Poly(d3, E, P, Q).total_degree())
                for (xe, xp_, xq) in pts:
                    got = ev(F2.coef_abs2, 20, xe, xp_, xq)[10 + m, f, pt, k]
                    want = float(absd3.subs({E: xe, P: xp_, Q: xq})) if absd3 != 0 else 0.0
                    worstC = max(worstC, abs(got - want) / (1 + abs(want)))
print(f'(C) coef_abs2 slots 0..9: {nterm} terms checked, mismatches {bad}; max terms per polynomial {maxterms}; '
      f'slots 10..19 vs symbolic third derivatives: max rel diff {worstC:.2e} (max degree {maxdeg})')

# (D) Taylor inequality per coefficient
import nb_cert2
worst_ratio = 0.0; viol = 0; ncheck = 0
for _ in range(NBOX):
    wd = [10 ** random.uniform(-4, -1) for _ in range(3)]
    e0 = random.uniform(0.9, 1.0 - wd[0]); l0 = random.uniform(0.0, 1.0 - wd[1]); w0 = random.uniform(0.0, math.pi / 2 - wd[2])
    bx = (e0, e0 + wd[0], l0, l0 + wd[1], w0, w0 + wd[2])
    cmid, crad, dcen, dwid, cmax, rr, hcen, rem = nb_cert2.box_coeffs2(*bx)
    xc = [0.5 * (bx[0] + bx[1]), 0.5 * (bx[2] + bx[3]), 0.5 * (bx[4] + bx[5])]
    for _ in range(4):
        x = [mp.mpf(random.uniform(bx[2 * i], bx[2 * i + 1])) for i in range(3)]
        d = [x[i] - mp.mpf(xc[i]) for i in range(3)]
        ce = ev(F1.coef_val, 4, x[0], x[1] * mp.sin(x[2]), x[1] * mp.cos(x[2]), obj=True)[0]
        for f in range(2):
            for pt in range(2):
                for k in range(7):
                    lin = sum(mp.mpf(float(dcen[v, f, pt, k])) * d[v] for v in range(3))
                    Hm = [hcen[m, f, pt, k] for m in range(6)]
                    quad = (Hm[0] * d[0] ** 2 + Hm[3] * d[1] ** 2 + Hm[5] * d[2] ** 2) / 2 + Hm[1] * d[0] * d[1] + \
                        Hm[2] * d[0] * d[2] + Hm[4] * d[1] * d[2]
                    r = abs(ce[f, pt, k] - mp.mpf(float(cmid[f, pt, k])) - lin - quad)
                    allowed = crad[f, pt, k] + rem[f, pt, k]
                    ncheck += 1
                    if r > allowed:
                        viol += 1
                        if viol <= 5:
                            print('   VIOLATION (D)', bx, [float(v) for v in x], f, pt, k, float(r), allowed)
                    if rem[f, pt, k] > 0:
                        worst_ratio = max(worst_ratio, float(r / allowed))
print(f'(D) Taylor inequality: {ncheck} coefficient checks on {NBOX} boxes x 4 points; violations {viol}; '
      f'max |c(x) - second-order model| / (crad + rem) = {worst_ratio:.3f}')

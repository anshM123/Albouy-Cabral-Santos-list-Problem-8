"""AUDIT-7 g17 (round-7 addendum A3): the second-order C3 / KQ polynomial engines (c3/gen_c3_2.py ->
nbcoef_c3_<C>_2.py, c3/make_c3_2.py -> nbeng_c3_<C>_2.py), charts KBQe (engine KQ, near part) and Ca.
  (A) coef_val2 slots 0..4 (value, first derivatives) == the first-order module nbcoef_c3_<C>.coef_val slots 0..4
      (checked against physics in AUDIT_3 c5);
  (B) slots 5..14 (second derivatives, sorted pairs) == 50-digit numerical derivatives of slots 1..4;
  (C) coef_abs2: slots 0..14 are the term-wise absolute values of coef_val2 (textual); slots 15..34 (sorted triples)
      equal the absolute-coefficient polynomials of the symbolic third derivatives of slot 0 (sympy), at random points;
  (D) per-coefficient Taylor inequality |c_k(x) - cmid_k - dcen_k.d - 1/2 d'H_k d| <= crad_k + rem_k at 50 digits
      (box_coeffs2 data; random boxes of the chart domain, widths 1e-5 .. 0.3 of the domain; random points);
  (E) point_enclosure2 containment, modes 1 (second order), 2 (mean value), 0 (min, used): the exact
      G^(j)(y; x) = sum_k c_k(x) (k)_j y^(k-j) at random y of the root interval and random x of the box (corners
      included) must lie in [vc - rad, vc + rad], j = 0..7.
usage: python g17_c3_second_order.py NBOX SEED CHART [CHART ...]"""
import os, re, sys, types, random, math, itertools
import numpy as np
import mpmath as mp
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
C3 = os.path.normpath(os.path.join(HERE, '..', '..', 'code', 'c3'))
sys.path.insert(0, C3); sys.path.insert(0, os.path.dirname(C3))
NBOX, SEED = int(sys.argv[1]), int(sys.argv[2]); charts = sys.argv[3:]
random.seed(SEED)
mp.mp.dps = 50
J = 6
NU1, R1 = 0.37, 1.001
DOMS = {'Ca': ((0, NU1), (0, R1 ** 2), (0, R1 ** 2), (0, R1 ** 3)), 'KBQe': ((0, 0.0501), (0, 1.0), (0, 1.0), (0, 1.0))}
ROOTS = {'Ca': (-3.0, 3.0), 'KBQe': (-0.4831, 0.4831)}
PAIRS = sorted(set(tuple(sorted(x)) for x in itertools.product(range(4), repeat=2)))
TRIPLES = sorted(set(tuple(sorted(x)) for x in itertools.product(range(4), repeat=3)))


def load_plain(path):
    src = open(path).read().replace('from numba import njit', '')
    src = re.sub(r'@njit\([^)]*\)\n', '', src)
    mod = types.ModuleType(os.path.basename(path)); exec(compile(src, path, 'exec'), mod.__dict__)
    return mod, src


def ev(fn, nslot, deg, x, obj=False):
    out = np.zeros((nslot, deg + 1), dtype=object if obj else float)
    if obj:
        out[...] = mp.mpf(0)
    fn(np.array(x, dtype=object) if obj else np.array(x, dtype=float), out)
    return out


for ch in charts:
    import importlib
    M1, _ = load_plain(os.path.join(C3, f'nbcoef_c3_{ch}.py'))
    M2, src2 = load_plain(os.path.join(C3, f'nbcoef_c3_{ch}_2.py'))
    eng = importlib.import_module(f'nbeng_c3_{ch}_2')
    DEG = M2.DEG; dom = DOMS[ch]; rlo, rhi = ROOTS[ch]
    assert M1.DEG == DEG
    # (A)
    wA = 0.0
    for _ in range(40):
        x = [random.uniform(a, b) for a, b in dom]
        o1 = ev(M1.coef_val, 5, DEG, x); o2 = ev(M2.coef_val2, 15, DEG, x)
        wA = max(wA, float(np.abs(o1 - o2[:5]).max() / (np.abs(o1).max() + 1e-300)))
    # (B)
    wB = mp.mpf(0)
    for _ in range(4):
        x0 = [mp.mpf(random.uniform(a, b)) for a, b in dom]
        o2 = ev(M2.coef_val2, 15, DEG, x0, obj=True)
        h = mp.mpf(10) ** -20
        for m, (a, b) in enumerate(PAIRS):
            xp = list(x0); xm = list(x0); xp[b] += h; xm[b] -= h
            num = (ev(M2.coef_val2, 15, DEG, xp, obj=True)[1 + a] - ev(M2.coef_val2, 15, DEG, xm, obj=True)[1 + a]) / (2 * h)
            for k in range(DEG + 1):
                wB = max(wB, abs(num[k] - o2[5 + m, k]) / (1 + abs(o2[5 + m, k])))
    # (C)
    body_v = src2.split('def coef_val2')[1].split('def coef_abs2')[0]; body_a = src2.split('def coef_abs2')[1]
    bv = dict(re.findall(r'out\[([^\]]+)\] = (.*)', body_v)); ba = dict(re.findall(r'out\[([^\]]+)\] = (.*)', body_a))
    bad = 0; nterm = 0
    for key, exv in bv.items():
        tv = exv.split(' + '); ta = ba[key].split(' + ')
        if len(tv) != len(ta):
            bad += 1; continue
        for x_, y_ in zip(tv, ta):
            nterm += 1
            cx, mx = x_.split(' * ', 1) if ' * ' in x_ else (x_, '')
            cy, my = y_.split(' * ', 1) if ' * ' in y_ else (y_, '')
            if mx != my or abs(float(cx)) != float(cy):
                bad += 1
    names = M2.BOX
    syms = sp.symbols(' '.join(names))
    loc = {}
    for n, s_ in zip(names, syms):
        for d in range(1, 30):
            loc[f'{n}{d}'] = s_ ** d
    wC = 0.0
    pts = [[random.uniform(0, 1.2 * b) for a, b in dom] for _ in range(2)]
    for k in range(DEG + 1):
        exs = bv[f'0, {k}']
        ex = sp.Integer(0) if exs == '0.0' else sp.sympify(exs.replace('1.0 * ', ''), locals=loc)
        for m, (a, b, c) in enumerate(TRIPLES):
            d3 = sp.expand(sp.diff(ex, syms[a], syms[b], syms[c]))
            for p in pts:
                got = ev(M2.coef_abs2, 35, DEG, p)[15 + m, k]
                if d3 == 0:
                    want = 0.0
                else:
                    want = float(sum(abs(cf) * sp.Mul(*[sp.Float(p[i]) ** mon[i] for i in range(4)])
                                     for mon, cf in sp.Poly(d3, *syms).terms()))
                wC = max(wC, abs(got - want) / (1 + abs(want)))
    print(f'{ch}: (A) slots 0..4 vs first-order module: max rel diff {wA:.2e}; (B) second derivatives vs numerical: '
          f'max rel err {mp.nstr(wB, 3)}; (C) abs-poly terms {nterm} (mismatches {bad}), third derivatives vs sympy: '
          f'max rel diff {wC:.2e}', flush=True)
    # (D), (E)
    nD = vD = 0; rD = 0.0; slackD = []; slackE = {}
    st = {0: [0, 0, 0.0], 1: [0, 0, 0.0], 2: [0, 0, 0.0]}
    for _ in range(NBOX):
        bx = []
        for a, b in dom:
            wd = (b - a) * 10 ** random.uniform(-5, math.log10(0.3))
            lo = random.uniform(a, b - wd) if random.random() > 0.15 else a
            bx += [lo, lo + wd]
        data = eng.box_coeffs2(*bx)
        cmid, crad, dcen, dwid, cmax, rr, hcen, rem, qw = data
        xc = [mp.mpf(0.5 * (bx[2 * i] + bx[2 * i + 1])) for i in range(4)]
        pts = [[bx[2 * i + ((m >> i) & 1)] for i in range(4)] for m in random.sample(range(16), 3)]
        pts += [[random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(4)] for _ in range(3)]
        exact = []
        for p in pts:
            cx = ev(M2.coef_val2, 15, DEG, [mp.mpf(v) for v in p], obj=True)[0]
            exact.append(cx)
            d = [mp.mpf(p[i]) - xc[i] for i in range(4)]
            for k in range(DEG + 1):
                lin = sum(mp.mpf(float(dcen[v, k])) * d[v] for v in range(4))
                quad = mp.mpf(0)
                for m, (a, b) in enumerate(PAIRS):
                    hm = mp.mpf(float(hcen[m, k]))
                    quad += (hm * d[a] * d[a] / 2) if a == b else (hm * d[a] * d[b])
                r = abs(cx[k] - mp.mpf(float(cmid[k])) - lin - quad)
                nD += 1
                allowed = crad[k] + rem[k]
                if r > allowed:
                    vD += 1
                    if vD <= 3:
                        print('   VIOLATION (D)', ch, bx, p, k, float(r), allowed, flush=True)
                if allowed > 0:
                    rD = max(rD, float(r / allowed))
                    if float(r / allowed) > 0.9999:
                        slackD.append(float(1 - r / allowed))
        for _y in range(2):
            y = random.uniform(rlo, rhi)
            for mode in (0, 1, 2):
                vc, rad = eng.point_enclosure2(cmid, crad, dcen, dwid, hcen, rem, qw, rr, y, J, mode)
                for j in range(J + 2):
                    for cx in exact:
                        exv = mp.mpf(0)
                        for k in range(j, DEG + 1):
                            ffk = mp.mpf(1)
                            for i in range(j):
                                ffk *= (k - i)
                            exv += cx[k] * ffk * mp.mpf(y) ** (k - j)
                        dd = abs(exv - mp.mpf(float(vc[j])))
                        s_ = st[mode]; s_[0] += 1
                        if dd > rad[j]:
                            s_[1] += 1
                            if s_[1] <= 3:
                                print('   VIOLATION (E)', ch, 'mode', mode, bx, y, j, float(dd), rad[j], flush=True)
                        if rad[j] > 0:
                            s_[2] = max(s_[2], float(dd / rad[j]))
                            if float(dd / rad[j]) > 0.9999:
                                slackE.setdefault(mode, []).append(float(1 - dd / rad[j]))
    print(f'{ch}: (D) Taylor inequality: {nD} checks, violations {vD}, worst |remainder| / (crad + rem) {rD:.3f}; '
          f'checks within 1e-4 of the bound {len(slackD)}, smallest relative slack {min(slackD) if slackD else float("nan"):.2e}', flush=True)
    for mode in sorted(slackE):
        print(f'{ch}: (E) mode {mode}: checks within 1e-4 of the radius {len(slackE[mode])}, smallest relative slack {min(slackE[mode]):.2e}', flush=True)
    for mode, nm in ((1, 'second-order only'), (2, 'mean-value only'), (0, 'min (used)')):
        s_ = st[mode]
        print(f'{ch}: (E) point_enclosure2 mode {mode} ({nm}): {s_[0]} checks, violations {s_[1]}, worst ratio {s_[2]:.6f}', flush=True)

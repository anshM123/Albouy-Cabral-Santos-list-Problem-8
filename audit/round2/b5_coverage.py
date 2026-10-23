"""AUDIT-2 b5: coverage of the physical parameter space by the production charts, using the ACTUAL run domains,
grids and exclusions (command lines of the running bb_run jobs, 2026-10-23 00:59):

  bulk  a-chart  eta in [0.20033740617957885, 1.2535563099909826]
  C1    eps in [0, 0.26], th_C1 in [0, 1.1071487177940904], grid 8x12x4x16, skip boxes inside
        [0, 0.0650001] x [0, 0.1846] x [0, 1.5708] x [0.375, 0.625]
  C1p   e in [0, 0.2001], l in [0.4, 1];  C1pp  e/l^2 in [0, 0.51], l in [0, 0.52]
  B1    e in [0, 0.07], ph, qh in [0, 1.75], tauh in [-2, 2], grid 7x7x7x8, skip boxes inside
        [0, 0.0100001] x [0, 0.2500001] x [0.7499999, 1.2500001] x [0.3749999, 0.6250001] (x3 units)
  B2    theta in [0, 0.2], e/theta <= 0.3675, |tau|/theta <= 0.6825
  B3+-  eta = |tau| in [0, 0.13], e/eta <= 0.5654, theta/eta <= 1.6154
  L2*   X = (e/0.01, ph/0.25, (qh-1)/0.25, tauh/0.5), nu = dominant |X_i| <= 1.01, other ratios in [-1.05, 1.05]
  F1    f-chart e in [0.95, 1], t_f <= 0.9975;   F2  e in [0.9499, 0.99], t_f >= 0.9975
  K     e in [0.99, 1], mu = q/R <= 0.0505, skip l >= 0.95 & w >= 1.4137 (C3 core, pending)
Physical space: e in (0, 1), l in [0, 1], w in [0, pi/2], t in [0, 1] (a-chart t = S/(1+S)).
Reports uncovered samples, separately for the declared C3 core.
"""
import math
import numpy as np

rng = np.random.default_rng(55)
ETA = (0.20033740617957885, 1.2535563099909826)


def grid_edges(lo, hi, n):
    return np.linspace(lo, hi, n + 1)


def skipped_boxes(dom, n, ex):
    E = [grid_edges(dom[2 * i], dom[2 * i + 1], n[i]) for i in range(4)]
    out = []
    for a in range(n[0]):
        for b in range(n[1]):
            for c in range(n[2]):
                for d in range(n[3]):
                    bx = (E[0][a], E[0][a + 1], E[1][b], E[1][b + 1], E[2][c], E[2][c + 1], E[3][d], E[3][d + 1])
                    if all(ex[2 * i] <= bx[2 * i] and bx[2 * i + 1] <= ex[2 * i + 1] for i in range(4)):
                        out.append(bx)
    return out


C1_SKIP = skipped_boxes((0, 0.26, 0, 1.1071487177940904, 0, 1.5707963267948966, 0, 1), (8, 12, 4, 16),
                        (0, 0.0650001, 0, 0.1846, 0, 1.5708, 0.375, 0.625))
B1_SKIP = skipped_boxes((0, 0.07, 0, 1.75, 0, 1.75, 0, 1), (7, 7, 7, 8),
                        (0, 0.0100001, 0, 0.2500001, 0.7499999, 1.2500001, 0.3749999, 0.6250001))
K_SKIP = skipped_boxes((0.99, 1.0, 0, 1, 0, 1.5707963267948966, 0, 0.0505), (2, 20, 10, 2),
                       (0.99, 1.0, 0.9499999, 1.0, 1.4137, 1.5708, 0, 0.0505))


def in_any(pt, boxes):
    return any(all(b[2 * i] <= pt[i] <= b[2 * i + 1] for i in range(4)) for b in boxes)


def covered(e, l, w, t):
    """returns list of chart names covering the physical point"""
    ch = []
    if ETA[0] <= math.asin(e) <= ETA[1]:
        ch.append('bulk')
    eps = math.hypot(e, l * l); thc = math.atan2(l * l, e)
    if eps <= 0.26 and thc <= 1.1071487177940904 and not in_any((eps, thc, w, t), C1_SKIP):
        ch.append('C1')
    if e <= 0.2001 and 0.4 <= l <= 1.0:
        ch.append('C1p')
    if l > 0 and e / (l * l) <= 0.51 and l <= 0.52:
        ch.append('C1pp')
    # near-coincident cube charts (only relevant inside the C1 skip block, but test everywhere)
    tau = t - 0.5
    theta = l * l / e
    if e <= 0.07:
        ph = l * math.sin(w) / e; qh = l * math.cos(w) / e; tauh = tau / e
        x3 = (tauh + 2) / 4
        if ph <= 1.75 and qh <= 1.75 and -2 <= tauh <= 2 and not in_any((e, ph, qh, x3), B1_SKIP):
            ch.append('B1')
        X = (e / 0.01, ph / 0.25, (qh - 1) / 0.25, tauh / 0.5)
        nu = max(abs(x) for x in X)
        if nu <= 1.01 and nu > 0 and all(abs(x) / nu <= 1.05 for x in X):
            ch.append('L2')
    if theta <= 0.2 and e / theta <= 0.3675 and abs(tau) / theta <= 0.6825:
        ch.append('B2')
    if 0 < abs(tau) <= 0.13 and e / abs(tau) <= 0.5654 and theta / abs(tau) <= 1.6154:
        ch.append('B3')
    # e -> 1 charts
    if e >= 0.9499:
        S = t / (1 - t) if t < 1 else math.inf
        kap = math.sqrt(S) / (1 - e)
        tf = kap * kap / (1 + kap * kap) if math.isfinite(kap) else 1.0
        if 0.95 <= e <= 1 and tf <= 0.9975:
            ch.append('F1')
        if 0.9499 <= e <= 0.99 and tf >= 0.9975:
            ch.append('F2')
        mu = 1 / kap if kap > 0 else math.inf
        if 0.99 <= e <= 1 and mu <= 0.0505:
            if in_any((e, l, w, mu), K_SKIP):
                ch.append('C3core')
            else:
                ch.append('K')
    return ch


def draw(kind):
    if kind == 'uniform':
        return rng.uniform(1e-9, 1 - 1e-9), rng.uniform(0, 1), rng.uniform(0, math.pi / 2), rng.uniform(0, 1)
    if kind == 'small_e':
        e = 10 ** rng.uniform(-7, math.log10(0.21)); l = rng.choice([rng.uniform(0, 1), 10 ** rng.uniform(-5, 0)])
        return e, l, rng.uniform(0, math.pi / 2), rng.choice([rng.uniform(0, 1), 0.5 + rng.uniform(-0.13, 0.13)])
    if kind == 'C1block':
        eps = rng.uniform(0, 0.0650001); thc = rng.uniform(0, 0.1846); t = rng.uniform(0.375, 0.625)
        return max(eps * math.cos(thc), 1e-12), math.sqrt(eps * math.sin(thc)), rng.uniform(0, math.pi / 2), t
    if kind == 'B1block':
        e = rng.uniform(1e-9, 0.0100001); ph = rng.uniform(0, 0.2500001); qh = rng.uniform(0.7499999, 1.2500001)
        tauh = rng.uniform(-0.5000004, 0.5000004)
        l = math.hypot(ph, qh) * e; w = math.atan2(ph, qh)
        return e, min(l, 1.0), w, 0.5 + e * tauh
    if kind == 'top':
        e = 1 - 10 ** rng.uniform(-9, math.log10(0.06)); t = rng.choice([rng.uniform(0, 1), 1 - 10 ** rng.uniform(-12, -1)])
        return e, rng.uniform(0, 1), rng.uniform(0, math.pi / 2), t
    if kind == 'seams':
        e = rng.choice([0.199, 0.2, 0.2001, 0.9499, 0.95, 0.9501, 0.99]) * (1 + rng.uniform(-2e-4, 2e-4))
        return e, rng.uniform(0, 1), rng.uniform(0, math.pi / 2), rng.uniform(0, 1)


for kind, n in (('uniform', 100000), ('small_e', 200000), ('C1block', 100000), ('B1block', 100000), ('top', 200000),
                ('seams', 100000)):
    unc = []; c3 = 0; hist = {}
    for _ in range(n):
        e, l, w, t = draw(kind)
        if not (0 < e < 1):
            continue
        ch = covered(e, l, w, t)
        for c in ch:
            hist[c] = hist.get(c, 0) + 1
        real = [c for c in ch if c != 'C3core']
        if not real:
            if 'C3core' in ch:
                c3 += 1
            else:
                unc.append((e, l, w, t))
    print(f"[{kind:8s}] {n} samples: uncovered {len(unc)} (+ {c3} in the declared C3 core); chart hits {dict(sorted(hist.items()))}")
    for x in unc[:5]:
        print("      uncovered e=%.10g l=%.6g w=%.4f t=%.6g" % x)
print("C1 skipped initial boxes:", len(C1_SKIP), "; B1 skipped:", len(B1_SKIP), "; K skipped:", len(K_SKIP))
print("C1 skip block max eps %.10f, max th_C1 %.10f -> max theta = tan = %.6f; t range [%.4f, %.4f]" % (
    max(b[1] for b in C1_SKIP), max(b[3] for b in C1_SKIP), math.tan(max(b[3] for b in C1_SKIP)),
    min(b[6] for b in C1_SKIP), max(b[7] for b in C1_SKIP)))
print("B1 skip block: e <= %.12g, ph <= %.12g, qh in [%.12g, %.12g], x3 in [%.12g, %.12g]" % (
    max(b[1] for b in B1_SKIP), max(b[3] for b in B1_SKIP), min(b[4] for b in B1_SKIP), max(b[5] for b in B1_SKIP),
    min(b[6] for b in B1_SKIP), max(b[7] for b in B1_SKIP)))

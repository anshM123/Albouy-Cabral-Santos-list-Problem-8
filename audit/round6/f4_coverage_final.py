"""AUDIT-6 f4: coverage of e >= 0.94 by the FINAL atlas (launch_runs.sh, 07:30): bulk (a-chart, e <= sin 1.2535563...),
F0 (f, e in [0.9499,1], t <= 0.5), KQ1 (e in [0.9499,1], t in [0.5,0.9975]), KQ2 (e in [0.9499,0.99], t >= 0.9975),
KQK1 (e >= 0.99, l <= 0.95, t >= 0.997456), KQK2 (e >= 0.99, l >= 0.95, w <= 3pi/20, t >= 0.997456),
KQB (e >= 0.99, l >= 0.95, 3pi/20 <= w <= 9pi/20, t >= 0.997456), C3 cell + KA (e >= 0.99, l >= 0.95, w >= 9pi/20,
mu = q/R <= 0.0505).  f-chart t = kappa^2/(1+kappa^2), kappa = s/(1-e) (s = R/a), mu = 1/kappa."""
import math, random
random.seed(606)
W3, W9, HP = 0.47123889803846897, 1.413716694115407, 1.5707963267948966
unc = []; hits = {}
for it in range(500000):
    e = 1 - 10 ** random.uniform(-11, math.log10(0.06)) if it % 2 else random.uniform(0.94, 1 - 1e-12)
    l = random.choice([random.uniform(0, 1), 1 - 10 ** random.uniform(-9, -1.2), 0.95 + random.uniform(-1e-7, 1e-7)])
    l = min(max(l, 0.0), 1.0)
    w = random.choice([random.uniform(0, HP), random.choice([W3, W9]) + random.uniform(-1e-7, 1e-7)])
    w = min(max(w, 0.0), HP)
    ta = random.choice([random.uniform(0, 1), 1 - 10 ** random.uniform(-14, -1)])
    S = ta / (1 - ta); s = math.sqrt(S); kap = s / (1 - e)
    tf = kap * kap / (1 + kap * kap) if kap < 1e150 else 1.0
    mu = 1 / kap if kap > 0 else math.inf
    tf = random.choice([tf, random.choice([0.5, 0.9975, 0.997456]) + random.uniform(-1e-9, 1e-9)])
    tf = min(max(tf, 0.0), 1.0)
    mu = math.sqrt((1 - tf) / tf) if tf > 0 else math.inf
    cov = []
    if math.asin(e) <= 1.2535563099909826: cov.append('bulk')
    if 0.9499 <= e <= 1 and tf <= 0.5: cov.append('F0')
    if 0.9499 <= e <= 1 and 0.5 <= tf <= 0.9975: cov.append('KQ1')
    if 0.9499 <= e <= 0.99 and tf >= 0.9975: cov.append('KQ2')
    if e >= 0.99 and tf >= 0.997456:
        if l <= 0.95: cov.append('KQK1')
        if l >= 0.95 and w <= W3: cov.append('KQK2')
        if l >= 0.95 and W3 <= w <= W9: cov.append('KQB')
    if e >= 0.99 and l >= 0.95 and w >= W9 and mu <= 0.0505: cov.append('C3+KA')
    for c in cov: hits[c] = hits.get(c, 0) + 1
    if not cov: unc.append((e, l, w, tf))
print(f"samples 500000 (e >= 0.94, seams 0.5/0.9975/0.997456, l = 0.95, w = 3pi/20, 9pi/20 stressed): uncovered {len(unc)}; hits {dict(sorted(hits.items()))}")
for u in unc[:5]: print('   uncovered', u)
print("1/(1 + 0.0505^2) =", 1 / (1 + 0.0505 ** 2), "(C3 cell t-threshold; KQK1/KQK2/KQB start at 0.997456 and KQ1 covers t <= 0.9975)")

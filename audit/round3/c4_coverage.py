"""AUDIT-3 c4: (a)+(b) C3 cone / zone membership for physical points of the C3 cell, and (f) coverage of e >= 0.95.

(a) For p = (e, l, w, mu) in the cell (e in [0.99, 1), l in [0.95, 1], w in [9pi/20, pi/2], mu in (0, 0.0505]) and tau:
    |tau| <= 3 delta(p): the dominant level-1 chart (ratios <= 1) must contain the ratio point (DOM of c3_run.py);
    3 delta <= |tau| <= 1.11: PhiP/PhiM with nu = |tau|/3: ratios a, K, E <= R1^k and root m <= R1^3;
    |tau| >= 1.11: away arc of KA (|f - pi| >= PHI_A) and p in the KA domain.
(b) zone: if the dominant chart is Ca/CK and E <= E2, m <= M2: inner roots |y| <= Y2 are in an L2*E/L2*m chart
    (|y| <= Y2 eta_p) or an L2*Y+- chart (Y2 eta_p <= |y| <= Y2): check membership of the corresponding chart points.
(f) e >= 0.95: (e, l, w, t) covered by bulk / F1 / F2 / K (minus its skipped initial boxes) / KB / C3 cell (KA domain),
    using the actual command-line parameters of the running jobs (01:33)."""
import math, random
import numpy as np
random.seed(404)
Y0, Y2, E2, M2, R1, R2, NU1, ETA1 = 3.0, 0.25, 0.01, 0.01, 1.001, 1.001, 0.37, 1.0001
DOM = {'Ca': (0, NU1, 0, R1**2, 0, R1**2, 0, R1**3), 'CK': (0, NU1, 0, R1, 0, R1**2, 0, R1**3),
       'CE': (0, NU1, 0, R1, 0, R1**2, 0, R1**3), 'Cm': (0, NU1, 0, R1, 0, R1**2, 0, R1**2),
       'Phi': (0, NU1, 0, R1, 0, R1**2, 0, R1**2), 'L2a': (0, NU1, 0, R1**2, 0, ETA1, 0, R2**2),
       'L2K': (0, NU1, 0, R1, 0, ETA1, 0, R2**2)}
W9 = 9 * math.pi / 20
PHI_A = 2 * math.atan(3.0 * 0.37) - 1e-9


def inside(pt, d):
    return all(d[2 * i] - 1e-15 <= pt[i] <= d[2 * i + 1] + 1e-15 for i in range(4))


bad = {'delta': 0, 'near': 0, 'inter': 0, 'away': 0, 'zone': 0}; n = 0; nz = 0; dmax = 0
for it in range(200000):
    e = 1 - 10 ** random.uniform(-9, -2); l = random.uniform(0.95, 1.0); w = random.uniform(W9, math.pi / 2)
    mu = 10 ** random.uniform(-9, math.log10(0.0505))
    if it % 3 == 0:
        mu = random.uniform(1e-9, 0.0505)
    qz, kap, eps = l * math.cos(w), 1 - l * math.sin(w), 1 - e
    d = max(qz, math.sqrt(kap), math.sqrt(eps), mu ** (1 / 3)); dmax = max(dmax, d)
    if d > 0.36963:
        bad['delta'] += 1
    comps = {'Ca': qz, 'CK': math.sqrt(kap), 'CE': math.sqrt(eps), 'Cm': mu ** (1 / 3)}
    ch = max(comps, key=comps.get); nu = d
    rat = {'Ca': (nu, kap / nu**2, eps / nu**2, mu / nu**3), 'CK': (nu, qz / nu, eps / nu**2, mu / nu**3),
           'CE': (nu, qz / nu, kap / nu**2, mu / nu**3), 'Cm': (nu, qz / nu, kap / nu**2, eps / nu**2)}[ch]
    n += 1
    if not inside(rat, DOM[ch]):
        bad['near'] += 1
    # intermediate: random tau in [3 delta, 1.11]
    tau = random.uniform(3 * d, 1.11)
    nu2 = tau / 3
    pt = (nu2, qz / nu2, kap / nu2**2, eps / nu2**2); m_root = mu / nu2**3
    if not (inside(pt, DOM['Phi']) and 0 <= m_root <= R1**3 + 1e-15):
        bad['inter'] += 1
    # away: tau >= 1.11 -> |f - pi| = 2 atan(tau) >= 2 atan(1.11) > PHI_A; p in KA domain
    if not (2 * math.atan(1.11) >= PHI_A and 0.99 <= e <= 1 and 0.95 <= l <= 1 and W9 <= w <= math.pi / 2 and mu <= 0.0505):
        bad['away'] += 1
    # zone
    if ch in ('Ca', 'CK'):
        E, m = rat[2], rat[3]
        if E <= E2 and m <= M2:
            nz += 1
            rv = rat[1]
            etap = max(math.sqrt(E / E2), math.sqrt(m / M2))
            pref = 'L2a' if ch == 'Ca' else 'L2K'
            ok = True
            if etap > 0:
                # E/m-dominant chart: (nu, rv, eta_p, ratio) with ratio = (m/M2)/eta^2 or (E/E2)/eta^2 <= 1, Y in [-1, 1]
                r4 = (m / M2) / etap**2 if E / E2 >= m / M2 else (E / E2) / etap**2
                ok &= inside((nu, rv, etap, r4), DOM[pref])
            for frac in (0.0, 0.3, 0.7, 1.0):           # |y| in [Y2 eta_p, Y2]: Y-dominant charts
                eta = etap + (1.0 - etap) * frac
                if eta <= 0:
                    continue
                ok &= inside((nu, rv, eta, (E / E2) / eta**2), DOM[pref]) and 0 <= (m / M2) / eta**2 <= R2**2 + 1e-15
            if not ok:
                bad['zone'] += 1
print(f"(a)/(b) C3 cell samples {n}: max delta {dmax:.6f} (NU1 = 0.37); membership failures {bad}; zone samples {nz}")

# (f) e >= 0.95 coverage with the running jobs' parameters
ETA_HI = 1.2535563099909826
gl = np.linspace(0, 1, 21); gw = np.linspace(0, math.pi / 2, 11); ge = np.linspace(0.99, 1.0, 3); gm = np.linspace(0, 0.0505, 3)
Kskip = [(ge[a], ge[a + 1], gl[b], gl[b + 1], gw[c], gw[c + 1], gm[d], gm[d + 1]) for a in range(2) for b in range(20)
         for c in range(10) for d in range(2)
         if 0.99 <= ge[a] and ge[a + 1] <= 1.0 and 0.9499999 <= gl[b] and gl[b + 1] <= 1.0 and 0.47 <= gw[c]
         and gw[c + 1] <= 1.5708]
print("K skipped initial boxes:", len(Kskip))


def in_K_skip(e, l, w, mu):
    return any(b[0] <= e <= b[1] and b[2] <= l <= b[3] and b[4] <= w <= b[5] and b[6] <= mu <= b[7] for b in Kskip)


unc = []; hits = {}
for it in range(400000):
    e = 1 - 10 ** random.uniform(-10, math.log10(0.06)) if it % 2 else random.uniform(0.94, 1 - 1e-12)
    l = random.choice([random.uniform(0, 1), 1 - 10 ** random.uniform(-9, -1.3)])
    w = random.choice([random.uniform(0, math.pi / 2), random.choice([3 * math.pi / 20, W9]) + random.uniform(-1e-6, 1e-6)])
    w = min(max(w, 0.0), math.pi / 2)
    t = random.choice([random.uniform(0, 1), 1 - 10 ** random.uniform(-12, -1)])
    S = t / (1 - t); s = math.sqrt(S); kap_f = s / (1 - e); tf = kap_f**2 / (1 + kap_f**2); mu = 1 / kap_f if kap_f > 0 else math.inf
    cov = []
    if math.asin(e) <= ETA_HI: cov.append('bulk')
    if 0.95 <= e <= 1 and tf <= 0.9975: cov.append('F1')
    if 0.9499 <= e <= 0.99 and tf >= 0.9975: cov.append('F2')
    if 0.99 <= e <= 1 and mu <= 0.0505:
        if not in_K_skip(e, l, w, mu): cov.append('K')
        if 0.95 <= l <= 1 and 3 * math.pi / 20 <= w <= W9: cov.append('KB')
        if 0.95 <= l <= 1 and W9 <= w <= math.pi / 2: cov.append('C3+KA')
    for c in cov: hits[c] = hits.get(c, 0) + 1
    if not cov: unc.append((e, l, w, t))
print(f"(f) e >= 0.94 samples 400000: uncovered {len(unc)}; hits {dict(sorted(hits.items()))}")
for u in unc[:5]: print("   uncovered", u)

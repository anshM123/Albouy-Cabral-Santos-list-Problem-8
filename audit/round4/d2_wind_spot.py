"""AUDIT-4 d2: 100-digit spot checks of winding leaves ('1'..'6' codes) with MY OWN physics and chart maps.

For a sampled winding leaf with recorded W and log-radius L (record['wL'], in tree order), at random box points:
  physical (e, p, q, t) via my own chart map  ->  f_t(u) = (1-t) A^2 rho^2 - t ((rho^2)')^2  (a-chart closed forms)
  ->  Fourier coefficients c_k by a 32-point DFT at 100 digits  ->  Q(z) = sum c_k z^{k+6}  ->  all roots (polyroots).
Checks:  #{|z| < e^L} == W  (the certificate's claim, uniform over the box);  no root within 1e-30 of |z| = e^L;
         #{|z| = 1} <= 12 - 2W.   The chart function is a positive multiple of f_t, so the roots are the same.
usage: python d2_wind_spot.py NLEAF NPT SEED RUN.json [...]
"""
import sys, json, math, random
import numpy as np
import mpmath as mp
mp.mp.dps = 100
NLEAF = int(sys.argv[1]); NPT = int(sys.argv[2]); random.seed(int(sys.argv[3])); files = sys.argv[4:]


def split(bx, S):
    w = [(bx[2 * i + 1] - bx[2 * i]) / S[i] for i in range(4)]
    i = int(np.argmax(w)); mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)


def phys(eng, x):
    x = [mp.mpf(v) for v in x]
    if eng == 'B1':
        e, ph, qh, x3 = x; tauh = -2 + 4 * x3
        return e, e * ph, e * qh, mp.mpf(1) / 2 + e * tauh
    if eng == 'L2E':
        nu, P, A, x3 = x; T = -mp.mpf(1.05) + 2 * mp.mpf(1.05) * x3
        e = mp.mpf(0.01) * nu; ph = mp.mpf('0.25') * P * nu; qh = 1 + mp.mpf('0.25') * A * nu; tauh = mp.mpf('0.5') * T * nu
        return e, e * ph, e * qh, mp.mpf(1) / 2 + e * tauh
    if eng == 'C1pp':
        epp, l, w, tau = x; e = epp * l * l
        return e, l * mp.sin(w), l * mp.cos(w), epp ** 2 * tau / (1 - tau)
    if eng == 'C1p':
        eps, l, w, tau = x
        return eps, l * mp.sin(w), l * mp.cos(w), eps ** 2 * tau / (1 - tau)
    if eng == 'C1':
        eps, th, w, t = x; e = eps * mp.cos(th); l = mp.sqrt(eps * mp.sin(th))
        return e, l * mp.sin(w), l * mp.cos(w), t
    if eng == 'a':
        eta, l, w, t = x
        return mp.sin(eta), l * mp.sin(w), l * mp.cos(w), t
    raise ValueError(eng)


def Q_roots(e, p, q, t):
    b = mp.sqrt(1 - e * e); N = 32
    vals = []
    for n in range(N):
        u = 2 * mp.pi * n / N
        r = 1 - e * mp.cos(u); Z = p * (mp.cos(u) - e) + q * b * mp.sin(u); Zp = -p * mp.sin(u) + q * b * mp.cos(u)
        A = 2 * e * mp.sin(u) * r; rho2 = r * r - Z * Z; d = A - 2 * Z * Zp
        vals.append((1 - t) * A * A * rho2 - t * d * d)
    c = [mp.fsum(vals[n] * mp.exp(-1j * k * 2 * mp.pi * n / N) for n in range(N)) / N for k in range(-6, 7)]
    sc = max(abs(v) for v in c)
    lo = 0
    while abs(c[lo]) < mp.mpf(10) ** -70 * sc:
        lo += 1
    hi = 12
    while abs(c[hi]) < mp.mpf(10) ** -70 * sc:
        hi -= 1
    z = mp.polyroots(c[lo:hi + 1][::-1], maxsteps=500, extraprec=400)
    return list(z) + [mp.mpf(0)] * lo, sc        # exact zeros at 0 (degree drop at infinity is symmetric)


for fn in files:
    D = json.load(open(fn)); H = D['header']; eng = H['engine']; S = H['scale']
    wl = []                      # (box, W, L)
    for r in D['results']:
        stack = [tuple(r['box'])]; wi = 0
        for ch in r['tree']:
            bx = stack.pop()
            if ch == '+':
                stack.extend(split(bx, S))
            elif '1' <= ch <= '6':
                wl.append((bx, int(ch), r['wL'][wi])); wi += 1
        assert wi == len(r['wL'])
    sample = random.sample(wl, min(NLEAF, len(wl)))
    viol = 0; near = 0; n = 0; maxreal = 0; mind = mp.inf
    for bx, W, L in sample:
        for k in range(NPT):
            x = [bx[2 * i] + (bx[2 * i + 1] - bx[2 * i]) * random.random() for i in range(4)]
            if x[0] <= 1e-12:
                x[0] = 1e-12
            if eng in ('C1pp', 'C1p') and x[3] > 1 - 1e-9:
                x[3] = 1 - 1e-9
            e, p, q, t = phys(eng, x)
            z, sc = Q_roots(e, p, q, t)
            rho = mp.e ** L
            inside = sum(1 for r in z if abs(r) < rho)
            dist = min(abs(abs(r) - rho) / rho for r in z)
            on = sum(1 for r in z if r != 0 and abs(mp.log(abs(r))) < mp.mpf(10) ** -25)
            mind = min(mind, dist); maxreal = max(maxreal, on); n += 1
            if inside != W or on > 12 - 2 * W:
                viol += 1; print('   VIOLATION', eng, bx, x, 'W', W, 'inside', inside, 'on circle', on, 'L', L)
            if dist < mp.mpf(10) ** -30:
                near += 1
    print(f"{fn.split('/')[-1]}: winding leaves {len(wl)}; checked {len(sample)} leaves x {NPT} pts = {n}; violations {viol}; "
          f"roots within 1e-30 of the test circle {near}; min relative distance of a root to |z| = e^L "
          f"{mp.nstr(mind, 3)}; max real roots {maxreal}", flush=True)

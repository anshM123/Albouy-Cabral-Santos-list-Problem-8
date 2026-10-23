"""Investigate the single a/f-chart mismatch of a2 with 50-digit arithmetic."""
import os, sys
import numpy as np
import mpmath as mp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'code'))
mp.mp.dps = 60
rng = np.random.default_rng(20261022)
# regenerate the exact parameters of the mismatch (same RNG stream as a2)
for it in range(12): rng.uniform(0.01, 0.99), rng.uniform(0, np.pi), rng.uniform(0, 2*np.pi)
for it in range(12): rng.uniform(0.01, 0.99), rng.uniform(0, np.pi), rng.uniform(0, 2*np.pi)
target = None
for it in range(400):
    e = rng.uniform(0.02, 0.98); l = rng.uniform(0, 1); w = rng.uniform(0, np.pi/2)
    s = np.exp(rng.uniform(np.log(0.05), np.log(20)))
    if abs(e - 0.9436) < 5e-5 and abs(l - 0.0061) < 5e-5:
        target = (e, l, w, s)
print("params", target)
e, l, w, s = [mp.mpf(x) for x in target]
p, q = l*mp.sin(w), l*mp.cos(w); b = mp.sqrt(1-e*e); S = s*s
def pol_u(u):
    r = 1 - e*mp.cos(u); Z = p*(mp.cos(u)-e) + b*q*mp.sin(u); Zp = -p*mp.sin(u) + b*q*mp.cos(u)
    rho2 = r*r - Z*Z; A = 2*e*mp.sin(u)*r; drho2 = A - 2*Z*Zp
    return A*A*rho2 - S*drho2*drho2
def pol_f(f, kap):
    zeta = p*mp.cos(f) + q*mp.sin(f); zp = -p*mp.sin(f) + q*mp.cos(f)
    sig2 = 1 - zeta**2; D = 1 + e*mp.cos(f); B = -2*zeta*zp + 2*e*(mp.sin(f) - q*zeta)
    return 4*(1+e)**2*e**2*mp.sin(f)**2*sig2 - kap**2*D**2*B**2
# dense sign-change count + local minima of |.| (possible double roots)
def scan(fn, n=200000):
    xs = [2*mp.pi*(k+mp.mpf(1)/3)/n for k in range(n)]
    v = [fn(x) for x in xs]
    sc = sum(1 for k in range(n) if mp.sign(v[k]) != mp.sign(v[(k+1)%n]))
    return sc, xs, v
kap = s/(1-e)
for name, fn in (("u-chart", pol_u), ("f-chart", lambda f: pol_f(f, kap))):
    sc, xs, v = scan(fn, 20000)
    print(name, "sign changes on 20000-point grid:", sc)
# near u = 0 and u = pi: zoom
for c0 in (mp.mpf(0), mp.pi):
    xs = [c0 + mp.mpf(k)/2000*mp.mpf('0.02') - mp.mpf('0.01') for k in range(2001)]
    v = [pol_u(x) for x in xs]
    sc = sum(1 for k in range(2000) if mp.sign(v[k]) != mp.sign(v[k+1]))
    print("u-chart zoom around u=%.4f: sign changes %d, min %.3e max %.3e" % (float(c0), sc, float(min(v)), float(max(v))))
    f0 = 2*mp.atan(mp.sqrt((1+e)/(1-e))*mp.tan(c0/2)) if c0 == 0 else mp.pi
    xs = [f0 + mp.mpf(k)/2000*mp.mpf('0.2') - mp.mpf('0.1') for k in range(2001)]
    v = [pol_f(x, kap) for x in xs]
    sc = sum(1 for k in range(2000) if mp.sign(v[k]) != mp.sign(v[k+1]))
    print("f-chart zoom around f=%.4f: sign changes %d, min %.3e max %.3e" % (float(f0), sc, float(min(v)), float(max(v))))

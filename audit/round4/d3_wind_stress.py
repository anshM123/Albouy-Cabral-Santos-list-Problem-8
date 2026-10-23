"""AUDIT-4 d3: adversarial stress test of wind_cert.wind_box on synthetic self-inversive polynomials with KNOWN roots.

Construction: Q(z) = mu * z^m0 * prod (z - e^{i phi_j}) * prod (z - w_j)(z - 1/conj(w_j)), deg = 12 - m0 (degree drop m0
matched by m0 zeros at 0), mu = unimodular constant making c_{-k} = conj(c_k) (so F(u) = e^{-6iu} Q(e^{iu}) is real).
Families: generic; multiple roots on the unit circle; roots at 0 (m0 = 1..3); pairs with |w| very close to rho = e^L
(adversarial: tiny min|G|); W up to 6.  The box is degenerate (zero parameter radii), t0 = t1 = 0.
Check: whenever wind_box returns ok, W == #{roots with |z| < e^L} (exact by construction).  Also Q == 0 must fail.
usage: python d3_wind_stress.py N SEED   (.venv python)"""
import os, sys, math, random, cmath
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(HERE, '..', '..', 'code'))
import wind_cert
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 1)


def build(roots, m0):
    """coefficients c_{-6..6} of the self-inversive Q with given nonzero roots and m0 zeros at 0"""
    p = np.poly1d([1.0 + 0j])
    for r in roots:
        p = p * np.poly1d([1.0, -r])
    coef = np.zeros(13, complex)
    lowcoef = p.c[::-1]                       # increasing powers
    coef[m0:m0 + len(lowcoef)] = lowcoef
    # z^12 conj(Q(1/conj z)) has coefficients conj(coef[12 - j]); find lambda with that = lambda * Q
    rev = np.conj(coef[::-1])
    k = np.argmax(np.abs(coef))
    lam = rev[k] / coef[k]
    mu = cmath.sqrt(lam)
    c = mu * coef                              # now conj(c[12 - j]) == c[j]
    assert np.allclose(np.conj(c[::-1]), c, atol=1e-9 * np.abs(c).max())
    return c                                   # c[j] multiplies z^j, i.e. c_{j-6}


def box_data(c):
    a = np.zeros(7); b = np.zeros(7)
    a[0] = c[6].real
    for k in range(1, 7):
        a[k] = 2 * c[6 + k].real; b[k] = -2 * c[6 + k].imag
    cmid = np.zeros((2, 2, 7)); cmid[0, 0] = a; cmid[0, 1] = b
    crad = 1e-13 * np.abs(cmid) + 1e-300
    dcen = np.zeros((3, 2, 2, 7)); dwid = np.zeros((3, 2, 2, 7))
    cmax = np.abs(cmid) * (1 + 1e-12) + crad
    rr = np.zeros(3)
    return cmid, crad, dcen, dwid, cmax, rr


st = dict(tests=0, ok=0, wrong=0, by_family={}); ex = []
for it in range(N):
    fam = it % 5
    m0 = random.choice([0, 0, 1, 2, 3]) if fam == 3 else 0
    nslots = 12 - 2 * m0
    roots = []
    # unit-circle roots (even multiplicity patterns allowed)
    ncirc = random.choice([0, 2, 4, 6]) if fam != 1 else random.choice([2, 4, 6, 8])
    ncirc = min(ncirc, nslots)
    for j in range(ncirc):
        if fam == 1 and j % 2 == 1:
            roots.append(roots[-1])            # double roots on the circle
        else:
            roots.append(cmath.exp(1j * random.uniform(0, 2 * math.pi)))
    npair = (nslots - ncirc) // 2
    L = -10 ** random.uniform(-2, 0.7)
    rho = math.exp(L)
    for j in range(npair):
        if fam == 2 and j == 0:
            r = rho * (1 + random.choice([-1, 1]) * 10 ** random.uniform(-12, -3))   # adversarial: near the test circle
        else:
            r = math.exp(-10 ** random.uniform(-2, 0.7))
        w = r * cmath.exp(1j * random.uniform(0, 2 * math.pi))
        roots += [w, 1 / np.conj(w)]
    if len(roots) + 2 * m0 != 12:
        continue
    c = build(roots, m0)
    truth = m0 + sum(1 for r in roots if abs(r) < rho)
    ok, W = wind_cert.wind_box(*box_data(c), 0.0, 0.0, L, 64, 3, 1024)
    st['tests'] += 1
    f = st['by_family'].setdefault(fam, [0, 0, 0])
    f[0] += 1
    if ok:
        st['ok'] += 1; f[1] += 1
        if W != truth:
            st['wrong'] += 1; f[2] += 1
            ex.append((fam, m0, L, W, truth))
print("stress:", st)
for x in ex[:10]:
    print("   WRONG W", x)
# Q == 0 must fail
z = np.zeros(13, complex)
print("Q == 0:", wind_cert.wind_box(*box_data(z), 0.0, 0.0, -0.5, 64, 3, 1024))

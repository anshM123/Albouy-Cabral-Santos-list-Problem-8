"""AUDIT-8 h4: the corner e = 1, mu = 0 (t = 1) of engine KQ.  There P = -D^2 B^2 with D = 1 + cos f (double zero at
f = pi) and B(pi) = -2 Pz Qz (1 - e) = 0, so P has a zero of multiplicity >= 6 at aphelion; B (trigonometric degree 2)
has 2 or 4 real zeros.  If B has 4 real zeros, N_mult(P) = 4 + 8 = 12 at the limit point, and no box containing it can
be certified with a total <= 11 by any sound zero count.  This script finds, for each omega, the set of l in
[0.9, 1] where B has 4 real zeros (B evaluated exactly in closed form at e = 1), and checks the l = 1 factorisation
B = -2 cos(f + w) (sin(f + w) + sin w)."""
import numpy as np
import mpmath as mp


def B(f, l, w, e=1.0):
    pz, qz = l * np.sin(w), l * np.cos(w)
    zeta = pz * np.cos(f) + qz * np.sin(f); dz = -pz * np.sin(f) + qz * np.cos(f)
    return -2 * zeta * dz + 2 * e * (np.sin(f) - qz * zeta)


def nreal(l, w):
    """number of zeros of B on the circle, from its degree-4 companion polynomial (B has trigonometric degree 2)"""
    n = 16
    f = 2 * np.pi * np.arange(n) / n
    c = np.fft.fft(B(f, l, w)) / n
    cc = np.array([c[(k - 2) % n] for k in range(5)])           # coefficients of z^-2 .. z^2
    z = np.roots(cc[::-1])
    return int(np.sum(np.abs(np.log(np.abs(z))) < 1e-7))


f = np.linspace(0, 2 * np.pi, 13)
for w in (0.1, 0.3, 0.6, 1.0):
    err = np.max(np.abs(B(f, 1.0, w) + 2 * np.cos(f + w) * (np.sin(f + w) + np.sin(w))))
    print(f'l = 1 factorisation check, w = {w}: max |B + 2 cos(f+w)(sin(f+w)+sin w)| = {err:.1e}')
print('omega   smallest l in [0.9, 1] with 4 real zeros of B (e = 1)')
for w in np.linspace(0.02, 1.40, 24):
    ls = np.linspace(0.9, 1.0, 20001)
    four = [l for l in ls if nreal(l, w) >= 4]
    print(f'{w:6.3f}   {min(four):.5f}' if four else f'{w:6.3f}   none')

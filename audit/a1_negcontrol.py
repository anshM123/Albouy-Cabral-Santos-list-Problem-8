"""Negative control for a1_symbolic: the comparator must flag planted errors."""
import runpy, sys, sympy as sp
ns = runpy.run_path('a1_symbolic.py', run_name='a1')
cmp = ns['cmp']; e, b, p, q = ns['e'], ns['b'], ns['p'], ns['q']
W = ns['W']; WF = ns['WF']
import coeffs3_ce
V, D = coeffs3_ce.allc(e, b, p, q)
n = 0
for f in range(2):
    for part in range(2):
        for k in range(7):
            if W[f][part][k] != 0:
                n += 1
print("nonzero exact coefficient polys (a-chart Phi/Psi):", n)
print("planted error 1e-6*e*q on V[1][0][3]:", cmp('x', V[1][0][3] + 1e-6 * e * q, W[1][0][3], (e, b, p, q)))
print("swapped Phi<->Psi a_2:", cmp('x', V[0][0][2], W[1][0][2], (e, b, p, q)))
print("sample exact poly Psi a_3:", W[1][0][3])

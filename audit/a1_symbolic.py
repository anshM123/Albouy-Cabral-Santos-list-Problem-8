"""AUDIT a1: independent symbolic re-derivation of every generated coefficient file.

Independent pipeline: real trigonometric-polynomial algebra (product-to-sum formulas on dicts
{k: coeff} for cos(ku) and sin(ku)), no z = e^{iu} substitution, no sympy Poly-in-z, no re/im.
Then every generated function is evaluated on sympy SYMBOLS and compared exactly (up to float
representation of rational constants, tolerance 1e-14 on each monomial coefficient):
  coeffs_ce.coeffs        (complex coefficients of Q(z) = z^6 Pol, with S)
  coeffs2_ce.phi_ab/psi_ab
  coeffs3_ce.allc         (values + partial derivatives in e, b, p, q)
  coeffsF.allcF           (f-chart values + partials in e, p, q)
  nbcoef_a / nbcoef_f     (numba code run as plain python on an object array: values, partials,
                           and the absolute-value polynomials, compared with sum |c_a| y^a)
Run with system python (sympy).  Usage: python a1_symbolic.py
"""
import os
import re
import sys
import types
import numpy as np
import sympy as sp

CODE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'code')
sys.path.insert(0, CODE)

e, b, p, q, S = sp.symbols('e b p q S', real=True)


class TP:
    """real trig polynomial  c[0] + sum_k c[k] cos(ku) + s[k] sin(ku)"""

    def __init__(self, c=None, s=None):
        self.c = {k: v for k, v in (c or {}).items()}
        self.s = {k: v for k, v in (s or {}).items()}

    @staticmethod
    def const(x):
        return TP({0: sp.sympify(x)})

    def _norm(self):
        self.c = {k: sp.expand(v) for k, v in self.c.items() if sp.expand(v) != 0}
        self.s = {k: sp.expand(v) for k, v in self.s.items() if k != 0 and sp.expand(v) != 0}
        return self

    def __add__(self, o):
        o = o if isinstance(o, TP) else TP.const(o)
        c = dict(self.c); s = dict(self.s)
        for k, v in o.c.items():
            c[k] = c.get(k, 0) + v
        for k, v in o.s.items():
            s[k] = s.get(k, 0) + v
        return TP(c, s)._norm()

    __radd__ = __add__

    def __neg__(self):
        return TP({k: -v for k, v in self.c.items()}, {k: -v for k, v in self.s.items()})

    def __sub__(self, o):
        o = o if isinstance(o, TP) else TP.const(o)
        return self + (-o)

    def __rsub__(self, o):
        return TP.const(o) - self

    def __mul__(self, o):
        if not isinstance(o, TP):
            o = sp.sympify(o)
            return TP({k: v * o for k, v in self.c.items()}, {k: v * o for k, v in self.s.items()})._norm()
        c = {}; s = {}

        def addc(k, v):
            k2 = abs(k)
            c[k2] = c.get(k2, 0) + v

        def adds(k, v):        # v sin(k u), sin(-ku) = -sin(ku)
            if k == 0:
                return
            if k < 0:
                k, v = -k, -v
            s[k] = s.get(k, 0) + v

        for m, x in self.c.items():
            for n, y in o.c.items():          # cos m cos n
                addc(m - n, x * y / 2); addc(m + n, x * y / 2)
            for n, y in o.s.items():          # cos m sin n = 1/2 [sin(m+n) - sin(m-n)]
                adds(m + n, x * y / 2); adds(m - n, -x * y / 2)
        for m, x in self.s.items():
            for n, y in o.c.items():          # sin m cos n = 1/2 [sin(m+n) + sin(m-n)]
                adds(m + n, x * y / 2); adds(m - n, x * y / 2)
            for n, y in o.s.items():          # sin m sin n = 1/2 [cos(m-n) - cos(m+n)]
                addc(m - n, x * y / 2); addc(m + n, -x * y / 2)
        return TP(c, s)._norm()

    __rmul__ = __mul__

    def d(self):
        """derivative in u"""
        c = {}; s = {}
        for k, v in self.c.items():
            if k:
                s[k] = s.get(k, 0) - k * v
        for k, v in self.s.items():
            c[k] = c.get(k, 0) + k * v
        return TP(c, s)._norm()

    def ab(self):
        a = [sp.expand(self.c.get(k, 0)) for k in range(7)]
        s = [sp.expand(self.s.get(k, 0)) if k else sp.Integer(0) for k in range(7)]
        assert max(list(self.c) + list(self.s) + [0]) <= 6
        return a, s


cu = TP({1: 1}); su = TP(s={1: 1})

# ---------------- a-chart (eccentric anomaly) ----------------
r = 1 - e * cu
Z = p * (cu - e) + b * q * su
rho2 = r * r - Z * Z
A = (r * r).d()
Phi = A * A * rho2
Psi = rho2.d() * rho2.d()
# independent check: A = 2 e sin u r
assert ((A - 2 * e * su * r).c == {} and (A - 2 * e * su * r).s == {})
PhiA, PhiS = Phi.ab()
PsiA, PsiS = Psi.ab()

# ---------------- f-chart (true anomaly) ----------------
cf, sf = cu, su
zeta = p * cf + q * sf
sig2 = 1 - zeta * zeta
Dd = 1 + e * cf
B = sig2.d() + 2 * e * (sf - q * zeta)
Bfull = sig2.d() * Dd + 2 * e * sf * sig2
dB = B - Bfull
assert dB.c == {} and dB.s == {}, "B simplification identity fails"
PhiF = 4 * (1 + e) ** 2 * e ** 2 * (sf * sf) * sig2
PsiF = Dd * Dd * B * B
PhiFA, PhiFS = PhiF.ab()
PsiFA, PsiFS = PsiF.ab()

TOL = 1e-13
report = []


def cmp(name, got, want, vars_):
    """compare polynomial `got` (may contain Floats) with exact `want`; return max abs coefficient diff"""
    diff = sp.expand(sp.sympify(got) - want)
    if diff == 0:
        return 0.0
    P = sp.Poly(diff, *vars_)
    m = max(abs(float(c)) for c in P.coeffs())
    return m


worst = {}


def note(tag, val):
    worst[tag] = max(worst.get(tag, 0.0), val)


# ---- coeffs2_ce ----
import coeffs2_ce
for fname, (WA, WS) in (("phi_ab", (PhiA, PhiS)), ("psi_ab", (PsiA, PsiS))):
    a, s = getattr(coeffs2_ce, fname)(e, b, p, q)
    for k in range(7):
        note('coeffs2 ' + fname, cmp(fname, a[k], WA[k], (e, b, p, q)))
        note('coeffs2 ' + fname, cmp(fname, s[k], WS[k], (e, b, p, q)))

# ---- coeffs3_ce ----
import coeffs3_ce
V, D = coeffs3_ce.allc(e, b, p, q)
W = [[PhiA, PhiS], [PsiA, PsiS]]
for f in range(2):
    for part in range(2):
        for k in range(7):
            note('coeffs3 V', cmp('V', V[f][part][k], W[f][part][k], (e, b, p, q)))
            for vi, var in enumerate((e, b, p, q)):
                note('coeffs3 D', cmp('D', D[vi][f][part][k], sp.diff(W[f][part][k], var), (e, b, p, q)))

# ---- coeffsF ----
import coeffsF
VF, DF = coeffsF.allcF(e, p, q)
WF = [[PhiFA, PhiFS], [PsiFA, PsiFS]]
for f in range(2):
    for part in range(2):
        for k in range(7):
            note('coeffsF V', cmp('VF', VF[f][part][k], WF[f][part][k], (e, p, q)))
            for vi, var in enumerate((e, p, q)):
                note('coeffsF D', cmp('DF', DF[vi][f][part][k], sp.diff(WF[f][part][k], var), (e, p, q)))

# ---- coeffs_ce: Q(z) = sum_k c_k z^k = z^6 Pol, Pol = Phi - S Psi ----
import coeffs_ce
re_, im_ = coeffs_ce.coeffs(e, b, p, q, S)
PolA = [sp.expand(PhiA[k] - S * PsiA[k]) for k in range(7)]
PolS = [sp.expand(PhiS[k] - S * PsiS[k]) for k in range(7)]
for kk in range(13):
    m = kk - 6
    if m == 0:
        wr, wi = PolA[0], sp.Integer(0)
    elif m > 0:          # a cos + s sin = Re[(a - i s) e^{imu}] * ... : c_m = (a_m - i s_m)/2
        wr, wi = PolA[m] / 2, -PolS[m] / 2
    else:
        wr, wi = PolA[-m] / 2, PolS[-m] / 2
    note('coeffs_ce re', cmp('re', re_[kk], wr, (e, b, p, q, S)))
    note('coeffs_ce im', cmp('im', im_[kk], wi, (e, b, p, q, S)))


# ---- nbcoef_a / nbcoef_f run as plain python on object arrays ----
def load_plain(path):
    src = open(path).read()
    src = src.replace('from numba import njit', '')
    src = re.sub(r'@njit\([^)]*\)\n', '', src)
    mod = types.ModuleType(os.path.basename(path))
    exec(compile(src, path, 'exec'), mod.__dict__)
    return mod


def abs_poly(expr, vars_, ys):
    if expr == 0:
        return sp.Integer(0)
    P = sp.Poly(sp.expand(expr), *vars_)
    return sp.expand(sum(abs(c) * sp.Mul(*[y ** d for y, d in zip(ys, mon)]) for mon, c in P.terms()))


for chart, vars_, WW in (('a', (e, b, p, q), W), ('f', (e, p, q), WF)):
    mod = load_plain(os.path.join(CODE, f'nbcoef_{chart}.py'))
    nv = len(vars_)
    out = np.empty((nv + 1, 2, 2, 7), dtype=object)
    mod.coef_val(*vars_, out)
    ys = sp.symbols('y0:%d' % nv, positive=True)
    outa = np.empty((nv + 1, 2, 2, 7), dtype=object)
    mod.coef_abs(*ys, outa)
    for f in range(2):
        for part in range(2):
            for k in range(7):
                want = WW[f][part][k]
                note(f'nbcoef_{chart} value', cmp('v', out[0, f, part, k], want, vars_))
                note(f'nbcoef_{chart} |value|', cmp('va', outa[0, f, part, k], abs_poly(want, vars_, ys), ys))
                for vi, var in enumerate(vars_):
                    dw = sp.diff(want, var)
                    note(f'nbcoef_{chart} partial', cmp('d', out[vi + 1, f, part, k], dw, vars_))
                    note(f'nbcoef_{chart} |partial|', cmp('da', outa[vi + 1, f, part, k], abs_poly(dw, vars_, ys), ys))

print("max |coefficient difference| per generated object (0 = identical):")
for k, v in worst.items():
    print(f"  {k:28s} {v:.3e}   {'OK' if v < TOL else 'MISMATCH'}")

# ---- extra structural facts used in STATEMENT / LOG ----
# planar case: Pol = A^2 (r^2 - S) at p = q = 0
pl = (Phi - S * Psi)
plA, plS = pl.ab()
target = (A * A * (r * r - S))
tA, tS = target.ab()
ok_planar = all(sp.expand((plA[k] - tA[k]).subs({p: 0, q: 0})) == 0 and sp.expand((plS[k] - tS[k]).subs({p: 0, q: 0})) == 0
                for k in range(7))
print("planar identity Pol = A^2 (r^2 - S) at p = q = 0:", ok_planar)
# e = 0, p = q = 0: Pol identically 0 (concentric coplanar circles)
print("Pol == 0 at e = 0, p = q = 0:", all(sp.expand(x.subs({e: 0, p: 0, q: 0})) == 0 for x in PolA + PolS))

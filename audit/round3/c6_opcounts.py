"""AUDIT-3 c6: operation counts in every generated coefficient module (code/nbcoef_*.py, code/c3/nbcoef_c3_*.py):
max terms per polynomial T, max total degree of a monomial (= multiplications per term incl. the constant), max power
chain; gamma_n with n = deg + T - 1 bounds the evaluation error relative to Pabs (Higham, sums of products)."""
import glob, os, re
u = 2.0 ** -53
rows = []
for path in sorted(glob.glob('../../code/nbcoef_*.py') + glob.glob('../../code/c3/nbcoef_c3_*.py')):
    src = open(path).read()
    body = src.split('def coef_val')[1].split('def coef_abs')[0]
    T = 0; D = 0; P = 0; worst_n = 0
    for rhs in re.findall(r'out\[[^\]]+\] = (.*)', body):
        terms = rhs.split(' + ')
        if terms == ['0.0']:
            continue
        T = max(T, len(terms))
        dmax = 0
        for t in terms:
            fac = t.split(' * ')[1:]
            deg = 0
            for f in fac:
                m = re.match(r'([A-Za-z_]+?)(\d+)$', f)
                if m and f != '1.0':
                    deg += int(m.group(2)); P = max(P, int(m.group(2)))
            dmax = max(dmax, deg)
        D = max(D, dmax)
        worst_n = max(worst_n, dmax + len(terms) - 1)
    rows.append((os.path.basename(path), T, D, P, worst_n))
for r in rows:
    print("%-26s terms %4d  degree %3d  max power %3d  n = %4d" % r)
n = max(r[4] for r in rows)
print("max n =", n, "; gamma_n =", n * u / (1 - n * u), "; + constant rounding u =", u, "; total < 1e-13:", n * u / (1 - n * u) + u < 1e-13)

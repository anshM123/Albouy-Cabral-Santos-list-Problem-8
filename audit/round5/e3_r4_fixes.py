"""AUDIT-5 e3: R4 fixes.  (i) crafted near-antipodal nodes (AUDIT_4 d4) through the new wind_cert.arg_in_branch loop;
(ii) guards; (iii) W reproduction with the FIXED wind_box on every winding leaf of the given runs."""
import os, sys, json, math, importlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.join(HERE, '..', '..', 'code'); sys.path.insert(0, CODE)
import wind_cert, nb_cert
eps = 1e-30
nodes = [(-eps, -1.0), (-eps, 1.0), (1.0, eps), (eps, -1.0)]; labs = [3, 2, 1, 4]
tot = 0.0
for i in range(4):
    j = (i + 1) % 4
    a0 = wind_cert.arg_in_branch(labs[i], *nodes[i]); a1 = wind_cert.arg_in_branch(labs[i], *nodes[j])
    tot += a1 - a0
print("(i) crafted nodes, new branch logic: total = %.15f -> winding %d (true -1)" % (tot, round(tot / (2 * math.pi))))
z = np.zeros((2, 2, 7)); z[0, 0, 0] = 1.0
cm = (z, 1e-13 * np.abs(z) + 1e-300, np.zeros((3, 2, 2, 7)), np.zeros((3, 2, 2, 7)), np.abs(z) * (1 + 1e-12), np.zeros(3))
print("(ii) guards: L = 0.1 ->", wind_cert.wind_box(*cm, 0.0, 0.0, 0.1, 64, 3, 1024), "; L = -20 ->", wind_cert.wind_box(*cm, 0.0, 0.0, -20.0, 64, 3, 1024),
      "; F = 1 (W = 6 + 0 = 6? no roots) at L = -1 ->", wind_cert.wind_box(*cm, 0.0, 0.0, -1.0, 64, 3, 1024))
def split(bx, S):
    w = [(bx[2 * i + 1] - bx[2 * i]) / S[i] for i in range(4)]
    i = int(np.argmax(w)); mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)
for fn in sys.argv[1:]:
    D = json.load(open(fn)); H = D['header']; eng = H['engine']; S = H['scale']
    if eng in ('a', 'f'):
        ch = 0 if eng == 'a' else 1; bc = lambda bx: nb_cert.box_coeffs(ch, *bx[:6])
    elif eng in ('C1', 'C1p', 'C1pp'):
        m = importlib.import_module('nbeng_' + eng); bc = lambda bx, m=m: m.box_coeffs(*bx[:6])
    else:
        m = importlib.import_module('nbeng_' + eng); bc = lambda bx, m=m: m.box_coeffs(*bx)
    n = 0; bad = 0
    for r in D['results']:
        stack = [tuple(r['box'])]; wi = 0
        for c in r['tree']:
            bx = stack.pop()
            if c == '+': stack.extend(split(bx, S))
            elif '1' <= c <= '6':
                L = r['wL'][wi]; wi += 1; n += 1
                ok, W = wind_cert.wind_box(*bc(bx), bx[6], bx[7], L, H.get('wind_m0', 64), H.get('wind_rounds', 3), H.get('wind_cap', 1024))
                if not ok or W != int(c): bad += 1
    print(f"(iii) {os.path.basename(fn)}: winding leaves re-certified with the fixed wind_box: {n}; mismatches {bad}", flush=True)

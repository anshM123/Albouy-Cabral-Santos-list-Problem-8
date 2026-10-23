"""Re-certify every winding leaf ('1'..'6') of a bb_run JSON with the current wind_cert (parallel).
Used for runs that were started with the first version of wind_cert (bulk, C1p), whose winding sum was later made
branch-exact (audit R4-1).  A leaf passes if wind_box(box, stored L) certifies the same W (hence the same bound).
usage: python recert_wind.py RUN.json WORKERS"""
import json
from runio import load_json
import math
import sys
import time
import numpy as np
from multiprocessing import Pool

ENG = None


def init(engine):
    global bcf, wb
    import wind_cert
    wb = wind_cert.wind_box
    if engine in ('a', 'f'):
        import nb_cert
        ch = 0 if engine == 'a' else 1
        bcf = lambda bx: nb_cert.box_coeffs(ch, *bx[:6])
    else:
        import importlib
        m = importlib.import_module('nbeng_' + engine)
        bcf = (lambda bx: m.box_coeffs(*bx[:6])) if engine in ('C1', 'C1p', 'C1pp') else (lambda bx: m.box_coeffs(*bx))


def check(args):
    items, m0, rounds, cap = args
    bad = []
    for bx, W, L in items:
        ok, W2 = wb(*bcf(bx), bx[6], bx[7], L, m0, rounds, cap)
        if not ok or W2 != W:
            bad.append((bx, W, L, bool(ok), int(W2)))
    return len(items), bad


def winding_leaves(D):
    H = D['header']; scale = H['scale']
    out = []
    for r in D['results']:
        stack = [tuple(r['box'])]; wi = 0; wl = r.get('wL', [])
        for ch in r['tree']:
            bx = stack.pop()
            if ch == '+':
                w = [(bx[2 * i + 1] - bx[2 * i]) / scale[i] for i in range(4)]
                i = int(np.argmax(w)); mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
                a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
                stack.extend([tuple(a), tuple(c)])
            elif '1' <= ch <= '6':
                out.append((bx, ord(ch) - ord('0'), wl[wi])); wi += 1
    return out


if __name__ == '__main__':
    fn = sys.argv[1]; workers = int(sys.argv[2])
    D = load_json(fn)
    H = D['header']
    leaves = winding_leaves(D)
    chunks = [leaves[i:i + 2000] for i in range(0, len(leaves), 2000)]
    t0 = time.time(); n = 0; bad = []
    with Pool(workers, initializer=init, initargs=(H['engine'],)) as pool:
        for cnt, b in pool.imap_unordered(check, [(c, H.get('wind_m0', 64), H.get('wind_rounds', 3), H.get('wind_cap', 1024))
                                                  for c in chunks]):
            n += cnt; bad += b
    print(f'{fn}: re-certified {n} winding leaves with the current wind_cert; problems {len(bad)}; {time.time() - t0:.0f}s')
    json.dump(dict(file=fn, checked=n, problems=bad), open(fn + '.recert_wind.json', 'w'))

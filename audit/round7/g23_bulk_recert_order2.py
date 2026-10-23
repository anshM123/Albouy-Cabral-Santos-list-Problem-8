"""AUDIT-7 g23: do leaves of the bulk run made with the FIRST-ORDER engine (the records completed before the 09:26
resume with --order2) re-certify with the second-order engine that bb_replay now uses (nb_cert2.certify_box_deep2 with
total <= 11; nb_cert2.wind_box2 on box_data2(0, ...) with the stored L and wind rounds 8 / cap 16384)?  And do leaves made
after the resume re-certify with the first-order engine (for information)?  Sample from the records file.
usage: python g23_bulk_recert_order2.py N SEED"""
import os, sys, json, math, random
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.normpath(os.path.join(HERE, '..', '..', 'code')); sys.path.insert(0, CODE)
import nb_cert, nb_cert2, wind_cert
N, SEED = int(sys.argv[1]), int(sys.argv[2]); random.seed(SEED)
S = (math.pi / 2, 1.0, math.pi / 2, 1.0)
recs = []
for line in open(os.path.join(CODE, 'runs', 'bulk.json.records.jsonl')):
    try:
        r = json.loads(line)
    except Exception:
        continue
    if r.get('done'):
        recs.append(r)
first, second = recs[:3825], recs[3825:]          # 3825 records were completed with the first-order engine


def leaves(rs, n):
    reg, wind, seen = [], [], [0, 0]
    for r in rs:
        stack = [tuple(r['box'])]; wl = r.get('wL', []); wi = 0
        for ch in r['tree']:
            bx = stack.pop()
            if ch == '+':
                w = [(bx[2 * i + 1] - bx[2 * i]) / S[i] for i in range(4)]
                i = max(range(4), key=lambda j: w[j]); mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
                a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
                stack += [tuple(a), tuple(c)]
            elif '1' <= ch <= '6':
                seen[1] += 1
                item = (bx, int(ch), wl[wi]); wi += 1
                if len(wind) < n: wind.append(item)
                else:
                    j = random.randrange(seen[1])
                    if j < n: wind[j] = item
            elif 'a' <= ch <= 'l' or 'A' <= ch <= 'L':
                seen[0] += 1
                item = (bx, ord(ch) - (ord('a') if ch >= 'a' else ord('A')))
                if len(reg) < n: reg.append(item)
                else:
                    j = random.randrange(seen[0])
                    if j < n: reg[j] = item
    return reg, wind, seen


for name, rs, order2 in (('first-order records', first, True), ('second-order records', second, False)):
    reg, wind, seen = leaves(rs, N)
    bad_r = bad_w = 0
    for bx, b in reg:
        if order2:
            ok, tot = nb_cert2.certify_box_deep2(0, *bx, 6, 64, 16, 8192)
        else:
            ok, tot = nb_cert.certify_box_deep(0, *bx, 6, 64, 16, 8192)
        bad_r += not (ok and tot <= 11)
    for bx, W, L in wind:
        if order2:
            ok, W2 = nb_cert2.wind_box2(*nb_cert2.box_data2(0, *bx[:6]), bx[6], bx[7], L, 64, 8, 16384)
        else:
            ok, W2 = wind_cert.wind_box(*nb_cert.box_coeffs(0, *bx[:6]), bx[6], bx[7], L, 64, 8, 16384)
        bad_w += not (ok and W2 == W)
    eng = 'second-order engine (bb_replay)' if order2 else 'first-order engine'
    print(f'{name} ({len(rs)}): {len(reg)} of {seen[0]} regular/deep leaves and {len(wind)} of {seen[1]} winding leaves '
          f're-certified with the {eng}: failures {bad_r} / {bad_w}', flush=True)

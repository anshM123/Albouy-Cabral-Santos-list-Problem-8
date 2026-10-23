"""AUDIT-6 f5: re-certify sampled KQ leaves with kq_run.certify (must return ok with total <= recorded bound)."""
import os, sys, json, math, random
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.join(HERE, '..', '..', 'code')
sys.path.insert(0, CODE); sys.path.insert(0, os.path.join(CODE, 'c3'))
os.chdir(CODE)
import kq_run
kq_run.init()
random.seed(int(sys.argv[2])); N = int(sys.argv[1])
S = (1.0, 1.0, math.pi / 2, 1.0)
def split(bx):
    w = [(bx[2 * i + 1] - bx[2 * i]) / S[i] for i in range(4)]
    i = int(np.argmax(w)); mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)
for fn in sys.argv[3:]:
    leaves = []
    for line in open(os.path.join(HERE, fn)):
        r = json.loads(line)
        if not r['done']: continue
        stack = [tuple(r['box'])]
        for ch in r['tree']:
            bx = stack.pop()
            if ch == '+': stack.extend(split(bx))
            elif 'a' <= ch <= 'l': leaves.append((bx, ord(ch) - ord('a')))
    prob = 0
    sample = random.sample(leaves, min(N, len(leaves)))
    for bx, b in sample:
        ok, tot = kq_run.certify(bx)
        if not ok or tot > b: prob += 1; print('  problem', bx, b, ok, tot)
    print(f"{fn}: leaves {len(leaves)}, re-certified {len(sample)}, problems {prob}", flush=True)

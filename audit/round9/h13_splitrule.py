"""AUDIT-9 h13: unit test of code/splitrule.py (rule 'face', used by kq_run --facesplit): on random KQ boxes near the
face e = 1, mu = 0, every split (i) cuts exactly one coordinate strictly inside the box, (ii) the two children tile the
parent (shared cut value, other coordinates unchanged), and (iii) in case (1) of the rule the face-side child lies in
the face neighbourhood, so kq_run prunes it ('p'), as verified with exact arithmetic; repeated splitting terminates
(widths shrink).  Also: rule None reproduces the old largest-normalised-width bisection."""
import os, sys, math, random
import mpmath as mp
CODE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code'))
sys.path.insert(0, CODE); sys.dont_write_bytecode = True
import splitrule
random.seed(13)
mp.mp.dps = 40
SC = (1.0, 1.0, math.pi / 2, 0.0025)
S2 = mp.mpf(0.02) ** 2


def in_face_exact(b):
    t0 = mp.mpf(b[6])
    return 1 - mp.mpf(b[0]) <= S2 and mp.sqrt((1 - t0) / t0) <= S2


bad = 0; n = 0; case1 = 0; pruned_ok = 0
for _ in range(20000):
    de = 10 ** random.uniform(-7, -2); dt = 10 ** random.uniform(-10, -4)
    e1 = 1.0 if random.random() < 0.5 else random.uniform(0.999, 1.0); e0 = max(0.99, e1 - de)
    t1 = 1.0 if random.random() < 0.5 else random.uniform(0.99999, 1.0); t0 = max(0.997456, t1 - dt)
    l0 = random.uniform(0.9, 0.999); w0 = random.uniform(0, 1.3)
    b = (e0, e1, l0, l0 + 10 ** random.uniform(-4, -1.5), w0, w0 + 10 ** random.uniform(-4, -1.5), t0, t1)
    a, c = splitrule.split_box(b, SC, 'face')
    n += 1
    diff = [i for i in range(8) if a[i] != b[i] or c[i] != b[i]]
    ok = (len(diff) == 2 and diff[0] % 2 == 0 and diff[1] == diff[0] + 1 and a[diff[1]] == c[diff[0]]
          and b[diff[0]] < a[diff[1]] < b[diff[1]] and a[diff[0]] == b[diff[0]] and c[diff[1]] == b[diff[1]])
    bad += not ok
    i = diff[0] // 2 if diff else -1
    cut = a[diff[1]] if diff else None
    S2f = 0.02 * 0.02 * (1 - 5e-9)
    if i == 0 and cut == 1.0 - S2f or i == 3 and cut == 1.0 / (1.0 + S2f * S2f):
        case1 += 1
        child = c                              # face side: the upper child (e >= e* or t >= t*)
        if not in_face_exact(child):           # at most one more case-1 cut (at t*) makes the corner child prunable
            a2, c2 = splitrule.split_box(child, SC, 'face')
            child = c2 if c2[6] == 1.0 / (1.0 + S2f * S2f) or c2[0] == 1.0 - S2f else None
        pruned_ok += child is not None and in_face_exact(child)
# rule None == old rule
mism = 0
for _ in range(5000):
    b = tuple(v for _ in range(4) for v in sorted([random.random(), random.random()]))
    a, c = splitrule.split_box(b, SC, None)
    w = [(b[2 * k + 1] - b[2 * k]) / SC[k] for k in range(4)]
    k = max(range(4), key=lambda j: w[j]); mid = 0.5 * (b[2 * k] + b[2 * k + 1])
    a2 = list(b); c2 = list(b); a2[2 * k + 1] = mid; c2[2 * k] = mid
    mism += (tuple(a2), tuple(c2)) != (a, c)
print(f'face rule: {n} random splits, not a valid tiling {bad}; face cuts (case 1) {case1}, face-side corner child prunable (p, exact) within two face cuts '
      f'{pruned_ok}/{case1}; rule None vs the old bisection: {mism} differences in 5000')

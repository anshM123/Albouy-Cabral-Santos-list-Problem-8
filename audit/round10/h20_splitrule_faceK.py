"""AUDIT-10b h20: splitrule 'face/K' (K = 30) and the unchanged rules 'face' / None against splitrule_v1_backup.py:
tiling of every split (exactly one coordinate cut strictly inside, children sharing the cut, other coordinates kept),
and bit-identical results for 'face' and None."""
import os, sys, math, random, importlib.util
CODE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code'))
sys.path.insert(0, CODE); sys.dont_write_bytecode = True
import splitrule
spec = importlib.util.spec_from_file_location('splitrule_v1', os.path.join(CODE, 'splitrule_v1_backup.py'))
old = importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
random.seed(20)
SC = (1.0, 1.0, math.pi / 2, 0.0025)


def tiles(b, a, c):
    d = [i for i in range(8) if a[i] != b[i] or c[i] != b[i]]
    return (len(d) == 2 and d[0] % 2 == 0 and d[1] == d[0] + 1 and a[d[1]] == c[d[0]] and b[d[0]] < a[d[1]] < b[d[1]]
            and a[d[0]] == b[d[0]] and c[d[1]] == b[d[1]])


bad = {'face/30': 0, 'face': 0, None: 0}; diff = {'face': 0, None: 0}; n = 0; finer = 0
for _ in range(20000):
    de = 10 ** random.uniform(-7, -2); dt = 10 ** random.uniform(-10, -4)
    e1 = 1.0 if random.random() < 0.5 else random.uniform(0.999, 1.0); e0 = max(0.99, e1 - de)
    t1 = 1.0 if random.random() < 0.5 else random.uniform(0.99999, 1.0); t0 = max(0.997456, t1 - dt)
    l0 = random.uniform(0.9, 0.999); w0 = random.uniform(0, 1.3)
    b = (e0, e1, l0, l0 + 10 ** random.uniform(-4, -1.5), w0, w0 + 10 ** random.uniform(-4, -1.5), t0, t1)
    n += 1
    for rule in ('face/30', 'face', None):
        a, c = splitrule.split_box(b, SC, rule)
        bad[rule] += not tiles(b, a, c)
        if rule != 'face/30':
            diff[rule] += (a, c) != old.split_box(b, SC, rule)
    a30, _ = splitrule.split_box(b, SC, 'face/30'); a1, _ = splitrule.split_box(b, SC, 'face')
    i30 = [i for i in range(8) if a30[i] != b[i]][0] // 2; i1 = [i for i in range(8) if a1[i] != b[i]][0] // 2
    finer += (i30 in (1, 2)) and (i1 in (0, 3))
print(f'{n} random boxes near the face: invalid tilings face/30 {bad["face/30"]}, face {bad["face"]}, None {bad[None]}; '
      f'differences to splitrule_v1: face {diff["face"]}, None {diff[None]}; boxes where face/30 cuts l or w but face cuts '
      f'e or t: {finer}')

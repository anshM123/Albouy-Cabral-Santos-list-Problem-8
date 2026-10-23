"""AUDIT-11 i1 (light, exact where it matters): coverage logic of the corner-K repair.
 1. K = {Pz >= 0.875, 0.125 <= Qz <= 0.25} inside K_lw = {l >= 0.88, w >= 1.29} (exact).
 2. KQF2 domain contains {eps, mu <= SF2^2} x K_lw, with the 'p' tolerance of KQS (exact, float values of the files).
 3. KQS domain contains {eps, mu <= SF^2 (1 - 1e-9)} x K_lw (the 'p' leaves of the other KQ runs) (exact).
 4. The seam tau = T_X2: header T_X of KQF2 versus 3 SF2 and the face-chart sigma range (exact).
 5. Owed boxes of the face base files and of KQF: which are excused by which rule; what is still owed.
 6. shell/30 versus face/30 (code/splitrule.py): equal splits off the cut, tiling."""
import sys, json, math, random
import mpmath as mp
sys.path.insert(0, r'<repo>\code')
import splitrule
mp.mp.dps = 50
CODE = r'<repo>\code'
X = lambda v: mp.mpf(v)                                   # exact value of a float
SF, SF2, HP = 0.02, 0.010005, 1.5707963267948966
W9 = 9 * math.pi / 20
print('1. min l on K =', mp.nstr(mp.sqrt(X(0.875) ** 2 + X(0.125) ** 2), 12), '>= fl(0.88) =', mp.nstr(X(0.88), 17),
      mp.sqrt(X(0.875) ** 2 + X(0.125) ** 2) >= X(0.88))
print('   min w on K = atan2(0.875, 0.25) =', mp.nstr(mp.atan2(X(0.875), X(0.25)), 12), '>= fl(1.29) =', mp.nstr(X(1.29), 17),
      mp.atan2(X(0.875), X(0.25)) >= X(1.29), '; max w on K < pi/2 (Qz >= 0.125 > 0)')


def hdr(path):
    with open(path, 'rb') as fh:
        head = fh.read(4000).decode('latin1')
    i = head.index('"header"'); j = head.index('{', i); depth = 0
    for k in range(j, len(head)):
        depth += head[k] == '{'; depth -= head[k] == '}'
        if depth == 0:
            return json.loads(head[j:k + 1])


Hf2 = hdr(CODE + r'\runs\KQF2.json'); Hs = hdr(CODE + r'\runs\KQS.json'); Hf = hdr(CODE + r'\runs\KQF.json')
d2, ds = Hf2['domain'], Hs['domain']
fs2 = SF2 * SF2 * (1 - 1e-9)                              # the float threshold of the 'p' test (kq_run.in_face, check_all)
fsF = SF * SF * (1 - 1e-9)
print('2. KQF2: 1 - e_lo =', mp.nstr(1 - X(d2[0]), 17), ' SF2^2 =', mp.nstr(X(SF2) ** 2, 17), ' margin', mp.nstr(1 - X(d2[0]) - X(SF2) ** 2, 4),
      '; mu_max =', mp.nstr(mp.sqrt((1 - X(d2[6])) / X(d2[6])), 17), ' margin', mp.nstr(mp.sqrt((1 - X(d2[6])) / X(d2[6])) - X(SF2) ** 2, 4))
print('   KQF2 l, w ranges', d2[2:6], 'contain K_lw:', d2[2] <= 0.88 and d2[3] >= 1.0 and d2[4] <= 1.29 and d2[5] >= HP)
print('   KQS p-test threshold (float)', repr(fs2), '<= SF2^2 exactly:', X(fs2) <= X(SF2) ** 2)
print('3. KQS: 1 - e_lo =', mp.nstr(1 - X(ds[0]), 17), ' vs SF^2 =', mp.nstr(X(SF) ** 2, 17), ' (', 1 - X(ds[0]) >= X(SF) ** 2,
      '), vs the other runs\' p threshold', repr(fsF), '(', 1 - X(ds[0]) >= X(fsF), ')')
print('   KQS mu_max =', mp.nstr(mp.sqrt((1 - X(ds[6])) / X(ds[6])), 12), '>= SF^2:', mp.sqrt((1 - X(ds[6])) / X(ds[6])) >= X(SF) ** 2,
      '; l, w ranges', ds[2:6])
print('   KQF (original) domain', Hf['domain'], ': 1 - e_lo =', mp.nstr(1 - X(Hf['domain'][0]), 17), '>= p threshold:', 1 - X(Hf['domain'][0]) >= X(fsF))
print('4. KQF2 header T_X', repr(Hf2['T_X']), ' 3 fl(SF2) exact =', mp.nstr(3 * X(SF2), 20), ' T_X - 3 SF2 =', mp.nstr(X(Hf2['T_X']) - 3 * X(SF2), 3))
for c in ('FE', 'FM', 'FIP', 'FIM'):
    h = hdr(CODE + r'\c3\runs\%s.json' % c)
    print('   %s header: chart %s, domain %s, n0 %s, SF %s, YF %s, R1 %s' % (c, h.get('chart'), h.get('domain'), h.get('n0'), h.get('SF'),
                                                                        h.get('YF', h.get('Y_F')), h.get('R1')))
# 5. owed boxes
K = lambda b: b[4] >= 0.875 and b[6] >= 0.125 and b[7] <= 0.25
for c in ('FE', 'FM', 'FIP', 'FIM'):
    D = json.load(open(CODE + r'\c3\runs\%s.json' % c))
    ow = [tuple(b) for r in D['results'] for b in r['left']] + [tuple(f[0]) for r in D['results'] for f in r['fails']]
    a = [b for b in ow if b[0] > SF]; rest = [b for b in ow if not b[0] > SF]
    k = [b for b in rest if K(b) and b[0] > SF2]; left = [b for b in rest if not (K(b) and b[0] > SF2)]
    inK_lo = [b for b in rest if K(b) and not b[0] > SF2]
    nrec = len(D['results']); inc = sum(1 for r in D['results'] if not r['done'])
    sl = min((b[0] for b in ow), default=None)
    print(f'5. {c}: records {nrec} (incomplete {inc}), owed {len(ow)}: s_lo > SF {len(a)}, inside K with s_lo > SF2 {len(k)}, '
          f'inside K but s_lo <= SF2 {len(inK_lo)}, outside K {len([b for b in rest if not K(b)])} -> still owed {len(left)}; min s_lo of owed {sl}')
    if left:
        print('     e.g.', left[:3])
    del D
D = json.load(open(CODE + r'\runs\KQF.json'))
ow = [tuple(b) for r in D['results'] for b in r['left']] + [tuple(f[0]) for r in D['results'] for f in r['fails']]
inside = [b for b in ow if b[2] >= 0.88 and b[4] >= 1.29]
print(f'   KQF base: owed {len(ow)}, inside K_lw {len(inside)}, still owed {len(ow) - len(inside)}')
try:
    B = json.load(open(CODE + r'\runs\KQF_r1_boxes.json'))['boxes']
    print(f'   KQF_r1_boxes.json: {len(B)} boxes; inside K_lw {sum(1 for b in B if b[2] >= 0.88 and b[4] >= 1.29)}; '
          f'equal to the owed set outside K_lw: {set(map(tuple, B)) == set(ow) - set(inside)}')
except Exception as ex:
    print('   KQF_r1_boxes.json:', ex)
# 6. shell/30 vs face/30
random.seed(5); S = [1.0, 1.0, HP, 0.0025]; neq = 0; ncut = 0; bad = 0
for _ in range(20000):
    e0 = 1 - 10 ** random.uniform(-9, math.log10(4e-4)); e1 = min(1.0, e0 + 10 ** random.uniform(-10, -4))
    t0 = 1 - 10 ** random.uniform(-14, math.log10(1.7e-7)); t1 = min(1.0, t0 + 10 ** random.uniform(-15, -8))
    l0 = random.uniform(0.88, 0.99); w0 = random.uniform(1.29, 1.55)
    bx = (e0, e1, l0, min(1.0, l0 + random.uniform(1e-4, 0.1)), w0, min(HP, w0 + random.uniform(1e-4, 0.1)), t0, t1)
    a, c = splitrule.split_box(bx, S, 'shell/30'); a2, c2 = splitrule.split_box(bx, S, 'face/30')
    S2 = SF * SF * (1 - 5e-9); cut = (1 - bx[1]) < S2 and splitrule.mu_of(bx[7]) < S2 and (bx[0] < 1 - S2 < bx[1] or bx[6] < 1 / (1 + S2 * S2) < bx[7])
    ncut += cut
    if (a, c) != (a2, c2) and not cut:
        neq += 1
    diff = [i for i in range(8) if a[i] != bx[i] or c[i] != bx[i]]
    if not (len(diff) == 2 and diff[0] % 2 == 0 and diff[1] == diff[0] + 1 and a[diff[1]] == c[diff[0]] and bx[diff[0]] <= a[diff[1]] <= bx[diff[1]]):
        bad += 1
print(f'6. shell/30 vs face/30 on 20000 random boxes: boxes where face/30 cuts (1): {ncut}; different splits off the cut: {neq}; '
      f'shell splits not tiling: {bad}')

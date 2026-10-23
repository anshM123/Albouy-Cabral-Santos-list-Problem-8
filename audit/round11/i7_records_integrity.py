"""AUDIT-11 i7: integrity of records that went through the concurrent-writer episode (R11-3) and the --resume repair.
For every parseable record: replay of the split tree (header rule and scale), and comparison of the replayed leaf-code
histogram with the record's own 'hist' (kq_run work(): hist[bound], 'o', 'p', 'w<W>'), nleaf, nwind, fails, left
and volume.  A line spliced from two records (head of one, tail of another) that still parses would almost surely
fail the histogram or the count checks.  Also: duplicates (same root) must be identical when both complete.
usage: python i7_records_integrity.py FILE[.json|.records.jsonl] RULE|- S0 S1 S2 S3   (rule/scale ignored for .json)"""
import sys, json, collections, math
sys.path.insert(0, r'<repo>\code')
import splitrule
fn = sys.argv[1]
if fn.endswith('.json'):
    D = json.load(open(fn)); H = D['header']; recs = D['results']; rule = H.get('split_rule'); S = H['scale']; bad_lines = 0
else:
    rule = None if sys.argv[2] == '-' else sys.argv[2]; S = [float(v) for v in sys.argv[3:7]]
    recs = []; bad_lines = 0
    for line in open(fn):
        try:
            recs.append(json.loads(line))
        except Exception:
            bad_lines += 1


def vol(b):
    return math.prod(b[2 * i + 1] - b[2 * i] for i in range(4))


def split(bx):
    if rule is not None:
        return splitrule.split_box(bx, S, rule)
    w = [(bx[2 * i + 1] - bx[2 * i]) / S[i] for i in range(4)]
    i = max(range(4), key=lambda j: w[j]); mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)


probs = collections.Counter(); ex = []
by_root = collections.defaultdict(list)
for r in recs:
    root = tuple(r['box']); by_root[root].append(r)
    stack = [root]; h = collections.Counter(); nl = 0; nd = 0; fl = []; v = 0.0; ok = True
    for ch in r['tree']:
        if not stack:
            probs['underflow'] += 1; ok = False; break
        bx = stack.pop()
        if ch == '+':
            stack.extend(split(bx)); continue
        v += vol(bx)
        if ch == '!':
            fl.append(list(bx)); continue
        nl += 1
        if ch in 'op':
            h[ch] += 1
        elif '1' <= ch <= '6':
            h['w' + ch] += 1; nd += 1
        elif 'a' <= ch <= 'l':
            h[str(ord(ch) - ord('a'))] += 1
        else:
            probs['bad code ' + ch] += 1; ok = False
    left = [list(b) for b in stack]
    checks = {'nleaf': nl == r['nleaf'], 'hist': dict(h) == {str(k): v_ for k, v_ in r.get('hist', {}).items()},
              'nwind': ('nwind' not in r) or nd == r['nwind'], 'wL': nd == len(r.get('wL', [])),
              'fails': [f[0] for f in r.get('fails', [])] == fl,
              'left': (not stack) if r['done'] else left == r.get('left', []),
              'volume': abs(v + sum(vol(b) for b in left if not r['done']) - vol(root)) <= 1e-9 * vol(root)}
    for k, good in checks.items():
        if not good:
            probs[k] += 1
            if len(ex) < 5:
                ex.append((k, root[:2], r['done'], len(r['tree'])))
dups = {k: v for k, v in by_root.items() if len(v) > 1}
dup_diff = sum(1 for v in dups.values() if len(set(json.dumps(x, sort_keys=True) for x in v if x['done'])) > 1)
print(f'{fn.split(chr(92))[-1]}: records {len(recs)}, unparsable lines {bad_lines}, roots {len(by_root)}, duplicated roots '
      f'{len(dups)} (complete duplicates that differ: {dup_diff}); rule {rule}; problems {dict(probs)} {ex}')

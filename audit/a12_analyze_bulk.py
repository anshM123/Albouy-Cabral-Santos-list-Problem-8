"""AUDIT a12: analyse an nb_run output (final JSON or the running .fails.jsonl): completeness, failures and their
location relative to the known degenerate sets:
  curve6   : w = pi/2, e = l^2, t -> 1  (6-fold root of f_1 at u = pi)
  faces    : w = 0, w = pi/2, l = 0, l = 1, t = 0, t = 1, eta at the domain ends.
usage: python a12_analyze_bulk.py FILE [domain eta_lo eta_hi]"""
import json, math, sys
import numpy as np
fn = sys.argv[1]
recs = []
if fn.endswith('.jsonl'):
    for line in open(fn):
        if line.strip():
            recs.append(json.loads(line))
    print(f"{fn}: {len(recs)} initial boxes with failures or incomplete (running log)")
else:
    d = json.load(open(fn)); recs = d['results']
    print(f"{fn}: domain {d['domain']}, n0 {d['n0']}, initial boxes {len(recs)}, complete {sum(r['done'] for r in recs)}")
    nb = int(np.prod(d['n0']))
    print("   all initial boxes present exactly once:", len(recs) == nb and len({tuple(r['box']) for r in recs}) == nb)
    print("   certified leaves:", sum(r['leaves'][0] if 'leaves' in r else r['ok'] for r in recs))
F = [f for r in recs for f in r['fails']]
inc = [r for r in recs if not r['done']]
print("   failed leaves:", len(F), "; incomplete initial boxes:", len(inc))
if F:
    B = np.array([f[0] for f in F]); tot = np.array([f[1] if not isinstance(f[1], list) else f[1][0] for f in F])
    mid = np.stack([(B[:, 2*i] + B[:, 2*i+1]) / 2 for i in range(4)], 1)
    e = np.sin(mid[:, 0]); l = mid[:, 1]; w = mid[:, 2]; t = mid[:, 3]
    print("   failed totals:", dict(zip(*[x.tolist() for x in np.unique(tot, return_counts=True)])))
    c6 = (np.abs(w - math.pi/2) < 0.01) & (np.abs(e - l*l) < 0.01) & (t > 0.99)
    print("   on curve6 (|w-pi/2|<.01, |e-l^2|<.01, t>.99): %d" % c6.sum())
    rest = ~c6
    if rest.any():
        # cluster the rest by rounding the midpoint
        keys = {}
        for i in np.flatnonzero(rest):
            k = (round(e[i], 2), round(l[i], 2), round(w[i], 2), round(t[i], 2))
            keys[k] = keys.get(k, 0) + 1
        print("   other failure clusters (e, l, w, t rounded to .01): count")
        for k, v in sorted(keys.items(), key=lambda kv: -kv[1])[:20]:
            print("      ", k, v)
for r in inc[:10]:
    print("   incomplete initial box:", np.round(r['box'], 4).tolist())

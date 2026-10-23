"""AUDIT-11 i0: headers and leaf-code sets of the round-11 runs (read-only)."""
import sys, json, os, glob, collections
sys.path.insert(0, r'<repo>\code')
from runio import load_json, exists
os.chdir(r'<repo>\code')
for fam in sys.argv[1:]:
    fs = [f for f in [fam + '.json'] + sorted(glob.glob(fam + '.r*.json')) if exists(f)]
    for fn in fs:
        try:
            D = load_json(fn)
        except Exception as ex:
            print(fn, 'cannot load:', ex); continue
        H = D['header']
        keys = {k: v for k, v in H.items() if k not in ('skipped', 'exclude')}
        print('==', fn, json.dumps(keys, default=str)[:3000])
        cs = collections.Counter(); nd = 0; nfail = 0; nleft = 0
        for r in D['results']:
            cs.update(r['tree']); nd += r['done']; nfail += len(r.get('fails', [])); nleft += len(r.get('left', []))
        print('   records', len(D['results']), 'done', nd, 'fails', nfail, 'left', nleft, 'codes', dict(sorted(cs.items())))

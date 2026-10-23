"""AUDIT-7 g2b: g2's replay applied to the completed records of the PARTIAL bulk run (records.jsonl snapshot; the grid
check is expected to fail until the run has finished)."""
import sys, json, os, math
sys.dont_write_bytecode = True
import g2_replay_all as g2
snap = sys.argv[1]
recs = []
for line in open(snap):
    try:
        r = json.loads(line)
    except Exception:
        continue
    if r.get('done'):
        recs.append(r)
dom = (0.20033740617957885, 1.2535563099909826, 0, 1, 0, g2.HP, 0, 1)
hdr = dict(engine='a', domain=list(dom), n0=[8, 10, 10, 10], exclude=[], skipped=[], scale=[math.pi / 2, 1.0, math.pi / 2, 1.0])
g2.load = lambda path: dict(header=hdr, results=recs) if path.endswith('bulk.json') else None
st, owed, _ = g2.replay_family('bulk', 'runs/bulk', 'bb', 'a', dom, (8, 10, 10, 10), [])
st['problems'] = [p for p in st['problems'] if not p.startswith('initial boxes != grid')] + ['(grid check skipped: partial)']
print(json.dumps(st, default=str))

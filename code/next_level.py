"""Close a continuation level of a run family and prepare the next one.

usage: python next_level.py KIND NAME LEVEL
  KIND  'kq' (kq_run.py, files runs/NAME.rLEVEL.json) or 'c3' (c3/c3_run.py, files c3/runs/NAME.rLEVEL.json)
  Reads the run file of level LEVEL (base file for LEVEL 0), lists its open boxes (left stacks of incomplete
  records and failed leaves) and writes them to <dir>/NAME_r<LEVEL+1>_boxes.json (key "boxes").  The file must have
  been closed first, i.e. written by the driver after every box returned (run the driver with TMAX = -1 and --resume
  to record the boxes not yet processed as open records).  Prints the number of open boxes."""
import json
import sys

kind, name, level = sys.argv[1], sys.argv[2], int(sys.argv[3])
d = 'runs' if kind == 'kq' else 'c3/runs'
fn = '%s/%s.json' % (d, name) if level == 0 else '%s/%s.r%d.json' % (d, name, level)
D = json.load(open(fn))
rs = D['results']
op = [b for r in rs for b in r['left']] + [x[0] for r in rs for x in r['fails']]
# face charts (c3_run FE, FM, FIP, FIM): only s <= SF is needed (the face neighbourhood has s = max(sqrt(eps), sqrt(mu))
# <= SF, and the gap tau in [YF s, YF SF] has sigma = tau/YF <= SF); open boxes with s_lo > SF are dropped (check_all.py
# excuses them by the same rule)
if kind == 'c3' and name in ('FE', 'FM', 'FIP', 'FIM'):
    n_ = len(op); op = [b for b in op if not b[0] > 0.02]
    print('face chart: dropped', n_ - len(op), 'open boxes with s_lo > SF = 0.02')
# KQF: boxes inside the corner K_lw = {l >= 0.88, w >= 1.29} are not needed (covered by KQF2 and KQS; check_all.py)
if kind == 'kq' and name == 'KQF':
    n_ = len(op); op = [b for b in op if not (b[2] >= 0.88 and b[4] >= 1.29)]
    print('KQF: dropped', n_ - len(op), 'open boxes inside the corner K_lw')
out = '%s/%s_r%d_boxes.json' % (d, name, level + 1)
json.dump(dict(boxes=op), open(out, 'w'))
print(fn, 'records', len(rs), 'complete', sum(r['done'] for r in rs), 'open boxes', len(op), '->', out)

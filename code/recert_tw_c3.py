"""Re-certify EVERY tau-winding leaf ('7') of the c3 chart runs (base and continuation files) from its stored contour
(tau_wind.check_tw_chart), replaying each split tree with the file's header scale (audit R13-1).
usage: python recert_tw_c3.py CHART [CHART ...]      (run from code/; files c3/runs/CHART.json and CHART.r*.json)"""
import glob
import sys
import numpy as np
from runio import load_json, exists
import tau_wind

sys.path.insert(0, 'c3')
tot = 0; bad = 0
for chart in sys.argv[1:]:
    files = ['c3/runs/%s.json' % chart] + sorted(set(f[:-3] if f.endswith('.gz') else f for f in
                                                     glob.glob('c3/runs/%s.r*.json' % chart) + glob.glob('c3/runs/%s.r*.json.gz' % chart)),
                                                 key=lambda s: int(s.split('.r')[-1].split('.')[0]))
    for fn in files:
        if not exists(fn):
            continue
        D = load_json(fn); H = D['header']; scale = H['scale']; n = 0; nb = 0
        for rec in D['results']:
            stack = [tuple(rec['box'])]; ti = 0
            for ch in rec['tree']:
                bx = stack.pop()
                if ch == '+':
                    w = [(bx[2 * i + 1] - bx[2 * i]) / scale[i] for i in range(4)]
                    i = int(np.argmax(w)); mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
                    a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
                    stack.extend([tuple(a), tuple(c)])
                    continue
                if ch == '7':
                    which, cre, cim, r = rec['tw'][ti]; ti += 1; n += 1
                    if which != 0 or not tau_wind.check_tw_chart(chart, bx, cre, cim, r):
                        nb += 1
                        print('TW re-certification problem', fn, bx, rec['tw'][ti - 1], flush=True)
            if ti != len(rec.get('tw', [])):
                nb += 1
                print('TW leaf count mismatch', fn, rec['box'], flush=True)
        if n:
            print(fn, 'tau-winding leaves re-certified', n, 'problems', nb, flush=True)
        tot += n; bad += nb
print('TOTAL tau-winding leaves re-certified', tot, 'problems', bad)

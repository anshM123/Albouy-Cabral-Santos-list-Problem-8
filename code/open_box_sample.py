"""60-digit sample of the OPEN boxes of unfinished KQ families (not a proof): for random points of the boxes still
owed by a family (the boxes listed for its next continuation level), count the real roots of the original a-chart
polynomial with multiplicity (verify_mp.pol_count) and report the maximum.  usage: python open_box_sample.py N FILE ..."""
import json
import math
import random
import sys
from verify_mp import pol_count

random.seed(7)
n = int(sys.argv[1])
for fn in sys.argv[2:]:
    B = json.load(open(fn))['boxes']
    worst = 0; cnt = 0; hist = {}
    for _ in range(n):
        b = random.choice(B)
        x = [random.uniform(b[2 * i], b[2 * i + 1]) for i in range(4)]
        e, l, w, t = x
        if not (0 < t < 1 and e < 1):
            continue
        s = math.sqrt(t / (1 - t)) * (1 - e)
        N = pol_count(e, l, w, s)
        cnt += 1; worst = max(worst, N); hist[N] = hist.get(N, 0) + 1
    print(fn, 'boxes', len(B), 'points', cnt, 'max N', worst, 'hist', dict(sorted(hist.items())), flush=True)

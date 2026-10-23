"""AUDIT a7: where do the v5 cost-map failures sit?  Test the conjecture that they lie on the degenerate curve
omega = pi/2, e = l^2, t -> 1 (where (rho^2)' = 2 l^2 (1-l^2) sin u (1 + cos u) has a triple zero at u = pi)."""
import json, math, sys
import numpy as np
fn = sys.argv[1] if len(sys.argv) > 1 else '../code/costmap.json'
d = json.load(open(fn))
res = d['results']
F = [f for r in res for f in r['fails']]
print("failed leaves:", len(F))
if F:
    B = np.array([f[0] for f in F]); tot = np.array([f[1] for f in F])
    mid = np.stack([(B[:, 2*i] + B[:, 2*i+1]) / 2 for i in range(4)], 1)
    e = np.sin(mid[:, 0]); l = mid[:, 1]; w = mid[:, 2]; t = mid[:, 3]
    print("failed totals histogram:", dict(zip(*np.unique(tot, return_counts=True))))
    near_w = np.abs(w - np.pi/2) < 0.02; near_t = t > 0.999; near_c = np.abs(e - l*l) < 0.01
    print("fraction with |w - pi/2| < 0.02: %.3f; t > 0.999: %.3f; |e - l^2| < 0.01: %.3f; all three: %.3f" % (
        near_w.mean(), near_t.mean(), near_c.mean(), (near_w & near_t & near_c).mean()))
    other = ~(near_w & near_t & near_c)
    if other.any():
        print("failures NOT on the curve (e, l, w, t, total):")
        for i in np.flatnonzero(other)[:25]:
            print("   e=%.4f l=%.4f w=%.4f t=%.6f  tot=%d  widths=%s" % (e[i], l[i], w[i], t[i], tot[i],
                  np.round([B[i,1]-B[i,0], B[i,3]-B[i,2], B[i,5]-B[i,4], B[i,7]-B[i,6]], 6)))

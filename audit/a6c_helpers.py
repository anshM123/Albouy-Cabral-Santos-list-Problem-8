import math
import numpy as np
import nb_cert
from nbcoef_a import coef_val as cval_a, coef_abs as cabs_a
from nbport import box_coeffs_port, arc_bounds_port
def port_trace(chart, box, J=6):
    pieces = box_coeffs_port(chart, *box[:6], cval_a, cabs_a)
    M0 = 64
    left = np.arange(M0) * (2 * math.pi / M0); width = np.full(M0, 2 * math.pi / M0)
    for rd in range(9):
        um = left + width / 2; hh = width / 2
        pos = np.ones((J + 1, um.size), bool); neg = np.ones((J + 1, um.size), bool)
        hop = np.ones(um.size, bool)
        for t in (box[6], box[7]):
            bd = arc_bounds_port(*pieces, t, um, hh, J)
            pos &= bd['lo'] > 0; neg &= bd['hi'] < 0
            hop &= np.all(np.abs(bd['vc'][:J+1]) <= bd['rad'][:J+1], axis=0)
        sg = np.where(pos, 1, np.where(neg, -1, 0)).astype(np.int64)
        tot = nb_cert.dp_count(sg)
        if tot <= 10 or hop.any() or rd == 8:
            break
        sel = sg[0] == 0
        nl, nw = [], []
        for i in range(left.size):
            if sel[i]: nl += [left[i], left[i] + width[i]/2]; nw += [width[i]/2]*2
            else: nl.append(left[i]); nw.append(width[i])
        left, width = np.array(nl), np.array(nw)
    print(f"trace: round {rd}, total {tot}, hopeless arcs {int(hop.sum())}, arcs {left.size}")
    bad = np.flatnonzero(np.all(sg == 0, axis=0))
    print("  arcs with no certified derivative sign (u-mid):", np.round(um[bad], 4)[:20])
    near = np.flatnonzero(np.abs(np.angle(np.exp(1j*(um - math.pi)))) < 0.3)
    for i in near[::max(1, len(near)//12)]:
        print("   u=%.4f signs j=0..6: %s" % (um[i], sg[:, i].tolist()))
    near0 = np.flatnonzero(np.abs(np.angle(np.exp(1j*um))) < 0.2)
    for i in near0[::max(1, len(near0)//6)]:
        print("   u=%.4f signs j=0..6: %s" % (um[i], sg[:, i].tolist()))

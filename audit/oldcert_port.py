"""port of the OLD (pre-23:36) nb_cert.certify_box refinement loop (no stall rule), for regression tests."""
import math
import numpy as np
import nb_cert
from nbcoef_a import coef_val as cval_a, coef_abs as cabs_a
from nbcoef_f import coef_val as cval_f, coef_abs as cabs_f
from nbport import box_coeffs_port, arc_bounds_port
J = 6


def port_certify(chart, box):
    """python replica of certify_box's refinement loop, returns (ok, total, left, width, sg)"""
    cval, cabs = (cval_a, cabs_a) if chart == 0 else (cval_f, cabs_f)
    pieces = box_coeffs_port(chart, *box[:6], cval, cabs)
    M0, rounds = 64, 8
    left = np.arange(M0) * (2 * math.pi / M0); width = np.full(M0, 2 * math.pi / M0)
    total = 999
    for rd in range(rounds + 1):
        um = left + width / 2; hh = width / 2
        pos = np.ones((J + 1, um.size), bool); neg = np.ones((J + 1, um.size), bool)
        hopeless = np.ones(um.size, bool)
        for t in (box[6], box[7]):
            bd = arc_bounds_port(*pieces, t, um, hh, J)
            pos &= bd['lo'] > 0; neg &= bd['hi'] < 0
            hopeless &= np.all(np.abs(bd['vc'][:J + 1]) <= bd['rad'][:J + 1], axis=0)
        sg = np.where(pos, 1, np.where(neg, -1, 0)).astype(np.int64)
        total = nb_cert.dp_count(sg)
        if total <= 10:
            return True, total, left, width, sg
        if rd == rounds or hopeless.any():
            break
        sel = sg[0] == 0
        if not sel.any():
            break
        nl = []; nw = []
        for i in range(left.size):
            if sel[i]:
                nl += [left[i], left[i] + width[i] / 2]; nw += [width[i] / 2, width[i] / 2]
            else:
                nl.append(left[i]); nw.append(width[i])
        left = np.array(nl); width = np.array(nw)
    return False, total, left, width, sg



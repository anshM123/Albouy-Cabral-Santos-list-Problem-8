"""Away-from-aphelion count for the C3 parameter box: real zeros of the f-chart polynomial on the arc
f in [pi + phi0, 3 pi - phi0] (mod 2 pi), uniformly over a box (e, l, w) x [t0, t1], t = 1/(1 + mu^2).
Same enclosures as nb_corner.away_count (nb_cert.box_coeffs / arc_signs, linear dp), but with an explicit target,
no early exit at 6, and an arc cap.  Arcs padded by ARC_PAD; the range is widened by PAD on both ends."""
import math
import numpy as np
from numba import njit
from nb_cert import box_coeffs, arc_signs, ARC_PAD
from nb_corner import dp_linear

PAD = 1e-12


@njit(cache=True)
def away_count_t(e0, e1, l0, l1, w0, w1, t0, t1, phi0, M0, rounds, J, target, cap):
    cmid, crad, dcen, dwid, cmax, rr = box_coeffs(1, e0, e1, l0, l1, w0, w1)
    a = math.pi + phi0 - PAD; L = 2 * math.pi - 2 * phi0 + 2 * PAD
    left = a + np.arange(M0) * (L / M0)
    width = np.full(M0, L / M0)
    best = 999.0
    for rd in range(rounds + 1):
        sg, hopeless = arc_signs(cmid, crad, dcen, dwid, cmax, rr, t0, t1, left + width / 2, width / 2 + ARC_PAD, J)
        best = min(best, dp_linear(sg))      # every round's bound is valid; keep the smallest
        if best <= target:
            return True, best
        if rd == rounds or hopeless or left.size > cap:
            break
        nsel = 0
        for i in range(left.size):
            if sg[0, i] == 0:
                nsel += 1
        if nsel == 0:
            break
        nl = np.empty(left.size + nsel); nw = np.empty(left.size + nsel)
        m = 0
        for i in range(left.size):
            if sg[0, i] == 0:
                nl[m] = left[i]; nw[m] = width[i] / 2; m += 1
                nl[m] = left[i] + width[i] / 2; nw[m] = width[i] / 2; m += 1
            else:
                nl[m] = left[i]; nw[m] = width[i]; m += 1
        left = nl; width = nw
    return False, best

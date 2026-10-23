"""Second-order parameter enclosure for the f-chart (chart 1 of nb_cert), used by the away count of engine KQ.

nb_cert bounds the variation of the f-chart polynomial over a parameter box x = (e, l, w) by the mean-value form
    |P(u; x) - P(u; x_c)| <= sum_v r_v (|dP/dx_v(u; x_c)| + w_v(u)),
where w_v bounds the variation of the gradient over the box by absolute-value polynomials (no cancellation).  Here the
same quantity is bounded by the second-order Taylor form
    |P(u; x) - P(u; x_c)| <= sum_v r_v |G_v(u)| + 1/2 sum_{v,w} r_v r_w |H_vw(u)| + sum_k |T_k(u)| R_k,
with G and H the gradient and the Hessian of P^(j)(u; .) at the box centre (exact chain rule from the derivatives of
the coefficient polynomials in (e, p, q), p = l sin w, q = l cos w), and R_k a bound for the third-order Taylor
remainder of coefficient k plus the evaluation errors of its centre derivatives.  The smaller of the two bounds is
used (both are valid).  Everything else (arcs, t-affine endpoints, Taylor tail in the anomaly, float padding) is as
in nb_cert.arc_signs; the extra sums are padded by 1e-12 of their absolute terms.

Third-order remainder of c(x) = C(e, l sin w, l cos w) along x_c + tau d, tau in [0, 1] (Lagrange form):
    |R| <= (1/6) sup |g'''|,  g''' = sum C_ijk g1_i g1_j g1_k + 3 sum C_ij g2_i g1_j + sum C_i g3_i,
with g1 = (de, dl s + l dw c, dl c - l dw s), |g1_p|, |g1_q| <= a1 = r_l + lmax r_w, |g1_e| <= r_e,
g2 = (0, 2 dl dw c - l dw^2 s, -2 dl dw s - l dw^2 c), |g2_p|, |g2_q| <= a2 = 2 r_l r_w + lmax r_w^2,
g3 = (0, -3 dl dw^2 s - l dw^3 c, -3 dl dw^2 c + l dw^3 s), |g3_p|, |g3_q| <= a3 = 3 r_l r_w^2 + lmax r_w^3,
and sup |C_i|, |C_ij|, |C_ijk| over the box bounded by the absolute-value polynomials at (|e_c| + r_e, |p_c| + dp,
|q_c| + dp), dp = r_l + lmax r_w (padded)."""
import math
import numpy as np
from numba import njit
from nb_cert import box_coeffs, ARC_PAD
from nbcoef_f2 import coef_val2, coef_abs2
from nb_corner import dp_linear
from wind_cert import G_centre, arg_in_branch

PAD = 1e-12
# index tables of nbcoef_f2: second derivatives (ee, ep, eq, pp, pq, qq) -> slots 4..9, third -> slots 10..19
SEC = np.array([[0, 0], [0, 1], [0, 2], [1, 1], [1, 2], [2, 2]])
THI = np.array([[0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 1, 1], [0, 1, 2], [0, 2, 2], [1, 1, 1], [1, 1, 2], [1, 2, 2],
                [2, 2, 2]])
MULT3 = np.array([1.0, 3.0, 3.0, 3.0, 6.0, 3.0, 1.0, 3.0, 3.0, 1.0])     # ordered triples per sorted triple


@njit(cache=True)
def sec_slot(i, j):
    a = min(i, j); b = max(i, j)
    if a == 0:
        return 4 + b                 # ee, ep, eq
    if a == 1:
        return 6 + b                 # pp (7), pq (8)
    return 9                         # qq


@njit(cache=True)
def box_coeffs2(x0lo, x0hi, llo, lhi, wlo, whi):
    """f-chart data of nb_cert.box_coeffs plus hcen[6, f, part, k] (Hessian in (e, l, w) at the centre, order
    ee, el, ew, ll, lw, ww) and rem[f, part, k] (third-order remainder + evaluation errors of the centre derivatives)."""
    cmid, crad, dcen, dwid, cmax, rr = box_coeffs(1, x0lo, x0hi, llo, lhi, wlo, whi)
    xc = 0.5 * (x0lo + x0hi); lc = 0.5 * (llo + lhi); wc = 0.5 * (wlo + whi)
    sw, cw = math.sin(wc), math.cos(wc)
    pc, qc = lc * sw, lc * cw
    re_, rl, rw = rr[0], rr[1], rr[2]
    lmax = lc + rl
    a1 = (rl + lmax * rw) * (1 + 1e-12) + 1e-15
    a2 = (2 * rl * rw + lmax * rw * rw) * (1 + 1e-12)
    a3 = (3 * rl * rw * rw + lmax * rw * rw * rw) * (1 + 1e-12)
    val = np.zeros((10, 2, 2, 7)); ab = np.zeros((20, 2, 2, 7))
    coef_val2(xc, pc, qc, val)
    coef_abs2(abs(xc) + re_, abs(pc) + a1, abs(qc) + a1, ab)
    A = np.array([re_, a1, a1])
    G2 = np.array([0.0, a2, a2])
    G3 = np.array([0.0, a3, a3])
    hcen = np.empty((6, 2, 2, 7)); rem = np.empty((2, 2, 7))
    for f in range(2):
        for pt in range(2):
            for k in range(7):
                Cp = val[2, f, pt, k]; Cq = val[3, f, pt, k]
                Cee = val[4, f, pt, k]; Cep = val[5, f, pt, k]; Ceq = val[6, f, pt, k]
                Cpp = val[7, f, pt, k]; Cpq = val[8, f, pt, k]; Cqq = val[9, f, pt, k]
                hcen[0, f, pt, k] = Cee
                hcen[1, f, pt, k] = Cep * sw + Ceq * cw
                hcen[2, f, pt, k] = lc * (Cep * cw - Ceq * sw)
                hcen[3, f, pt, k] = Cpp * sw * sw + 2 * Cpq * sw * cw + Cqq * cw * cw
                hcen[4, f, pt, k] = lc * (Cpp * sw * cw + Cpq * (cw * cw - sw * sw) - Cqq * sw * cw) + (Cp * cw - Cq * sw)
                hcen[5, f, pt, k] = lc * lc * (Cpp * cw * cw - 2 * Cpq * sw * cw + Cqq * sw * sw) - lc * (Cp * sw + Cq * cw)
                # third-order remainder (1/6) sup |g'''|
                t3 = 0.0
                for m in range(10):
                    t3 += MULT3[m] * ab[10 + m, f, pt, k] * A[THI[m, 0]] * A[THI[m, 1]] * A[THI[m, 2]]
                t2 = 0.0
                for i in range(1, 3):
                    for j in range(3):
                        t2 += ab[sec_slot(i, j), f, pt, k] * G2[i] * A[j]
                t1 = (ab[2, f, pt, k] + ab[3, f, pt, k]) * a3
                r3 = (t3 + 3.0 * t2 + t1) / 6.0
                # evaluation errors of the centre gradient and Hessian (chain rule included), generously padded
                s1 = ab[1, f, pt, k] + ab[2, f, pt, k] + ab[3, f, pt, k]
                s2 = 0.0
                for m in range(6):
                    s2 += ab[4 + m, f, pt, k]
                e1 = 1e-13 * (s1 * (1.0 + lmax)) * (re_ + rl + rw)
                e2 = 1e-13 * ((s2 + s1) * (1.0 + lmax) * (1.0 + lmax)) * (re_ + rl + rw) * (re_ + rl + rw)
                rem[f, pt, k] = (r3 + e1 + e2) * (1 + 1e-12)
    return cmid, crad, dcen, dwid, cmax, rr, hcen, rem


@njit(cache=True)
def enclose_point2(c, cr, dc, dw, hc, rm, DA, HA, RA, qw, rr, u, J, mode):
    """value vc[j] of P^(j)(u; x_c) and radius rad[j] >= |P^(j)(u; x) - vc[j]| for all x in the box, j = 0..J+1
    (one t-endpoint; c, cr, dc, dw, hc, rm are the t-combined data).  mode 0: the smaller of the two complete radii
    (used by the certificates); 1: the second-order radius only; 2: the mean-value radius only (tests).  The mean-value
    radius is computed with the same operations and the same padding as nb_cert.arc_signs."""
    vc = np.empty(J + 2); rad = np.empty(J + 2)
    for j in range(J + 2):
        s = 0.0; mag = 0.0; rd = 0.0
        g = np.zeros(3); w = np.zeros(3); H = np.zeros(6)
        aab = 0.0; aab_new = 0.0; pabs = 0.0; rsum = 0.0
        for k in range(7):
            if j > 0 and k == 0:
                continue
            kj = float(k) ** j
            # trig-argument budget of the mean-value radius: exactly as nb_cert.arc_signs
            aab += kj * (abs(c[0, k]) + abs(c[1, k]) + cr[0, k] + cr[1, k]
                         + rr[0] * (abs(dc[0, 0, k]) + abs(dc[0, 1, k]) + dw[0, 0, k] + dw[0, 1, k])
                         + rr[1] * (abs(dc[1, 0, k]) + abs(dc[1, 1, k]) + dw[1, 0, k] + dw[1, 1, k])
                         + rr[2] * (abs(dc[2, 0, k]) + abs(dc[2, 1, k]) + dw[2, 0, k] + dw[2, 1, k]))
            T = math.cos(k * u + j * math.pi / 2) * kj
            U = math.sin(k * u + j * math.pi / 2) * kj
            aT = abs(T); aU = abs(U)
            s += c[0, k] * T + c[1, k] * U
            mag += abs(c[0, k]) * aT + abs(c[1, k]) * aU
            rd += cr[0, k] * aT + cr[1, k] * aU
            hk = 0.0
            for v in range(3):
                g[v] += dc[v, 0, k] * T + dc[v, 1, k] * U
                w[v] += dw[v, 0, k] * aT + dw[v, 1, k] * aU
                pabs += rr[v] * (DA[v, 0, k] * aT + DA[v, 1, k] * aU)
            for m in range(6):
                H[m] += hc[m, 0, k] * T + hc[m, 1, k] * U
                pabs += qw[m] * (HA[m, 0, k] * aT + HA[m, 1, k] * aU)
                hk += qw[m] * (HA[m, 0, k] + HA[m, 1, k])
            rsum += rm[0, k] * aT + rm[1, k] * aU
            pabs += RA[0, k] * aT + RA[1, k] * aU
            # trig-argument budget of the second-order radius: all magnitudes entering s, rd, g, H and the remainder
            aab_new += kj * (abs(c[0, k]) + abs(c[1, k]) + cr[0, k] + cr[1, k]
                             + rr[0] * (abs(dc[0, 0, k]) + abs(dc[0, 1, k]))
                             + rr[1] * (abs(dc[1, 0, k]) + abs(dc[1, 1, k]))
                             + rr[2] * (abs(dc[2, 0, k]) + abs(dc[2, 1, k]))
                             + hk + RA[0, k] + RA[1, k])
        # mean-value radius (nb_cert.arc_signs: prad = sum_v r_v (|g_v| + w_v), then the same float padding)
        old = rr[0] * (abs(g[0]) + w[0]) + rr[1] * (abs(g[1]) + w[1]) + rr[2] * (abs(g[2]) + w[2])
        t_old = rd + old
        t_old = t_old + 1e-12 * (mag + t_old) + 2e-14 * aab
        # second-order radius
        first = rr[0] * abs(g[0]) + rr[1] * abs(g[1]) + rr[2] * abs(g[2])
        second = 0.0
        for m in range(6):
            second += qw[m] * abs(H[m])
        t_new = rd + first + second + rsum + 1e-12 * pabs
        t_new = t_new + 1e-12 * (mag + t_new) + 2e-14 * aab_new
        tot = min(t_old, t_new) if mode == 0 else (t_new if mode == 1 else t_old)
        vc[j] = s; rad[j] = tot
    return vc, rad


@njit(cache=True)
def arc_signs2(cmid, crad, dcen, dwid, cmax, rr, hcen, rem, t0, t1, um, hh, J):
    M = um.size
    pos = np.ones((J + 1, M), dtype=np.bool_)
    neg = np.ones((J + 1, M), dtype=np.bool_)
    hopeless = np.ones(M, dtype=np.bool_)
    # weights of the quadratic form bound: 1/2 r_v^2 for v = w, r_v r_w for v < w (order ee, el, ew, ll, lw, ww)
    qw = np.array([0.5 * rr[0] * rr[0], rr[0] * rr[1], rr[0] * rr[2], 0.5 * rr[1] * rr[1], rr[1] * rr[2],
                   0.5 * rr[2] * rr[2]])
    # absolute magnitudes (both t-parts) of the centre gradient, Hessian and remainder data: the extra sums of the
    # second-order bound are padded by 1e-12 of these (t-combination, summation and product roundings)
    DA = np.abs(dcen[:, 0]) + np.abs(dcen[:, 1])
    HA = np.abs(hcen[:, 0]) + np.abs(hcen[:, 1])
    RA = rem[0] + rem[1]
    for ti in range(2):
        t = t0 if ti == 0 else t1
        c = (1 - t) * cmid[0] - t * cmid[1]
        cr = (1 - t) * crad[0] + t * crad[1] + 1e-15 * (np.abs(cmid[0]) + np.abs(cmid[1]))
        dc = (1 - t) * dcen[:, 0] - t * dcen[:, 1]
        dw = (1 - t) * dwid[:, 0] + t * dwid[:, 1] + 1e-15 * (np.abs(dcen[:, 0]) + np.abs(dcen[:, 1]))
        hc = (1 - t) * hcen[:, 0] - t * hcen[:, 1]
        rm = (1 - t) * rem[0] + t * rem[1]
        cm = (1 - t) * cmax[0] + t * cmax[1]
        for i in range(M):
            u = um[i]; h = hh[i]
            vc, rad = enclose_point2(c, cr, dc, dw, hc, rm, DA, HA, RA, qw, rr, u, J, 0)
            for j in range(J + 1):
                if abs(vc[j]) > rad[j]:
                    hopeless[i] = False
                Bj2 = 0.0
                for k in range(1, 7):
                    Bj2 += float(k) ** (j + 2) * (cm[0, k] + cm[1, k])
                Bj2 *= 1 + 1e-12
                m1 = abs(vc[j + 1]) + rad[j + 1]
                tail = rad[j] + h * m1 + 0.5 * h * h * Bj2
                tail = tail + 1e-15 * (abs(vc[j]) + tail)          # rounding of the final sums
                lo = vc[j] - tail
                hi = vc[j] + tail
                if not (lo > 0):
                    pos[j, i] = False
                if not (hi < 0):
                    neg[j, i] = False
    sg = np.zeros((J + 1, M), dtype=np.int64)
    for j in range(J + 1):
        for i in range(M):
            if pos[j, i]:
                sg[j, i] = 1
            elif neg[j, i]:
                sg[j, i] = -1
    any_hopeless = False
    for i in range(M):
        if hopeless[i]:
            any_hopeless = True
    return sg, any_hopeless


@njit(cache=True)
def away_count_t2(e0, e1, l0, l1, w0, w1, t0, t1, phi0, M0, rounds, J, target, cap):
    """as c3_away.away_count_t (zeros of the f-chart polynomial on f in [pi + phi0, 3 pi - phi0], uniformly over the
    box), with the second-order enclosure of arc_signs2."""
    cmid, crad, dcen, dwid, cmax, rr, hcen, rem = box_coeffs2(e0, e1, l0, l1, w0, w1)
    a = math.pi + phi0 - PAD; L = 2 * math.pi - 2 * phi0 + 2 * PAD
    left = a + np.arange(M0) * (L / M0)
    width = np.full(M0, L / M0)
    best = 999.0
    for rd in range(rounds + 1):
        sg, hopeless = arc_signs2(cmid, crad, dcen, dwid, cmax, rr, hcen, rem, t0, t1, left + width / 2,
                                  width / 2 + ARC_PAD, J)
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


def away_count_arc2(e0, e1, l0, l1, w0, w1, t0, t1, a0, L0, M0, rounds, J, target, cap):
    """as away_count_t2, but on the arc f in [a0, a0 + L0] (padded by PAD at both ends): zeros of the f-chart polynomial
    there, uniformly over the box (kq_run --peri uses two such arcs between the aphelion and perihelion parts)."""
    cmid, crad, dcen, dwid, cmax, rr, hcen, rem = box_coeffs2(e0, e1, l0, l1, w0, w1)
    a = a0 - PAD; L = L0 + 2 * PAD
    left = a + np.arange(M0) * (L / M0)
    width = np.full(M0, L / M0)
    best = 999.0
    for rd in range(rounds + 1):
        sg, hopeless = arc_signs2(cmid, crad, dcen, dwid, cmax, rr, hcen, rem, t0, t1, left + width / 2,
                                  width / 2 + ARC_PAD, J)
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


# ---------------------------------------------------------------------------------------------------------------------
# Winding certificate (wind_cert) with the second-order enclosure: comp_bounds2 / labels2 / wind_box2 are
# wind_cert.comp_bounds / labels / wind_box with the radius min(mean-value bound, second-order bound) as above.

@njit(cache=True)
def comp_bounds2(c, cr, dc, dw, hc, rm, DA, HA, RA, qw, cm, rr, wgt, phase, u, h):
    """enclosure [lo, hi] over the arc [u - h, u + h] and the parameter box of
       S(th) = sum_k wgt_k (A_k cos(k th + phase) + B_k sin(k th + phase))   (A, B = coefficient data c[0], c[1]).
    The radius of each derivative is the smaller of the mean-value radius (computed exactly as wind_cert.comp_bounds)
    and the second-order radius, each with its own complete padding."""
    vals = np.empty(2); rads = np.empty(2)
    for j in range(2):
        s = 0.0; mag = 0.0; rd = 0.0; aab = 0.0; aab_new = 0.0; pabs = 0.0; rsum = 0.0
        g = np.zeros(3); w = np.zeros(3); H = np.zeros(6)
        for k in range(7):
            if j > 0 and k == 0:
                continue
            kj = wgt[k] * float(k) ** j
            T = math.cos(k * u + phase + j * math.pi / 2) * kj
            U = math.sin(k * u + phase + j * math.pi / 2) * kj
            aT = abs(T); aU = abs(U)
            s += c[0, k] * T + c[1, k] * U
            mag += abs(c[0, k]) * aT + abs(c[1, k]) * aU
            rd += cr[0, k] * aT + cr[1, k] * aU
            ab_k = abs(c[0, k]) + abs(c[1, k]) + cr[0, k] + cr[1, k]
            ab_n = abs(c[0, k]) + abs(c[1, k]) + cr[0, k] + cr[1, k]
            for v in range(3):
                g[v] += dc[v, 0, k] * T + dc[v, 1, k] * U
                w[v] += dw[v, 0, k] * aT + dw[v, 1, k] * aU
                ab_k += rr[v] * (abs(dc[v, 0, k]) + abs(dc[v, 1, k]) + dw[v, 0, k] + dw[v, 1, k])
                ab_n += rr[v] * (abs(dc[v, 0, k]) + abs(dc[v, 1, k]))
                pabs += rr[v] * (DA[v, 0, k] * aT + DA[v, 1, k] * aU)
            for m in range(6):
                H[m] += hc[m, 0, k] * T + hc[m, 1, k] * U
                pabs += qw[m] * (HA[m, 0, k] * aT + HA[m, 1, k] * aU)
                ab_n += qw[m] * (HA[m, 0, k] + HA[m, 1, k])
            rsum += rm[0, k] * aT + rm[1, k] * aU
            pabs += RA[0, k] * aT + RA[1, k] * aU
            ab_n += RA[0, k] + RA[1, k]
            aab += abs(kj) * ab_k
            aab_new += abs(kj) * ab_n
        # mean-value radius, as wind_cert.comp_bounds
        prad = 0.0
        for v in range(3):
            prad += rr[v] * (abs(g[v]) + w[v])
        t_old = rd + prad
        t_old = t_old + 1e-12 * (mag + t_old) + 2e-14 * aab
        # second-order radius
        first = rr[0] * abs(g[0]) + rr[1] * abs(g[1]) + rr[2] * abs(g[2])
        second = 0.0
        for m in range(6):
            second += qw[m] * abs(H[m])
        t_new = rd + first + second + rsum + 1e-12 * pabs
        t_new = t_new + 1e-12 * (mag + t_new) + 2e-14 * aab_new
        vals[j] = s; rads[j] = min(t_old, t_new)
    B2 = 0.0
    for k in range(1, 7):
        B2 += abs(wgt[k]) * float(k) ** 2 * (cm[0, k] + cm[1, k])
    B2 *= 1 + 1e-12
    tail = rads[0] + h * (abs(vals[1]) + rads[1]) + 0.5 * h * h * B2
    tail = tail + 1e-15 * (abs(vals[0]) + tail)
    return vals[0] - tail, vals[0] + tail


@njit(cache=True)
def labels2(cmid, crad, dcen, dwid, cmax, rr, hcen, rem, t0, t1, L, um, hh):
    """label per arc: 1 Re>0, 2 Im>0, 3 Re<0, 4 Im<0, 0 none (certified for all params and t in [t0, t1])."""
    M = um.size
    wR = np.empty(7); wI = np.empty(7)
    for k in range(7):
        wR[k] = math.cosh(k * L) * (1 + 1e-12)
        wI[k] = math.sinh(k * L) * (1 + 1e-12)
    qw = np.array([0.5 * rr[0] * rr[0], rr[0] * rr[1], rr[0] * rr[2], 0.5 * rr[1] * rr[1], rr[1] * rr[2],
                   0.5 * rr[2] * rr[2]])
    DA = np.abs(dcen[:, 0]) + np.abs(dcen[:, 1])
    HA = np.abs(hcen[:, 0]) + np.abs(hcen[:, 1])
    RA = rem[0] + rem[1]
    lab = np.zeros(M, dtype=np.int64)
    okR_p = np.ones(M, dtype=np.bool_); okR_n = np.ones(M, dtype=np.bool_)
    okI_p = np.ones(M, dtype=np.bool_); okI_n = np.ones(M, dtype=np.bool_)
    for ti in range(2):
        t = t0 if ti == 0 else t1
        c = (1 - t) * cmid[0] - t * cmid[1]
        cr = (1 - t) * crad[0] + t * crad[1] + 1e-15 * (np.abs(cmid[0]) + np.abs(cmid[1]))
        dc = (1 - t) * dcen[:, 0] - t * dcen[:, 1]
        dw = (1 - t) * dwid[:, 0] + t * dwid[:, 1] + 1e-15 * (np.abs(dcen[:, 0]) + np.abs(dcen[:, 1]))
        hc = (1 - t) * hcen[:, 0] - t * hcen[:, 1]
        rm = (1 - t) * rem[0] + t * rem[1]
        cm = (1 - t) * cmax[0] + t * cmax[1]
        for i in range(M):
            loR, hiR = comp_bounds2(c, cr, dc, dw, hc, rm, DA, HA, RA, qw, cm, rr, wR, 0.0, um[i], hh[i])
            loI, hiI = comp_bounds2(c, cr, dc, dw, hc, rm, DA, HA, RA, qw, cm, rr, wI, -math.pi / 2, um[i], hh[i])
            if not (loR > 0):
                okR_p[i] = False
            if not (hiR < 0):
                okR_n[i] = False
            if not (loI > 0):
                okI_p[i] = False
            if not (hiI < 0):
                okI_n[i] = False
    for i in range(M):
        if okR_p[i]:
            lab[i] = 1
        elif okI_p[i]:
            lab[i] = 2
        elif okR_n[i]:
            lab[i] = 3
        elif okI_n[i]:
            lab[i] = 4
    return lab


@njit(cache=True)
def wind_box2(cmid, crad, dcen, dwid, cmax, rr, hcen, rem, t0, t1, L, M0, rounds, cap):
    """wind_cert.wind_box with labels2: returns (ok, W), W = number of zeros of Q in |z| < e^L, certified uniformly
    over the box.  Requires -15 <= L < 0."""
    if not (L < 0.0 and L >= -15.0):
        return False, -1
    left = np.arange(M0) * (2 * math.pi / M0)
    width = np.full(M0, 2 * math.pi / M0)
    for rd in range(rounds + 1):
        lab = labels2(cmid, crad, dcen, dwid, cmax, rr, hcen, rem, t0, t1, L, left + width / 2, width / 2 + ARC_PAD)
        nbad = 0
        for i in range(left.size):
            if lab[i] == 0:
                nbad += 1
        if nbad == 0:
            M = left.size
            vre = np.empty(M); vim = np.empty(M)
            for i in range(M):
                vre[i], vim[i] = G_centre(cmid, t0, L, left[i])
            tot = 0.0
            for i in range(M):
                j = i + 1 if i + 1 < M else 0
                a0 = arg_in_branch(lab[i], vre[i], vim[i])
                a1 = arg_in_branch(lab[i], vre[j], vim[j])
                if a0 > 50.0 or a1 > 50.0:
                    return False, -1
                tot += a1 - a0
            wind = int(round(tot / (2 * math.pi)))
            if abs(tot - 2 * math.pi * wind) > 0.5:
                return False, -1
            W = 6 + wind
            if W < 1 or W > 6:
                return False, -1
            return True, W
        if rd == rounds or left.size > cap:
            break
        nl = np.empty(left.size + nbad); nw = np.empty(left.size + nbad)
        m = 0
        for i in range(left.size):
            if lab[i] == 0:
                nl[m] = left[i]; nw[m] = width[i] / 2; m += 1
                nl[m] = left[i] + width[i] / 2; nw[m] = width[i] / 2; m += 1
            else:
                nl[m] = left[i]; nw[m] = width[i]; m += 1
        left = nl; width = nw
    return False, -1


# ---------------------------------------------------------------------------------------------------------------------
# a-chart (chart 0 of nb_cert): box variables (eta, l, w) with e = sin eta, b = cos eta, p = l sin w, q = l cos w;
# coefficient polynomials C(e, b, p, q) (gen_nb2a.py -> nbcoef_a2.py).  Same data layout as box_coeffs2.

from nbcoef_a2 import coef_val2 as coef_val2a, coef_abs2 as coef_abs2a
from nb_cert import dp_count, LIMIT

PAIRS4 = np.array([[0, 0], [0, 1], [0, 2], [0, 3], [1, 1], [1, 2], [1, 3], [2, 2], [2, 3], [3, 3]])
TRI4 = np.array([[a, b, c] for a in range(4) for b in range(a, 4) for c in range(b, 4)])
MULT4 = np.array([1.0 if (a == b and b == c) else (3.0 if (a == b or b == c) else 6.0) for a, b, c in TRI4])


@njit(cache=True)
def pair4(i, j):
    a = min(i, j); b = max(i, j)
    for m in range(10):
        if PAIRS4[m, 0] == a and PAIRS4[m, 1] == b:
            return m
    return -1


@njit(cache=True)
def box_coeffs2a(x0lo, x0hi, llo, lhi, wlo, whi):
    """a-chart data of nb_cert.box_coeffs(0, ...) plus hcen (Hessian in (eta, l, w) at the centre; order
    ee, el, ew, ll, lw, ww with e standing for eta) and rem (third-order remainder + evaluation errors).
    Along x_c + tau d: g1 = (d_eta cos, -d_eta sin, dl s + l dw c, dl c - l dw s), |g1| <= (r_eta, r_eta, a1, a1);
    |g2| <= (r_eta^2, r_eta^2, a2, a2); |g3| <= (r_eta^3, r_eta^3, a3, a3)."""
    cmid, crad, dcen, dwid, cmax, rr = box_coeffs(0, x0lo, x0hi, llo, lhi, wlo, whi)
    xc = 0.5 * (x0lo + x0hi); lc = 0.5 * (llo + lhi); wc = 0.5 * (wlo + whi)
    ec, bc = math.sin(xc), math.cos(xc)
    sw, cw = math.sin(wc), math.cos(wc)
    pc, qc = lc * sw, lc * cw
    re_, rl, rw = rr[0], rr[1], rr[2]
    lmax = lc + rl
    a1 = (rl + lmax * rw) * (1 + 1e-12) + 1e-15
    a2 = (2 * rl * rw + lmax * rw * rw) * (1 + 1e-12)
    a3 = (3 * rl * rw * rw + lmax * rw * rw * rw) * (1 + 1e-12)
    re1 = re_ * (1 + 1e-12) + 1e-15
    re2 = re1 * re1 * (1 + 1e-12)
    re3 = re2 * re1 * (1 + 1e-12)
    val = np.zeros((15, 2, 2, 7)); ab = np.zeros((35, 2, 2, 7))
    coef_val2a(ec, bc, pc, qc, val)
    coef_abs2a(abs(ec) + re1, abs(bc) + re1, abs(pc) + a1, abs(qc) + a1, ab)
    A = np.array([re1, re1, a1, a1])
    G2 = np.array([re2, re2, a2, a2])
    G3 = np.array([re3, re3, a3, a3])
    hcen = np.empty((6, 2, 2, 7)); rem = np.empty((2, 2, 7))
    for f in range(2):
        for pt in range(2):
            for k in range(7):
                Ce = val[1, f, pt, k]; Cb = val[2, f, pt, k]; Cp = val[3, f, pt, k]; Cq = val[4, f, pt, k]
                Cee = val[5, f, pt, k]; Ceb = val[6, f, pt, k]; Cep = val[7, f, pt, k]; Ceq = val[8, f, pt, k]
                Cbb = val[9, f, pt, k]; Cbp = val[10, f, pt, k]; Cbq = val[11, f, pt, k]
                Cpp = val[12, f, pt, k]; Cpq = val[13, f, pt, k]; Cqq = val[14, f, pt, k]
                dEp = Cep * bc - Cbp * ec
                dEq = Ceq * bc - Cbq * ec
                hcen[0, f, pt, k] = Cee * bc * bc - 2 * Ceb * ec * bc + Cbb * ec * ec - Ce * ec - Cb * bc
                hcen[1, f, pt, k] = dEp * sw + dEq * cw
                hcen[2, f, pt, k] = lc * (dEp * cw - dEq * sw)
                hcen[3, f, pt, k] = Cpp * sw * sw + 2 * Cpq * sw * cw + Cqq * cw * cw
                hcen[4, f, pt, k] = lc * (Cpp * sw * cw + Cpq * (cw * cw - sw * sw) - Cqq * sw * cw) + (Cp * cw - Cq * sw)
                hcen[5, f, pt, k] = lc * lc * (Cpp * cw * cw - 2 * Cpq * sw * cw + Cqq * sw * sw) - lc * (Cp * sw + Cq * cw)
                t3 = 0.0
                for m in range(20):
                    t3 += MULT4[m] * ab[15 + m, f, pt, k] * A[TRI4[m, 0]] * A[TRI4[m, 1]] * A[TRI4[m, 2]]
                t2 = 0.0
                for i in range(4):
                    for j in range(4):
                        t2 += ab[5 + pair4(i, j), f, pt, k] * G2[i] * A[j]
                t1 = 0.0
                for i in range(4):
                    t1 += ab[1 + i, f, pt, k] * G3[i]
                r3 = (t3 + 3.0 * t2 + t1) / 6.0
                s1 = ab[1, f, pt, k] + ab[2, f, pt, k] + ab[3, f, pt, k] + ab[4, f, pt, k]
                s2 = 0.0
                for m in range(10):
                    s2 += ab[5 + m, f, pt, k]
                e1 = 1e-13 * (s1 * (1.0 + lmax)) * (re_ + rl + rw)
                e2 = 1e-13 * ((s2 + s1) * (1.0 + lmax) * (1.0 + lmax)) * (re_ + rl + rw) * (re_ + rl + rw)
                rem[f, pt, k] = (r3 + e1 + e2) * (1 + 1e-12)
    return cmid, crad, dcen, dwid, cmax, rr, hcen, rem


@njit(cache=True)
def box_data2(chart, x0lo, x0hi, llo, lhi, wlo, whi):
    if chart == 0:
        return box_coeffs2a(x0lo, x0hi, llo, lhi, wlo, whi)
    return box_coeffs2(x0lo, x0hi, llo, lhi, wlo, whi)


@njit(cache=True)
def certify_box2(chart, x0lo, x0hi, llo, lhi, wlo, whi, t0, t1, J, M0, rounds):
    """nb_cert.certify_box (chart 0 or 1) with the second-order enclosure arc_signs2."""
    cmid, crad, dcen, dwid, cmax, rr, hcen, rem = box_data2(chart, x0lo, x0hi, llo, lhi, wlo, whi)
    left = np.arange(M0) * (2 * math.pi / M0)
    width = np.full(M0, 2 * math.pi / M0)
    total = 999
    prev = 10 ** 9
    stall = 0
    for rd in range(rounds + 1):
        sg, hopeless = arc_signs2(cmid, crad, dcen, dwid, cmax, rr, hcen, rem, t0, t1, left + width / 2,
                                  width / 2 + ARC_PAD, J)
        total = dp_count(sg)
        if total <= LIMIT:
            return True, total
        if rd == rounds or hopeless or left.size > 4096:
            break
        # stop refining once the count has not improved for two rounds: the box must be split instead
        if total >= prev:
            stall += 1
            if stall >= 2:
                break
        else:
            stall = 0
        prev = total
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
    return False, total


@njit(cache=True)
def certify_box_deep2(chart, x0lo, x0hi, llo, lhi, wlo, whi, t0, t1, J, M0, rounds, arc_cap):
    """nb_cert.certify_box_deep (chart 0 or 1) with the second-order enclosure arc_signs2."""
    cmid, crad, dcen, dwid, cmax, rr, hcen, rem = box_data2(chart, x0lo, x0hi, llo, lhi, wlo, whi)
    left = np.arange(M0) * (2 * math.pi / M0)
    width = np.full(M0, 2 * math.pi / M0)
    total = 999
    for rd in range(rounds + 1):
        sg, hopeless = arc_signs2(cmid, crad, dcen, dwid, cmax, rr, hcen, rem, t0, t1, left + width / 2,
                                  width / 2 + ARC_PAD, J)
        total = dp_count(sg)
        if total <= LIMIT:
            return True, total
        if rd == rounds or hopeless or left.size > arc_cap:
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
    return False, total

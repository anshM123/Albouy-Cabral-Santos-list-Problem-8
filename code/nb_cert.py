"""Numba certificate engine (both charts).  Same mathematics as certify5_ce, faster.

Box variables X = (x0, l, w, t): chart 'a': x0 = eta (e = sin eta, b = cos eta), eccentric anomaly;
                                 chart 'f': x0 = e, true anomaly, pericentre distance 1.
p = l sin w, q = l cos w, t in [0,1] with f_t = (1-t) Phi - t Psi.
Certified claim per box: for every parameter in the box, f_t has at most LIMIT real zeros on the circle
(with multiplicity), via Rolle runs over arcs.

Enclosures (rigorous up to the stated float padding):
  centre coefficients c(x_c) with radius 1e-13 * Pabs(|x_c|)            [rounding of inputs and evaluation]
  parameter variation, pointwise in u:  sum_v r_v ( |G_v^c(u)| + sum_k w_{v,k} |T_k(u)| )
      G_v^c = centre derivative (chain rule), w_{v,k} = bound on the variation of d c_k/d v over the box
      from the absolute-polynomial inequality |P(x) - P(x0)| <= Pabs(|x0| + r) - Pabs(|x0|);
  arc remainder  h |f^{(j+1)}| + h^2/2 sup |f^{(j+2)}|;  float padding 1e-12 x (sum of absolute terms).
"""
import math
import numpy as np
from numba import njit
from nbcoef_a import coef_val as cval_a, coef_abs as cabs_a
from nbcoef_f import coef_val as cval_f, coef_abs as cabs_f
from nbcoef_fe import coef_val as cval_fe, coef_abs as cabs_fe
from nbcoef_fm5 import coef_val as cval_fm5, coef_abs as cabs_fm5
from nbcoef_fm8 import coef_val as cval_fm8, coef_abs as cabs_fm8

LIMIT = 11          # parity lemma (PROOF.md section 2): an even count <= 11 is <= 10
ARC_PAD = 1e-12     # absolute inflation of every arc half-width: consecutive arcs overlap, the union covers [-pad, 2 pi + pad]


@njit(cache=True)
def box_coeffs(chart, x0lo, x0hi, llo, lhi, wlo, whi):
    """returns cmid[f,part,k], crad[f,part,k], dcen[v,f,part,k], dwid[v,f,part,k], cmax[f,part,k], r[3]"""
    xc = 0.5 * (x0lo + x0hi); lc = 0.5 * (llo + lhi); wc = 0.5 * (wlo + whi)
    r0 = 0.5 * (x0hi - x0lo); rl = 0.5 * (lhi - llo); rw = 0.5 * (whi - wlo)
    sw, cw = math.sin(wc), math.cos(wc)
    pc, qc = lc * sw, lc * cw
    lmax = lc + rl
    dp = rl + lmax * rw + 1e-15
    dq = rl + lmax * rw + 1e-15
    if chart == 0:
        ec, bc = math.sin(xc), math.cos(xc)
        de = r0 + 1e-15; db = r0 + 1e-15
        nv = 4
        val = np.zeros((5, 2, 2, 7)); a0 = np.zeros((5, 2, 2, 7)); a1 = np.zeros((5, 2, 2, 7))
        cval_a(ec, bc, pc, qc, val)
        cabs_a(abs(ec), abs(bc), abs(pc), abs(qc), a0)
        cabs_a(abs(ec) + de, abs(bc) + db, abs(pc) + dp, abs(qc) + dq, a1)
    elif chart == 1:
        ec, bc = xc, 0.0
        de = r0 + 1e-15; db = 0.0
        nv = 3
        val = np.zeros((4, 2, 2, 7)); a0 = np.zeros((4, 2, 2, 7)); a1 = np.zeros((4, 2, 2, 7))
        cval_f(ec, pc, qc, val)
        cabs_f(abs(ec), abs(pc), abs(qc), a0)
        cabs_f(abs(ec) + de, abs(pc) + dp, abs(qc) + dq, a1)
    elif chart == 3 or chart == 4:
        # Moebius-preconditioned f-charts (gen_nb fm5: c = 1/2, fm8: c = 4/5); variables (e, p, q) as in chart 1
        ec, bc = xc, 0.0
        de = r0 + 1e-15; db = 0.0
        nv = 3
        val = np.zeros((4, 2, 2, 7)); a0 = np.zeros((4, 2, 2, 7)); a1 = np.zeros((4, 2, 2, 7))
        if chart == 3:
            cval_fm5(ec, pc, qc, val); cabs_fm5(abs(ec), abs(pc), abs(qc), a0)
            cabs_fm5(abs(ec) + de, abs(pc) + dp, abs(qc) + dq, a1)
        else:
            cval_fm8(ec, pc, qc, val); cabs_fm8(abs(ec), abs(pc), abs(qc), a0)
            cabs_fm8(abs(ec) + de, abs(pc) + dp, abs(qc) + dq, a1)
    else:
        # chart 2 ('fe'): f-chart polynomial in eps = 1 - e; box variable x0 = e, so d/de = -d/deps
        ec, bc = 1.0 - xc, 0.0
        de = r0 + 1e-15 + 1e-16 * (abs(xc) + 1.0); db = 0.0
        nv = 3
        val = np.zeros((4, 2, 2, 7)); a0 = np.zeros((4, 2, 2, 7)); a1 = np.zeros((4, 2, 2, 7))
        cval_fe(ec, pc, qc, val)
        cabs_fe(abs(ec), abs(pc), abs(qc), a0)
        cabs_fe(abs(ec) + de, abs(pc) + dp, abs(qc) + dq, a1)
    cmid = np.empty((2, 2, 7)); crad = np.empty((2, 2, 7))
    dcen = np.empty((3, 2, 2, 7)); dwid = np.empty((3, 2, 2, 7)); cmax = np.empty((2, 2, 7))
    emax = abs(ec) + de; bmax = abs(bc) + db
    for f in range(2):
        for pt in range(2):
            for k in range(7):
                cmid[f, pt, k] = val[0, f, pt, k]
                crad[f, pt, k] = 1e-13 * a0[0, f, pt, k] + 1e-300
                cmax[f, pt, k] = a1[0, f, pt, k] * (1 + 1e-12)          # |c| over the box
                if chart == 0:
                    De, Db, Dp, Dq = val[1, f, pt, k], val[2, f, pt, k], val[3, f, pt, k], val[4, f, pt, k]
                    vDe = a1[1, f, pt, k] - a0[1, f, pt, k]; vDb = a1[2, f, pt, k] - a0[2, f, pt, k]
                    vDp = a1[3, f, pt, k] - a0[3, f, pt, k]; vDq = a1[4, f, pt, k] - a0[4, f, pt, k]
                    # d/d eta = De*b - Db*e
                    dcen[0, f, pt, k] = De * bc - Db * ec
                    dwid[0, f, pt, k] = (vDe * bmax + abs(De) * db + vDb * emax + abs(Db) * de) * (1 + 1e-12) + 1e-13 * (a1[1, f, pt, k] + a1[2, f, pt, k])
                else:
                    De, Dp, Dq = val[1, f, pt, k], val[2, f, pt, k], val[3, f, pt, k]
                    vDe = a1[1, f, pt, k] - a0[1, f, pt, k]
                    vDp = a1[2, f, pt, k] - a0[2, f, pt, k]; vDq = a1[3, f, pt, k] - a0[3, f, pt, k]
                    dcen[0, f, pt, k] = -De if chart == 2 else De
                    dwid[0, f, pt, k] = vDe * (1 + 1e-12) + 1e-13 * a1[1, f, pt, k]
                Dpa = abs(Dp); Dqa = abs(Dq)
                # d/dl = Dp sin w + Dq cos w ; d/dw = l (Dp cos w - Dq sin w)
                dcen[1, f, pt, k] = Dp * sw + Dq * cw
                dwid[1, f, pt, k] = (vDp + vDq + (Dpa + Dqa) * rw) * (1 + 1e-12) + 1e-13 * (Dpa + Dqa + vDp + vDq)
                dcen[2, f, pt, k] = lc * (Dp * cw - Dq * sw)
                dwid[2, f, pt, k] = (lmax * (vDp + vDq + (Dpa + Dqa) * rw) + rl * (Dpa + Dqa + vDp + vDq)) * (1 + 1e-12) \
                    + 1e-13 * lmax * (Dpa + Dqa + vDp + vDq)
    # mean-value radii padded for the float centre (|fl(mid) - mid| <= 1 ulp) and the float half-width
    rr = np.array([r0 * (1 + 1e-12) + 1e-15 * (abs(xc) + 1.0), rl * (1 + 1e-12) + 1e-15 * (abs(lc) + 1.0),
                   rw * (1 + 1e-12) + 1e-15 * (abs(wc) + 1.0)])
    return cmid, crad, dcen, dwid, cmax, rr


@njit(cache=True)
def arc_signs(cmid, crad, dcen, dwid, cmax, rr, t0, t1, um, hh, J):
    M = um.size
    pos = np.ones((J + 1, M), dtype=np.bool_)
    neg = np.ones((J + 1, M), dtype=np.bool_)
    hopeless = np.ones(M, dtype=np.bool_)
    for ti in range(2):
        t = t0 if ti == 0 else t1
        c = (1 - t) * cmid[0] - t * cmid[1]
        cr = (1 - t) * crad[0] + t * crad[1] + 1e-15 * (np.abs(cmid[0]) + np.abs(cmid[1]))
        dc = (1 - t) * dcen[:, 0] - t * dcen[:, 1]
        dw = (1 - t) * dwid[:, 0] + t * dwid[:, 1] + 1e-15 * (np.abs(dcen[:, 0]) + np.abs(dcen[:, 1]))
        cm = (1 - t) * cmax[0] + t * cmax[1]
        for i in range(M):
            u = um[i]; h = hh[i]
            vc = np.empty(J + 2); rad = np.empty(J + 2)
            for j in range(J + 2):
                s = 0.0; mag = 0.0; rd = 0.0
                g0 = 0.0; g1 = 0.0; g2 = 0.0; w0 = 0.0; w1 = 0.0; w2 = 0.0
                aab = 0.0      # sum over k of k^j x (all absolute coefficient magnitudes entering s, rd, prad)
                for k in range(7):
                    if j > 0 and k == 0:
                        continue
                    kj = float(k) ** j
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
                    g0 += dc[0, 0, k] * T + dc[0, 1, k] * U
                    g1 += dc[1, 0, k] * T + dc[1, 1, k] * U
                    g2 += dc[2, 0, k] * T + dc[2, 1, k] * U
                    w0 += dw[0, 0, k] * aT + dw[0, 1, k] * aU
                    w1 += dw[1, 0, k] * aT + dw[1, 1, k] * aU
                    w2 += dw[2, 0, k] * aT + dw[2, 1, k] * aU
                prad = rr[0] * (abs(g0) + w0) + rr[1] * (abs(g1) + w1) + rr[2] * (abs(g2) + w2)
                tot = rd + prad
                # float padding: summation (relative) + trig-argument rounding (|arg| < 50: abs. error < 2e-14 k^j)
                tot = tot + 1e-12 * (mag + tot) + 2e-14 * aab
                vc[j] = s; rad[j] = tot
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
def dp_count(sg):
    J1, M = sg.shape
    st = -1
    for i in range(M):
        if sg[0, i] != 0:
            st = i
            break
    if st < 0:
        return 999
    dp = np.empty(M - 1)
    rs = np.zeros(J1, dtype=np.int64)
    for i in range(M - 1):
        a = (st + 1 + i) % M
        best = 999.0
        for j in range(J1):
            v = sg[j, a]
            if v == 0:
                continue
            if i == 0:
                rs[j] = 0
            else:
                ap = (st + i) % M
                if sg[j, ap] != v:
                    rs[j] = i
            prev = dp[rs[j] - 1] if rs[j] > 0 else 0.0
            if prev + j < best:
                best = prev + j
        dp[i] = best
    return int(dp[M - 2]) if M > 1 else 0


@njit(cache=True)
def certify_box(chart, x0lo, x0hi, llo, lhi, wlo, whi, t0, t1, J, M0, rounds):
    cmid, crad, dcen, dwid, cmax, rr = box_coeffs(chart, x0lo, x0hi, llo, lhi, wlo, whi)
    left = np.arange(M0) * (2 * math.pi / M0)
    width = np.full(M0, 2 * math.pi / M0)
    total = 999
    prev = 10 ** 9
    stall = 0
    for rd in range(rounds + 1):
        sg, hopeless = arc_signs(cmid, crad, dcen, dwid, cmax, rr, t0, t1, left + width / 2, width / 2 + ARC_PAD, J)
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
def certify_box_deep(chart, x0lo, x0hi, llo, lhi, wlo, whi, t0, t1, J, M0, rounds, arc_cap):
    """certify_box without the stall break and with a configurable arc cap (used to re-check failed leaves)."""
    cmid, crad, dcen, dwid, cmax, rr = box_coeffs(chart, x0lo, x0hi, llo, lhi, wlo, whi)
    left = np.arange(M0) * (2 * math.pi / M0)
    width = np.full(M0, 2 * math.pi / M0)
    total = 999
    for rd in range(rounds + 1):
        sg, hopeless = arc_signs(cmid, crad, dcen, dwid, cmax, rr, t0, t1, left + width / 2, width / 2 + ARC_PAD, J)
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

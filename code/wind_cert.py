"""Winding-number certificate: the chart polynomial has at least W >= 1 roots strictly inside the disk |z| < rho < 1.

Setting: on a parameter box the chart function is f_t(u) = sum_k a_k cos ku + b_k sin ku (k = 0..6), affine in t in
[t0, t1], given by the coefficient enclosures (cmid, crad, dcen, dwid, cmax, rr) of nb_cert.box_coeffs or of a
generated engine (same layout: [f, part, k], 3 mean-value variables).  Put c_0 = a_0, c_k = (a_k - i b_k)/2,
c_{-k} = conj(c_k), and Q(z) = sum_k c_k z^{k+6} (degree <= 12; z = e^{iu}).  Q is self-inversive (its zeros are
symmetric under z -> 1/conj(z), zeros at 0 matching a degree drop), and f has the same real zeros, with the same
multiplicities, as Q on |z| = 1.

On |z| = rho = e^L (L < 0):  Q(rho e^{i th}) = rho^6 e^{6 i th} G(th) with
  G(th) = sum_k cosh(kL)(a_k cos k th + b_k sin k th)  +  i sum_k sinh(kL)(a_k sin k th - b_k cos k th).
If every arc of a cover of the circle carries a certified half-plane label for G (Re G > 0, Re G < 0, Im G > 0 or
Im G < 0, uniformly over the box and t in [t0, t1] -- G is affine in t), then G has no zero on |z| = rho anywhere in the
box, the winding number of G is constant on the box, and the number of zeros of Q in |z| < rho is
W = 6 + wind(G).  Then the number of zeros of f on the unit circle, with multiplicity, is at most 12 - 2W.
wind(G) is computed at the box centre from G at the arc boundaries: consecutive boundaries lie in one labelled arc,
so each argument step is < pi in modulus and the principal values add up exactly.

Enclosures per component use the scheme of nb_cert.arc_signs (value at the arc centre, parameter radius, Taylor
remainder, float padding including trig-argument rounding), with the harmonic weights cosh(kL) / sinh(kL) applied to
all coefficient data (the weights are computed in floating point and padded by a relative 1e-12)."""
import math
import numpy as np
from numba import njit

ARC_PAD = 1e-12


@njit(cache=True)
def comp_bounds(c, cr, dc, dw, cm, rr, wgt, phase, u, h):
    """enclosure [lo, hi] over the arc [u - h, u + h] and the parameter box of
       S(th) = sum_k wgt_k (A_k cos(k th + phase) + B_k sin(k th + phase))   (A, B = coefficient data c[0], c[1])."""
    vals = np.empty(2); rads = np.empty(2)
    for j in range(2):
        s = 0.0; mag = 0.0; rd = 0.0; aab = 0.0
        g = np.zeros(3); w = np.zeros(3)
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
            for v in range(3):
                g[v] += dc[v, 0, k] * T + dc[v, 1, k] * U
                w[v] += dw[v, 0, k] * aT + dw[v, 1, k] * aU
                ab_k += rr[v] * (abs(dc[v, 0, k]) + abs(dc[v, 1, k]) + dw[v, 0, k] + dw[v, 1, k])
            aab += abs(kj) * ab_k
        prad = 0.0
        for v in range(3):
            prad += rr[v] * (abs(g[v]) + w[v])
        tot = rd + prad
        tot = tot + 1e-12 * (mag + tot) + 2e-14 * aab
        vals[j] = s; rads[j] = tot
    B2 = 0.0
    for k in range(1, 7):
        B2 += abs(wgt[k]) * float(k) ** 2 * (cm[0, k] + cm[1, k])
    B2 *= 1 + 1e-12
    tail = rads[0] + h * (abs(vals[1]) + rads[1]) + 0.5 * h * h * B2
    tail = tail + 1e-15 * (abs(vals[0]) + tail)
    return vals[0] - tail, vals[0] + tail


@njit(cache=True)
def labels(cmid, crad, dcen, dwid, cmax, rr, t0, t1, L, um, hh):
    """label per arc: 1 Re>0, 2 Im>0, 3 Re<0, 4 Im<0, 0 none (certified for all params and t in [t0, t1])."""
    M = um.size
    wR = np.empty(7); wI = np.empty(7)
    for k in range(7):
        wR[k] = math.cosh(k * L) * (1 + 1e-12)
        wI[k] = math.sinh(k * L) * (1 + 1e-12)
    lab = np.zeros(M, dtype=np.int64)
    okR_p = np.ones(M, dtype=np.bool_); okR_n = np.ones(M, dtype=np.bool_)
    okI_p = np.ones(M, dtype=np.bool_); okI_n = np.ones(M, dtype=np.bool_)
    for ti in range(2):
        t = t0 if ti == 0 else t1
        c = (1 - t) * cmid[0] - t * cmid[1]
        cr = (1 - t) * crad[0] + t * crad[1] + 1e-15 * (np.abs(cmid[0]) + np.abs(cmid[1]))
        dc = (1 - t) * dcen[:, 0] - t * dcen[:, 1]
        dw = (1 - t) * dwid[:, 0] + t * dwid[:, 1] + 1e-15 * (np.abs(dcen[:, 0]) + np.abs(dcen[:, 1]))
        cm = (1 - t) * cmax[0] + t * cmax[1]
        # widen the weights' rounding into the radii: |wgt_computed - wgt| <= 1e-12 |wgt| is absorbed by
        # scaling cr, dw by (1 + 1e-12) relative to the weighted magnitudes (mag term) inside comp_bounds
        for i in range(M):
            loR, hiR = comp_bounds(c, cr, dc, dw, cm, rr, wR, 0.0, um[i], hh[i])
            # Im G = sum sinh(kL)(a_k sin k th - b_k cos k th) = sum sinh(kL)(a_k cos(k th - pi/2) + b_k sin(k th - pi/2))
            loI, hiI = comp_bounds(c, cr, dc, dw, cm, rr, wI, -math.pi / 2, um[i], hh[i])
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
def G_centre(cmid, t, L, th):
    c = (1 - t) * cmid[0] - t * cmid[1]
    re = 0.0; im = 0.0
    for k in range(7):
        ch = math.cosh(k * L); sh = math.sinh(k * L)
        re += ch * (c[0, k] * math.cos(k * th) + c[1, k] * math.sin(k * th))
        im += sh * (c[0, k] * math.sin(k * th) - c[1, k] * math.cos(k * th))
    return re, im


@njit(cache=True)
def arg_in_branch(lb, re, im):
    """argument of re + i im measured in the branch of half-plane label lb (audit R4-1):
    1 Re>0: (-pi/2, pi/2);  2 Im>0: (0, pi);  3 Re<0: (pi/2, 3pi/2);  4 Im<0: (-pi, 0).  Returns 99 if the value is
    not strictly inside the half-plane."""
    if lb == 1:
        if not re > 0:
            return 99.0
        return math.atan(im / re)
    if lb == 2:
        if not im > 0:
            return 99.0
        return math.atan2(im, re)            # in (0, pi) since im > 0
    if lb == 3:
        if not re < 0:
            return 99.0
        return math.pi + math.atan(im / re)  # in (pi/2, 3pi/2)
    if lb == 4:
        if not im < 0:
            return 99.0
        return math.atan2(im, re)            # in (-pi, 0) since im < 0
    return 99.0


@njit(cache=True)
def wind_box(cmid, crad, dcen, dwid, cmax, rr, t0, t1, L, M0, rounds, cap):
    """returns (ok, W) with W = number of zeros of Q in |z| < e^L, certified uniformly over the box (ok = True).
    Requires -15 <= L < 0 (audit R4-2/R4-3)."""
    if not (L < 0.0 and L >= -15.0):
        return False, -1
    left = np.arange(M0) * (2 * math.pi / M0)
    width = np.full(M0, 2 * math.pi / M0)
    for rd in range(rounds + 1):
        lab = labels(cmid, crad, dcen, dwid, cmax, rr, t0, t1, L, left + width / 2, width / 2 + ARC_PAD)
        nbad = 0
        for i in range(left.size):
            if lab[i] == 0:
                nbad += 1
        if nbad == 0:
            # winding at the centre (t = t0).  Boundary i (= left[i]) and boundary i+1 both lie in the closed arc i,
            # so G there lies in the arc's open half-plane; measuring both arguments in that half-plane's branch gives
            # the exact continuous argument change along the arc, a number in (-pi, pi).
            M = left.size
            vre = np.empty(M); vim = np.empty(M)
            for i in range(M):
                vre[i], vim[i] = G_centre(cmid, t0, L, left[i])
            tot = 0.0
            for i in range(M):
                j = i + 1 if i + 1 < M else 0          # wrap-around reuses the value at left[0] (audit R4-2)
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


def choose_L(cmid, t):
    """log-radius for the test circle, from the roots of Q at the box centre: maximise min|G| / scale over candidate
    circles that enclose at least one root (a heuristic choice; any L < 0 gives a valid certificate)."""
    c = (1 - t) * cmid[0] - t * cmid[1]
    cc = np.zeros(13, dtype=complex)
    cc[6] = c[0, 0]
    for k in range(1, 7):
        cc[6 + k] = (c[0, k] - 1j * c[1, k]) / 2
        cc[6 - k] = (c[0, k] + 1j * c[1, k]) / 2
    nz = np.nonzero(np.abs(cc) > 1e-300)[0]
    if nz.size == 0:
        return None
    hi = nz.max()
    poly = cc[:hi + 1][::-1]
    if len(poly) < 2:
        return None
    try:
        r = np.roots(poly)
    except Exception:
        return None
    lg = np.sort(np.log(np.maximum(np.abs(r), 1e-300)))
    inside = lg[lg < -1e-3]
    if nz.min() > 0:
        inside = np.concatenate([[-np.inf] * int(nz.min()), inside])      # exact zeros at 0
    if inside.size == 0:
        return None
    fin = np.sort(inside[np.isfinite(inside)])
    cands = []
    for i in range(fin.size):
        nxt = fin[i + 1] if i + 1 < fin.size else 0.0
        cands.append(0.5 * (fin[i] + min(nxt, 0.0)))
    if np.isinf(inside).any():
        cands.append((fin[0] - 0.5) if fin.size else -1.0)
    th = np.linspace(0, 2 * np.pi, 256, endpoint=False)
    best = None
    for L in cands:
        if not np.isfinite(L) or L >= -1e-3:
            continue
        k = np.arange(7)
        ch = np.cosh(k * L); sh = np.sinh(k * L)
        re = (ch[None, :] * (c[0][None, :] * np.cos(np.outer(th, k)) + c[1][None, :] * np.sin(np.outer(th, k)))).sum(1)
        im = (sh[None, :] * (c[0][None, :] * np.sin(np.outer(th, k)) - c[1][None, :] * np.cos(np.outer(th, k)))).sum(1)
        scale = (np.abs(ch) * (np.abs(c[0]) + np.abs(c[1]))).sum() + 1e-300
        score = np.sqrt(re * re + im * im).min() / scale
        if best is None or score > best[0]:
            best = (score, float(L))
    if best is None:
        return None
    return max(best[1], -15.0)

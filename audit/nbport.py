"""Faithful numpy port of nb_cert.box_coeffs / arc_signs that RETURNS the enclosure pieces
(vc, rad, m1, Bj2, lo, hi) instead of only signs.  Used by a5 to test the enclosures; validated against the
numba signs (a5 part A).  Independent evaluator `ft_derivs` (closed forms + FFT, no generated code) is used
as ground truth."""
import math
import numpy as np


def box_coeffs_port(chart, x0lo, x0hi, llo, lhi, wlo, whi, cval, cabs):
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
        val = np.zeros((5, 2, 2, 7)); a0 = np.zeros((5, 2, 2, 7)); a1 = np.zeros((5, 2, 2, 7))
        cval(ec, bc, pc, qc, val)
        cabs(abs(ec), abs(bc), abs(pc), abs(qc), a0)
        cabs(abs(ec) + de, abs(bc) + db, abs(pc) + dp, abs(qc) + dq, a1)
    else:
        ec, bc = xc, 0.0
        de = r0 + 1e-15; db = 0.0
        val = np.zeros((4, 2, 2, 7)); a0 = np.zeros((4, 2, 2, 7)); a1 = np.zeros((4, 2, 2, 7))
        cval(ec, pc, qc, val)
        cabs(abs(ec), abs(pc), abs(qc), a0)
        cabs(abs(ec) + de, abs(pc) + dp, abs(qc) + dq, a1)
    cmid = val[0].copy()
    crad = 1e-13 * a0[0] + 1e-300
    cmax = a1[0] * (1 + 1e-12)
    emax = abs(ec) + de; bmax = abs(bc) + db
    dcen = np.empty((3, 2, 2, 7)); dwid = np.empty((3, 2, 2, 7))
    if chart == 0:
        De, Db, Dp, Dq = val[1], val[2], val[3], val[4]
        vDe = a1[1] - a0[1]; vDb = a1[2] - a0[2]; vDp = a1[3] - a0[3]; vDq = a1[4] - a0[4]
        dcen[0] = De * bc - Db * ec
        dwid[0] = (vDe * bmax + np.abs(De) * db + vDb * emax + np.abs(Db) * de) * (1 + 1e-12) + 1e-13 * (a1[1] + a1[2])
    else:
        De, Dp, Dq = val[1], val[2], val[3]
        vDe = a1[1] - a0[1]; vDp = a1[2] - a0[2]; vDq = a1[3] - a0[3]
        dcen[0] = De
        dwid[0] = vDe * (1 + 1e-12) + 1e-13 * a1[1]
    Dpa = np.abs(Dp); Dqa = np.abs(Dq)
    dcen[1] = Dp * sw + Dq * cw
    dwid[1] = (vDp + vDq + (Dpa + Dqa) * rw) * (1 + 1e-12) + 1e-13 * (Dpa + Dqa + vDp + vDq)
    dcen[2] = lc * (Dp * cw - Dq * sw)
    dwid[2] = (lmax * (vDp + vDq + (Dpa + Dqa) * rw) + rl * (Dpa + Dqa + vDp + vDq)) * (1 + 1e-12) \
        + 1e-13 * lmax * (Dpa + Dqa + vDp + vDq)
    rr = np.array([r0, rl, rw])
    return cmid, crad, dcen, dwid, cmax, rr


def arc_bounds_port(cmid, crad, dcen, dwid, cmax, rr, t, um, hh, J):
    """returns dict of arrays [j, i]: vc, rad (pointwise-in-u parameter enclosure at um), lo, hi (arc x box)"""
    c = (1 - t) * cmid[0] - t * cmid[1]
    cr = (1 - t) * crad[0] + t * crad[1] + 1e-15 * (np.abs(cmid[0]) + np.abs(cmid[1]))
    dc = (1 - t) * dcen[:, 0] - t * dcen[:, 1]
    dw = (1 - t) * dwid[:, 0] + t * dwid[:, 1] + 1e-15 * (np.abs(dcen[:, 0]) + np.abs(dcen[:, 1]))
    cm = (1 - t) * cmax[0] + t * cmax[1]
    M = um.size
    ks = np.arange(7, dtype=float)
    VC = np.zeros((J + 2, M)); RAD = np.zeros((J + 2, M))
    for j in range(J + 2):
        ang = np.outer(um, ks) + j * math.pi / 2
        T = np.cos(ang) * ks ** j; U = np.sin(ang) * ks ** j
        if j > 0:
            T[:, 0] = 0.0; U[:, 0] = 0.0
        s = T @ c[0] + U @ c[1]
        mag = np.abs(T) @ np.abs(c[0]) + np.abs(U) @ np.abs(c[1])
        rd = np.abs(T) @ cr[0] + np.abs(U) @ cr[1]
        prad = np.zeros(M)
        for v in range(3):
            g = T @ dc[v, 0] + U @ dc[v, 1]
            w = np.abs(T) @ dw[v, 0] + np.abs(U) @ dw[v, 1]
            prad += rr[v] * (np.abs(g) + w)
        tot = rd + prad
        tot = tot + 1e-12 * (mag + tot)
        VC[j] = s; RAD[j] = tot
    LO = np.zeros((J + 1, M)); HI = np.zeros((J + 1, M)); BJ2 = np.zeros(J + 1)
    for j in range(J + 1):
        Bj2 = sum(float(k) ** (j + 2) * (cm[0, k] + cm[1, k]) for k in range(1, 7)) * (1 + 1e-12)
        m1 = np.abs(VC[j + 1]) + RAD[j + 1]
        LO[j] = VC[j] - RAD[j] - hh * m1 - 0.5 * hh * hh * Bj2
        HI[j] = VC[j] + RAD[j] + hh * m1 + 0.5 * hh * hh * Bj2
        BJ2[j] = Bj2
    return dict(vc=VC, rad=RAD, lo=LO, hi=HI, bj2=BJ2)


def signs_port(chart, box, um, hh, J, cval, cabs):
    x0lo, x0hi, llo, lhi, wlo, whi, t0, t1 = box
    pieces = box_coeffs_port(chart, x0lo, x0hi, llo, lhi, wlo, whi, cval, cabs)
    pos = np.ones((J + 1, um.size), bool); neg = np.ones((J + 1, um.size), bool)
    for t in (t0, t1):
        bd = arc_bounds_port(*pieces, t, um, hh, J)
        pos &= bd['lo'] > 0
        neg &= bd['hi'] < 0
    return np.where(pos, 1, np.where(neg, -1, 0))


# ---------------- independent ground truth ----------------
def phipsi_coeffs(chart, x0, l, w, N=32):
    """real Fourier coefficients (a[0..6], s[0..6]) of Phi and Psi from the closed forms (no generated code)."""
    u = 2 * np.pi * np.arange(N) / N
    p, q = l * np.sin(w), l * np.cos(w)
    if chart == 0:
        e, b = np.sin(x0), np.cos(x0)
        r = 1 - e * np.cos(u)
        Z = p * (np.cos(u) - e) + b * q * np.sin(u)
        Zp = -p * np.sin(u) + b * q * np.cos(u)
        A = 2 * e * np.sin(u) * r
        rho2 = r * r - Z * Z
        d = A - 2 * Z * Zp
        Phi = A * A * rho2; Psi = d * d
    else:
        e = x0
        zeta = p * np.cos(u) + q * np.sin(u); zp = -p * np.sin(u) + q * np.cos(u)
        sig2 = 1 - zeta * zeta
        D = 1 + e * np.cos(u)
        B = -2 * zeta * zp + 2 * e * (np.sin(u) - q * zeta)
        Phi = 4 * (1 + e) ** 2 * e ** 2 * np.sin(u) ** 2 * sig2
        Psi = D * D * B * B
    out = []
    for F in (Phi, Psi):
        c = np.fft.rfft(F) / N
        a = np.zeros(7); s = np.zeros(7)
        a[0] = c[0].real
        a[1:] = 2 * c[1:7].real
        s[1:] = -2 * c[1:7].imag
        out.append((a, s))
    return out


def ft_derivs(chart, x0, l, w, t, u, J):
    """f_t^{(j)}(u) for j = 0..J at points u (array); returns array [j, len(u)] and the coefficient magnitude."""
    (aP, sP), (aS, sS) = phipsi_coeffs(chart, x0, l, w)
    a = (1 - t) * aP - t * aS; s = (1 - t) * sP - t * sS
    ks = np.arange(7, dtype=float)
    out = np.zeros((J + 1, np.size(u)))
    for j in range(J + 1):
        ang = np.outer(np.atleast_1d(u), ks) + j * np.pi / 2
        out[j] = (np.cos(ang) * ks ** j) @ a + (np.sin(ang) * ks ** j) @ s
        if j > 0:
            out[j] -= 0.0
    mag = np.sum(np.abs(a)) + np.sum(np.abs(s))
    return out, mag

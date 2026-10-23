"""e -> 1 corner certificate (CORNER_E1.md), R-normalised f-chart: circle radius 1, pericentre distance mu, e in
[1-eps0, 1].

  F_-+(f) = mu k sin f sigma  -+  D B,   k = 2e(1+e),  D = 1 + e cos f,  zeta = Pz cos f + Qz sin f,
  sigma = sqrt(1 - zeta^2),  B = (sigma^2)' + 2e (sin f - Qz zeta).
  Zeros of F_- (F_+) are the stationary points of g_- (g_+).  P = F_- F_+ = mu^2 k^2 sin^2 f sigma^2 - D^2 B^2 is
  the f-chart polynomial with t = 1/(1 + mu^2).

Per parameter box (e, l, w) and mu-interval [mu0, mu1]:
  (a) near:  #zeros(F_-) on |f - pi| <= phi0 is <= 3       (interval Rolle runs, linear dp)
  (b) near:  #zeros(F_+) on |f - pi| <= phi0 is <= 1       (Lemma N2: B' < 0 there; |phi_B| <= 2|PzQz| eps/beta0;
             k(cos f sigma + sin f sigma') < 0 on the small interval; eps (2 C^2 Bmax eps) < beta0 eps; see check_plus)
  (c) away:  #real zeros of P on |f - pi| >= phi0 is <= 6  (nb_cert f-chart enclosures, linear dp)
  Hence at most 10 stationary points.
All enclosures: outward-rounded interval arithmetic (np.nextafter) in (a), (b); nb_cert enclosures in (c).
"""
import math
import numpy as np
from numba import njit
from nb_cert import box_coeffs, arc_signs, ARC_PAD

INF = np.inf
PAD = 1e-12      # every sub-range is widened by PAD on both sides (ranges overlap; true pi vs fl(pi) absorbed)


@njit(cache=True)
def dn(x):
    return np.nextafter(x, -INF)


@njit(cache=True)
def up(x):
    return np.nextafter(x, INF)


@njit(cache=True)
def iadd(a0, a1, b0, b1):
    return dn(a0 + b0), up(a1 + b1)


@njit(cache=True)
def isub(a0, a1, b0, b1):
    return dn(a0 - b1), up(a1 - b0)


@njit(cache=True)
def imul(a0, a1, b0, b1):
    p1 = a0 * b0; p2 = a0 * b1; p3 = a1 * b0; p4 = a1 * b1
    return dn(min(min(p1, p2), min(p3, p4))), up(max(max(p1, p2), max(p3, p4)))


@njit(cache=True)
def iscale(c, a0, a1):
    if c >= 0:
        return dn(c * a0), up(c * a1)
    return dn(c * a1), up(c * a0)


@njit(cache=True)
def idiv(a0, a1, b0, b1):
    # requires 0 < b0
    q1 = a0 / b0; q2 = a0 / b1; q3 = a1 / b0; q4 = a1 / b1
    return dn(min(min(q1, q2), min(q3, q4))), up(max(max(q1, q2), max(q3, q4)))


@njit(cache=True)
def isqrt(a0, a1):
    return dn(math.sqrt(max(a0, 0.0))), up(math.sqrt(a1))


@njit(cache=True)
def trig_arc(f0, f1):
    """enclosures of cos f, sin f for f in [f0, f1] subset [pi/2, 3pi/2] (cos <= 0 there, sin decreasing)."""
    c0 = math.cos(f0); c1 = math.cos(f1)
    clo = -1.0 if (f0 <= math.pi <= f1) else min(c0, c1)
    chi = max(c0, c1)
    slo = math.sin(f1); shi = math.sin(f0)
    return dn(dn(clo)), up(up(chi)), dn(dn(slo)), up(up(shi))


@njit(cache=True)
def near_derivs(f0, f1, e0, e1, pz0, pz1, qz0, qz1):
    """Interval enclosures over f in [f0,f1] and the parameter box of
       T^(j)/k and (DB)^(j) for j = 0..4, where T = k sin f sigma.  Returns arrays (5,2) for Ts = (sin f sigma)^(j),
       DB^(j), plus sigma enclosure, B and B' enclosures."""
    c0, c1, s0, s1 = trig_arc(f0, f1)
    # zeta and derivatives: zeta = Pz c + Qz s ; zeta' = -Pz s + Qz c ; zeta'' = -zeta ; zeta''' = -zeta'
    z0 = iadd(*imul(pz0, pz1, c0, c1), *imul(qz0, qz1, s0, s1))
    z1 = iadd(*imul(-pz1, -pz0, s0, s1), *imul(qz0, qz1, c0, c1))
    z2 = (-z0[1], -z0[0]); z3 = (-z1[1], -z1[0]); z4 = z0
    zz = imul(z0[0], z0[1], z0[0], z0[1])
    if z0[0] <= 0.0 <= z0[1]:
        zz = (0.0, zz[1])
    S2_0 = isub(1.0, 1.0, zz[0], zz[1])                                  # sigma^2
    S2_1 = iscale(-2.0, *imul(z0[0], z0[1], z1[0], z1[1]))              # -2 z z'
    t_a = imul(z1[0], z1[1], z1[0], z1[1])
    if z1[0] <= 0.0 <= z1[1]:
        t_a = (0.0, t_a[1])
    S2_2 = iscale(-2.0, *iadd(*t_a, *imul(z0[0], z0[1], z2[0], z2[1])))  # -2 (z'^2 + z z'')
    S2_3 = iscale(-2.0, *iadd(*iscale(3.0, *imul(z1[0], z1[1], z2[0], z2[1])), *imul(z0[0], z0[1], z3[0], z3[1])))
    t_b = imul(z2[0], z2[1], z2[0], z2[1])
    if z2[0] <= 0.0 <= z2[1]:
        t_b = (0.0, t_b[1])
    S2_4 = iscale(-2.0, *iadd(*iadd(*iscale(3.0, *t_b), *iscale(4.0, *imul(z1[0], z1[1], z3[0], z3[1]))),
                               *imul(z0[0], z0[1], z4[0], z4[1])))
    sg = isqrt(*S2_0)
    two_s = iscale(2.0, *sg)
    sg1 = idiv(*S2_1, *two_s)
    sg2 = idiv(*isub(*S2_2, *iscale(2.0, *imul(sg1[0], sg1[1], sg1[0], sg1[1]))), *two_s)
    sg3 = idiv(*isub(*S2_3, *iscale(6.0, *imul(sg1[0], sg1[1], sg2[0], sg2[1]))), *two_s)
    sg4 = idiv(*isub(*isub(*S2_4, *iscale(8.0, *imul(sg1[0], sg1[1], sg3[0], sg3[1]))),
                     *iscale(6.0, *imul(sg2[0], sg2[1], sg2[0], sg2[1]))), *two_s)
    SG = np.empty((5, 2))
    SG[0] = sg; SG[1] = sg1; SG[2] = sg2; SG[3] = sg3; SG[4] = sg4
    # sin f derivatives: s, c, -s, -c, s
    SN = np.empty((5, 2))
    SN[0, 0], SN[0, 1] = s0, s1
    SN[1, 0], SN[1, 1] = c0, c1
    SN[2, 0], SN[2, 1] = -s1, -s0
    SN[3, 0], SN[3, 1] = -c1, -c0
    SN[4, 0], SN[4, 1] = s0, s1
    binom = np.array([[1, 0, 0, 0, 0], [1, 1, 0, 0, 0], [1, 2, 1, 0, 0], [1, 3, 3, 1, 0], [1, 4, 6, 4, 1]], dtype=np.float64)
    Ts = np.zeros((5, 2))
    for j in range(5):
        lo = 0.0; hi = 0.0
        for i in range(j + 1):
            pr = imul(SN[i, 0], SN[i, 1], SG[j - i, 0], SG[j - i, 1])
            pr = iscale(binom[j, i], pr[0], pr[1])
            lo, hi = iadd(lo, hi, pr[0], pr[1])
        Ts[j, 0] = lo; Ts[j, 1] = hi
    # D derivatives: 1 + e c, -e s, -e c, e s, e c
    D = np.empty((5, 2))
    D[0] = iadd(1.0, 1.0, *imul(e0, e1, c0, c1))
    D[1] = iscale(-1.0, *imul(e0, e1, s0, s1))
    D[2] = iscale(-1.0, *imul(e0, e1, c0, c1))
    D[3] = imul(e0, e1, s0, s1)
    D[4] = imul(e0, e1, c0, c1)
    # B derivatives: B^(j) = (sigma^2)^(j+1) + 2e ( sin^(j) f - Qz zeta^(j) )
    S2d = np.empty((5, 2))
    S2d[0] = S2_1; S2d[1] = S2_2; S2d[2] = S2_3; S2d[3] = S2_4
    Zd = np.empty((5, 2))
    Zd[0] = z0; Zd[1] = z1; Zd[2] = z2; Zd[3] = z3; Zd[4] = z4
    Bd = np.empty((4, 2))
    for j in range(4):
        inner = isub(SN[j, 0], SN[j, 1], *imul(qz0, qz1, Zd[j, 0], Zd[j, 1]))
        Bd[j] = iadd(S2d[j, 0], S2d[j, 1], *iscale(2.0, *imul(e0, e1, inner[0], inner[1])))
    DB = np.zeros((4, 2))
    for j in range(4):
        lo = 0.0; hi = 0.0
        for i in range(j + 1):
            pr = imul(D[i, 0], D[i, 1], Bd[j - i, 0], Bd[j - i, 1])
            pr = iscale(binom[j, i], pr[0], pr[1])
            lo, hi = iadd(lo, hi, pr[0], pr[1])
        DB[j, 0] = lo; DB[j, 1] = hi
    return Ts, DB, SG, Bd, S2_0


@njit(cache=True)
def kint(e0, e1):
    return dn(2 * e0 * (1 + e0)), up(2 * e1 * (1 + e1))


@njit(cache=True)
def dp_linear(sg):
    """minimal Rolle count over a linear sequence of arcs (no wrap-around)."""
    J1, M = sg.shape
    dp = np.full(M + 1, 999.0)
    dp[0] = 0.0
    rs = np.zeros(J1, dtype=np.int64)
    for i in range(M):
        best = 999.0
        for j in range(J1):
            v = sg[j, i]
            if v == 0:
                continue
            if i == 0 or sg[j, i - 1] != v:
                rs[j] = i
            cand = dp[rs[j]] + j
            if cand < best:
                best = cand
        dp[i + 1] = best
    return dp[M]


@njit(cache=True)
def near_branch(sgn, fa, fb, e0, e1, pz0, pz1, qz0, qz1, mu0, mu1, M0, rounds, target):
    """Rolle-run bound for #zeros of F = mu k sin f sigma + sgn * DB on [fa, fb] (factor-based interval enclosures).
    sgn = -1: F_- (nearest point); sgn = +1: F_+ (farthest point)."""
    k0, k1 = kint(e0, e1)
    edges = np.linspace(fa, fb, M0 + 1)
    lefts = edges[:-1].copy(); rights = edges[1:].copy()
    best = 999.0
    for rd in range(rounds + 1):
        M = lefts.size
        sg = np.zeros((4, M), dtype=np.int64)
        for i in range(M):
            Ts, DB, SG, Bd, S2 = near_derivs(lefts[i], rights[i], e0, e1, pz0, pz1, qz0, qz1)
            if not (S2[0] > 0):
                return 999.0
            for j in range(4):
                Tk = imul(k0, k1, Ts[j, 0], Ts[j, 1])
                pos = True; neg = True
                for mu in (mu0, mu1):
                    v = iadd(*iscale(mu, Tk[0], Tk[1]), *iscale(float(sgn), DB[j, 0], DB[j, 1]))
                    if not (v[0] > 0):
                        pos = False
                    if not (v[1] < 0):
                        neg = False
                sg[j, i] = 1 if pos else (-1 if neg else 0)
        best = dp_linear(sg)
        if best <= target:
            return best
        if rd == rounds:
            break
        nl = np.empty(2 * M); nr = np.empty(2 * M)
        for i in range(M):
            mid = 0.5 * (lefts[i] + rights[i])
            nl[2 * i] = lefts[i]; nr[2 * i] = mid; nl[2 * i + 1] = mid; nr[2 * i + 1] = rights[i]
        lefts = nl; rights = nr
    return best


@njit(cache=True)
def near_minus(e0, e1, pz0, pz1, qz0, qz1, mu0, mu1, phi0, M0, rounds):
    """Rolle-run bound for #zeros of F_- = mu k sin f sigma - DB on [pi - phi0, pi + phi0]."""
    k0, k1 = kint(e0, e1)
    edges = np.linspace(math.pi - phi0, math.pi + phi0, M0 + 1)
    lefts = edges[:-1].copy(); rights = edges[1:].copy()
    best = 999.0
    for rd in range(rounds + 1):
        M = lefts.size
        sg = np.zeros((4, M), dtype=np.int64)
        for i in range(M):
            Ts, DB, SG, Bd, S2 = near_derivs(lefts[i], rights[i], e0, e1, pz0, pz1, qz0, qz1)
            if not (S2[0] > 0):
                return 999.0
            for j in range(4):
                Tk = imul(k0, k1, Ts[j, 0], Ts[j, 1])
                pos = True; neg = True
                for mu in (mu0, mu1):
                    v = isub(*iscale(mu, Tk[0], Tk[1]), DB[j, 0], DB[j, 1])
                    if not (v[0] > 0):
                        pos = False
                    if not (v[1] < 0):
                        neg = False
                sg[j, i] = 1 if pos else (-1 if neg else 0)
        best = dp_linear(sg)
        if best <= 3:
            return best
        if rd == rounds:
            break
        nl = np.empty(2 * M); nr = np.empty(2 * M)
        for i in range(M):
            mid = 0.5 * (lefts[i] + rights[i])
            nl[2 * i] = lefts[i]; nr[2 * i] = mid; nl[2 * i + 1] = mid; nr[2 * i + 1] = rights[i]
        lefts = nl; rights = nr
    return best


@njit(cache=True)
def check_plus(e0, e1, pz0, pz1, qz0, qz1, phi0):
    """Lemma N2 on [pi - phi0, pi + phi0]: returns True if
       (1) B' <= -beta0 < 0 on the arc, (2) with C = 2 max|Pz Qz| / beta0 the interval |phi| <= C eps1 lies in the arc,
       (3) k (cos f sigma + sin f sigma') < 0 on |phi| <= C eps1, and (4) 2 C^2 Bp_abs eps1 < beta0,
       where Bp_abs bounds |B'| on the arc and eps1 = 1 - e0.  Then F_+ has at most one zero on the arc (CORNER_E1.md)."""
    Ts, DB, SG, Bd, S2 = near_derivs(math.pi - phi0, math.pi + phi0, e0, e1, pz0, pz1, qz0, qz1)
    if not (S2[0] > 0):
        return False
    if not (Bd[1, 1] < 0):
        return False
    beta0 = -Bd[1, 1]
    bpabs = max(abs(Bd[1, 0]), abs(Bd[1, 1]))
    pq = max(abs(pz0), abs(pz1)) * max(abs(qz0), abs(qz1))
    eps1 = up(1.0 - e0)
    C = up(up(2 * pq) / beta0)
    w = up(C * eps1)
    if not (w < phi0):
        return False
    Ts2, DB2, SG2, Bd2, S22 = near_derivs(math.pi - w - PAD, math.pi + w + PAD, e0, e1, pz0, pz1, qz0, qz1)
    if not (Ts2[1, 1] < 0):           # (sin f sigma)' < 0 on the small interval
        return False
    if not (up(up(2 * C * C) * up(bpabs * eps1)) < beta0):
        return False
    return True


@njit(cache=True)
def dp_away(sg, lefts, phi0):
    """arcs cover [0, 2pi) minus (pi - phi0, pi + phi0): run the linear dp on the sequence starting after pi + phi0."""
    return dp_linear(sg)


@njit(cache=True)
def away_count(e0, e1, l0, l1, w0, w1, t0, t1, phi0, M0, rounds, J):
    """Rolle-run bound for real zeros of the f-chart polynomial on f in [pi + phi0, 3 pi - phi0] (mod 2 pi)."""
    cmid, crad, dcen, dwid, cmax, rr = box_coeffs(1, e0, e1, l0, l1, w0, w1)
    a = math.pi + phi0 - PAD; L = 2 * math.pi - 2 * phi0 + 2 * PAD
    left = a + np.arange(M0) * (L / M0)
    width = np.full(M0, L / M0)
    best = 999.0
    for rd in range(rounds + 1):
        sg, hopeless = arc_signs(cmid, crad, dcen, dwid, cmax, rr, t0, t1, left + width / 2, width / 2 + ARC_PAD, J)
        best = dp_linear(sg)
        if best <= 6:
            return best
        if rd == rounds or hopeless:
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
    return best


@njit(cache=True)
def certify_corner_box(e0, e1, l0, l1, w0, w1, mu0, mu1, phi0, phi_in):
    """Total = #F_- zeros on |phi| <= phi0  +  1 (Lemma N2: F_+ on |phi| <= phi_in)
             + #F_+ zeros on phi_in <= |phi| <= phi0  +  #real zeros of P on |phi| >= phi0   <= 10 ?
    returns (ok, total, parts...)"""
    # libm sin/cos padded by two ulps before the outward-rounded products (sin increasing, cos decreasing on [0, pi/2])
    pz0 = dn(l0 * dn(dn(math.sin(w0)))); pz1 = up(l1 * up(up(math.sin(w1))))
    qz0 = dn(l0 * dn(dn(math.cos(w1)))); qz1 = up(l1 * up(up(math.cos(w0))))
    pz0 = max(pz0, 0.0); qz0 = max(qz0, 0.0)
    pl = check_plus(e0, e1, pz0, pz1, qz0, qz1, phi_in + PAD)
    if not pl:
        return False, 999.0, 999.0, 999.0, 999.0
    nm = near_branch(-1, math.pi - phi0 - PAD, math.pi + phi0 + PAD, e0, e1, pz0, pz1, qz0, qz1, mu0, mu1, 32, 6, 3)
    po = near_branch(1, math.pi + phi_in - PAD, math.pi + phi0 + PAD, e0, e1, pz0, pz1, qz0, qz1, mu0, mu1, 16, 6, 0) + \
        near_branch(1, math.pi - phi0 - PAD, math.pi - phi_in + PAD, e0, e1, pz0, pz1, qz0, qz1, mu0, mu1, 16, 6, 0)
    t0 = 1.0 / (1.0 + mu1 * mu1); t1 = 1.0 / (1.0 + mu0 * mu0)
    # parity lemma (PROOF.md section 2): the real-root count of P with multiplicity is even, so a total <= 11 suffices
    budget = 11.0 - nm - 1.0 - po
    if budget < 0:
        return False, nm + 1 + po, nm, po, 999.0
    aw = away_count(e0, e1, l0, l1, w0, w1, t0, t1, phi0, 64, 8, 6)
    tot = nm + 1.0 + po + aw
    return (tot <= 11.0), tot, nm, po, aw

"""tau-winding certificate (TW): at most 10 real zeros of P, uniformly over a KQ parameter box (C3_CORE.md section 8).

Q(th) = P(f)(1 + th^2)^6, th = tan((f - pi)/2), is a real polynomial of degree <= 12 in th (c3 chart KBQe; its
coefficients are polynomials in (eps, Pz, Qz, mu)).  Its reverse G(tp) = tp^12 Q(-1/tp), tp = tan(f/2) = -1/th, is chart
KBPe.  The real zeros of P on the circle counted with multiplicity are the real zeros of Q plus a zero of order
12 - deg Q at f = 0 (th = infinity); with the non-real zeros of Q they number exactly 12.  So if Q (or G) has a zero
strictly inside a closed disc that lies in the open upper half-plane, its conjugate is a second non-real zero and P has
at most 12 - 2 = 10 real zeros (and the same for G: a non-real zero tp of G gives the non-real zero -1/tp of Q).

Certificate.  F = Q or G.  Contour: the circle |z - c| = r with Im c > r > 0, covered by closed arcs [a_j - h_j, a_j + h_j]
of the angle (contiguous, each sharing its end points with its neighbours); z_j = c + r e^{i a_j} (float), and every
point z of arc j satisfies |z - z_j| <= D_j = r h_j (1 + 1e-12) + 1e-12 (|c| + r).  For every x of the box and every z
with |z - z_j| <= D_j,
    |F(z; x) - F(z_j; x_c)| <= dev_j + T_j,
    dev_j = sum_{m>=1} |F^(m)(z_j; x_c)|/m! D_j^m          (exact Taylor expansion of the polynomial F(.; x_c) at z_j)
    T_j  >= |F(z; x) - F(z; x_c)|                         (parameter variation on the disc)
where T_j is the engine's parameter-variation bound (c3/make_c3_2.py point_enclosure2: the smaller of the mean-value
form and the second-order form, plus the coefficient-evaluation pad crad and the 1e-12 summation pads), with every
modulus of a complex polynomial p in z replaced by its maximum over the disc, bounded by sum_m |p^(m)(z_j)|/m! D_j^m,
and every power |z|^k by (|z_j| + D_j)^k.  Floating-point errors of the centre data are covered by 1e-12 of the
corresponding sums of absolute terms (as in FLOAT_RIGOUR.md section 5).  If |F(z_j; x_c)| > R_j = (dev_j + T_j +
pads)(1 + 1e-12) on every arc, F has no zero on the contour for any x of the box; the number of zeros inside is then
constant on the (connected) box, equal to its value at x_c, the winding number of F(.; x_c) along the contour.  Each
disc of radius R_j < |F(z_j; x_c)| about F(z_j; x_c) keeps the argument within pi/2 of the centre value, so the winding
number is the sum of the principal arguments of F(z_{j+1}; x_c)/F(z_j; x_c), divided by 2 pi.  The choice of the
contour (from the zeros at the box centre) affects only success, never validity."""
import cmath
import importlib
import math
import sys

import numpy as np
from numba import njit

sys.path.insert(0, 'c3')
_ENG = {}

BINOM = np.array([[math.comb(k, m) if m <= k else 0 for m in range(13)] for k in range(13)], dtype=np.float64)


def engine(which):
    """which 0: KBQe (th around aphelion), 1: KBPe (tp around perihelion)"""
    if which not in _ENG:
        _ENG[which] = importlib.import_module('nbeng_c3_KBQe_2' if which == 0 else 'nbeng_c3_KBPe_2')
    return _ENG[which]


def dn(x):
    return np.nextafter(np.nextafter(x, -np.inf), -np.inf)


def up(x):
    return np.nextafter(np.nextafter(x, np.inf), np.inf)


def params(bx):
    """the (eps, Pz, Qz, mu) box of a KQ box (e, l, w, t), outward-rounded exactly as in kq_run.count_total"""
    e0, e1, l0, l1, w0, w1, t0, t1 = bx
    pz0 = dn(l0 * dn(math.sin(w0))); pz1 = up(l1 * up(math.sin(w1)))
    qz0 = dn(l0 * dn(math.cos(w1))); qz1 = up(l1 * up(math.cos(w0)))
    ep0 = max(dn(1.0 - e1), 0.0); ep1 = up(1.0 - e0)
    t1c = min(t1, 1.0)
    mu0 = max(dn(math.sqrt(max((1.0 - t1c) / t1c, 0.0))), 0.0); mu1 = up(math.sqrt((1.0 - t0) / t0))
    return ep0, ep1, pz0, pz1, qz0, qz1, mu0, mu1


@njit(cache=True)
def disc_max(a, z, D, binom):
    """(|p(z)|, sum_{m>=1} |p^(m)(z)|/m! D^m) for the polynomial p = sum a_k t^k (complex a, k <= 12)"""
    v, b0, dev = disc_full(a, z, D, binom)
    return v, dev


@njit(cache=True)
def disc_full(a, z, D, binom):
    """(|p(z)|, p(z), sum_{m>=1} |p^(m)(z)|/m! D^m): the Taylor coefficients b_m of p at z, p(z + d) = sum b_m d^m"""
    n = a.shape[0]
    zp = np.empty(n, dtype=np.complex128)
    zp[0] = 1.0
    for k in range(1, n):
        zp[k] = zp[k - 1] * z
    v0 = 0.0
    b0 = 0.0 + 0.0j
    dev = 0.0
    Dm = 1.0
    for m in range(n):
        b = 0.0 + 0.0j
        for k in range(m, n):
            b += a[k] * binom[k, m] * zp[k - m]
        if m == 0:
            v0 = abs(b); b0 = b
        else:
            dev += abs(b) * Dm
        Dm *= D
    return v0, b0, dev


@njit(cache=True)
def arc_radius(cmid, crad, dcen, dwid, cmax, rr, hcen, rem, qw, z, D, binom):
    """(|F(z; x_c)|, R): F(w; x) lies within R of F(z; x_c) for every |w - z| <= D and every x of the box"""
    n = cmid.shape[0]
    rho = abs(z) + D
    rp = np.empty(n)
    rp[0] = 1.0
    for k in range(1, n):
        rp[k] = rp[k - 1] * rho
    cm = cmid.astype(np.complex128)
    v0, b0, dev = disc_full(cm, z, D, binom)
    mag = 0.0; rd = 0.0; rs = 0.0
    for k in range(n):
        mag += abs(cmid[k]) * rp[k]
        rd += crad[k] * rp[k]
        rs += rem[k] * rp[k]
    first = 0.0; old = 0.0; pabs = rs
    for v in range(4):
        gv, gd = disc_max(dcen[v].astype(np.complex128), z, D, binom)
        wv = 0.0; pa = 0.0
        for k in range(n):
            wv += dwid[v, k] * rp[k]
            pa += abs(dcen[v, k]) * rp[k]
        first += rr[v] * (gv + gd)
        old += rr[v] * (gv + gd + wv)
        pabs += rr[v] * pa
    second = 0.0
    for m in range(10):
        hv, hd = disc_max(hcen[m].astype(np.complex128), z, D, binom)
        second += qw[m] * (hv + hd)
        pa = 0.0
        for k in range(n):
            pa += abs(hcen[m, k]) * rp[k]
        pabs += qw[m] * pa
    new = first + second + rs + 1e-12 * pabs
    T = rd + min(old, new)
    T = T + 1e-12 * (mag + T)
    R = (dev + T + 1e-12 * mag) * (1.0 + 1e-12)
    return v0, b0, R


@njit(cache=True)
def centre_value(cmid, z):
    s = 0.0 + 0.0j
    zp = 1.0 + 0.0j
    for k in range(cmid.shape[0]):
        s += cmid[k] * zp
        zp *= z
    return s


@njit(cache=True)
def contour_ok(cmid, crad, dcen, dwid, cmax, rr, hcen, rem, qw, c, r, M0, rounds, cap, binom):
    """(ok, W): F has no zero on |z - c| = r for any x of the box, and W zeros inside at x_c"""
    a = np.empty(M0); h = np.empty(M0)
    for j in range(M0):
        a[j] = 2.0 * np.pi * (j + 0.5) / M0
        h[j] = np.pi / M0
    for rd in range(rounds + 1):
        n = a.shape[0]
        good = np.empty(n, dtype=np.bool_)
        vals = np.empty(n, dtype=np.complex128)
        nbad = 0
        for j in range(n):
            zj = c + r * np.exp(1j * a[j])
            D = r * h[j] * (1.0 + 1e-12) + 1e-12 * (abs(c) + r)
            v0, b0, R = arc_radius(cmid, crad, dcen, dwid, cmax, rr, hcen, rem, qw, zj, D, binom)
            good[j] = v0 > R
            vals[j] = b0
            if not good[j]:
                nbad += 1
        if nbad == 0:
            # winding number at x_c: every arc's values lie in a disc about vals[j] that excludes 0, and consecutive
            # arcs share an end point, so each step arg(vals[j+1]/vals[j]) is the principal value
            tot = 0.0
            for j in range(n):
                q = vals[(j + 1) % n] / vals[j]
                tot += math.atan2(q.imag, q.real)
            return True, int(np.round(tot / (2.0 * np.pi)))
        if rd == rounds or n + nbad > cap:
            return False, -1
        a2 = np.empty(n + nbad); h2 = np.empty(n + nbad)
        m = 0
        for j in range(n):
            if good[j]:
                a2[m] = a[j]; h2[m] = h[j]; m += 1
            else:
                a2[m] = a[j] - h[j] / 2; h2[m] = h[j] / 2; m += 1
                a2[m] = a[j] + h[j] / 2; h2[m] = h[j] / 2; m += 1
        a = a2; h = h2
    return False, -1


def candidates(cmid):
    """contours (c, r) from the zeros of F(.; x_c): discs in the open upper half-plane about a zero with Im > 0, about
    the midpoint of two such zeros, or about the mean of all of them (a cluster of zeros may be enclosed together: only
    W >= 1 is needed), with the radius halfway between two consecutive zero distances and below Im c; largest margin
    (distance of the circle from the zeros and from the real axis) first"""
    co = cmid[::-1]
    nz = np.nonzero(np.abs(co) > 0)[0]
    if nz.size == 0:
        return []
    co = co[nz[0]:]
    if co.size < 2:
        return []
    rts = [complex(z) for z in np.roots(co)]
    upz = [z for z in rts if z.imag > 0]
    cen = list(upz)
    for i in range(len(upz)):
        for j in range(i + 1, len(upz)):
            cen.append(0.5 * (upz[i] + upz[j]))
    if len(upz) > 2:
        cen.append(sum(upz) / len(upz))
    out = []
    for c in cen:
        ds = sorted(abs(c - w) for w in rts) + [10.0]
        for k in range(1, len(ds)):
            lo = ds[k - 1]; hi = min(ds[k], c.imag)
            if not hi > lo:
                break
            r = 0.5 * (lo + hi)
            out.append((0.5 * (hi - lo) / max(1.0, abs(c) + r), c, r))
    out.sort(key=lambda u: -u[0])
    res = []
    for _, c, r in out:
        if all(abs(c - c2) > 1e-3 * r or abs(r - r2) > 1e-3 * r for c2, r2 in res):
            res.append((c, r))
    return res


DISC_W = 6     # disc-winding outer count: W >= DISC_W zeros of Q in |th| < T_X -> at most 12 - W <= 6 real zeros of P
               # with |th| >= T_X (f = 0 included), the outer count of the face construction (C3_CORE.md section 8)


def certify_disc(bx, tx, M0=64, rounds=8, cap=16384):
    """(ok, W): Q = F of chart KBQe has no zero on the circle |th| = tx for any x of the box, and W >= DISC_W zeros in
    the disc |th| < tx at x_c (hence for every x of the box): at most 12 - W real zeros of P on |th| >= tx"""
    data = engine(0).box_coeffs2(*params(bx))
    ok, W = contour_ok(*data, 0j, tx, M0, rounds, cap, BINOM)
    return (ok and W >= DISC_W), W


def check_disc(bx, tx, M0=64, rounds=10, cap=65536):
    """re-certification of a disc-winding leaf (True if it certifies at most 6 real zeros on |th| >= tx)"""
    return certify_disc(bx, tx, M0, rounds, cap)[0]


def certify_tw(bx, M0=64, rounds=6, cap=4096, tries=3):
    """(ok, which, c_re, c_im, r): a contour certifying <= 10 real zeros on the whole box, or ok = False"""
    p = params(bx)
    for which in (0, 1):
        data = engine(which).box_coeffs2(*p)
        cmid = data[0]
        for c, r in candidates(cmid)[:tries]:
            if abs(c) > 1.5 or not (c.imag > r > 0):
                continue
            ok, W = contour_ok(*data, c, r, M0, rounds, cap, BINOM)
            if ok and W >= 1:
                return True, which, c.real, c.imag, r
    return False, -1, 0.0, 0.0, 0.0


TW_CHARTS = ('Ca', 'FE', 'FM')   # c3 charts with a second-order engine whose root variable is a positive multiple of th:
                                  # Ca: G(y) = Q(nu y)/nu^10, FE/FM: G(y) = Q(s y)/s^6 (nu, s > 0 at physical points)


def certify_tw_chart(eng, bx, M0=64, rounds=6, cap=4096, tries=3):
    """TW on a c3 chart box (box variables of the chart, contour in the chart's root variable): (ok, c_re, c_im, r)"""
    data = eng.box_coeffs2(*bx)
    for c, r in candidates(data[0])[:tries]:
        if not (c.imag > r > 0):
            continue
        ok, W = contour_ok(*data, c, r, M0, rounds, cap, BINOM)
        if ok and W >= 1:
            return True, c.real, c.imag, r
    return False, 0.0, 0.0, 0.0


def check_tw_chart(chart, bx, cre, cim, r, M0=64, rounds=10, cap=65536):
    """re-certification of a stored contour of a c3 chart leaf"""
    if chart not in TW_CHARTS or not (cim > r > 0):
        return False
    eng = importlib.import_module('nbeng_c3_%s_2' % chart)
    ok, W = contour_ok(*eng.box_coeffs2(*bx), complex(cre, cim), r, M0, rounds, cap, BINOM)
    return ok and W >= 1


def check_tw(bx, which, cre, cim, r, M0=64, rounds=10, cap=65536):
    """re-certification of a stored contour (True if it certifies <= 10 real zeros on the box)"""
    data = engine(which).box_coeffs2(*params(bx))
    if not (cim > r > 0):
        return False
    ok, W = contour_ok(*data, complex(cre, cim), r, M0, rounds, cap, BINOM)
    return ok and W >= 1

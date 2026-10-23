"""Write nbeng_<chart>.py for the blow-up charts of gen_bu.py.

Affine charts (B1, B2, L2E, L2P, L2Ap, L2Am): box (x0, x1, x2) handled by the mean-value enclosure (chain rule through
the atoms), x3 in [x3lo, x3hi] enters affinely (f = (1 - x3) F_A - x3 F_B): same certificate as nb_cert / make_engines.
4-variable charts (B3p, B3m, L2Tp, L2Tm): all four box variables by the mean-value enclosure, one function F.
Certificate logic (arcs, Rolle runs, refinement, padding) is identical to nb_cert.py; ARC_PAD overlaps the arcs.
usage: python make_engines_bu.py CHART [CHART ...]"""
import sys

S0 = 0.01   # level-2 scale of e (gen_bu.SC2[0])

ATOMS = {
    'B1': '''
@njit(cache=True)
def atoms_box(x0lo, x0hi, x1lo, x1hi, x2lo, x2hi, x3lo, x3hi):
    # atoms e, ph, qh, b, ib ; box (e, ph, qh)
    lo = np.empty(5); hi = np.empty(5); ac = np.empty(5)
    ec = 0.5 * (x0lo + x0hi)
    lo[0] = x0lo; hi[0] = x0hi; ac[0] = ec
    lo[1] = x1lo; hi[1] = x1hi; ac[1] = 0.5 * (x1lo + x1hi)
    lo[2] = x2lo; hi[2] = x2hi; ac[2] = 0.5 * (x2lo + x2hi)
    lo[3] = math.sqrt(1 - x0hi * x0hi); hi[3] = math.sqrt(1 - x0lo * x0lo); ac[3] = math.sqrt(1 - ec * ec)
    lo[4] = 1 / hi[3]; hi[4] = 1 / lo[3]; ac[4] = 1 / ac[3]
    return ac, lo, hi
''',
    'B2': '''
@njit(cache=True)
def atoms_box(x0lo, x0hi, x1lo, x1hi, x2lo, x2hi, x3lo, x3hi):
    # atoms th, eh, b, ib, sw, cw ; box (th, eh, w); e = th eh >= 0
    lo = np.empty(6); hi = np.empty(6); ac = np.empty(6)
    tc = 0.5 * (x0lo + x0hi); hc = 0.5 * (x1lo + x1hi); wc = 0.5 * (x2lo + x2hi)
    lo[0] = x0lo; hi[0] = x0hi; ac[0] = tc
    lo[1] = x1lo; hi[1] = x1hi; ac[1] = hc
    emax = x0hi * x1hi; emin = x0lo * x1lo; ecen = tc * hc
    lo[2] = math.sqrt(1 - emax * emax); hi[2] = math.sqrt(1 - emin * emin); ac[2] = math.sqrt(1 - ecen * ecen)
    lo[3] = 1 / hi[2]; hi[3] = 1 / lo[2]; ac[3] = 1 / ac[2]
    lo[4] = math.sin(x2lo); hi[4] = math.sin(x2hi); ac[4] = math.sin(wc)
    lo[5] = math.cos(x2hi); hi[5] = math.cos(x2lo); ac[5] = math.cos(wc)
    return ac, lo, hi
''',
    'B3': '''
@njit(cache=True)
def atoms_box(x0lo, x0hi, x1lo, x1hi, x2lo, x2hi, x3lo, x3hi):
    # atoms eta, eh, thh, b, ib, sw, cw ; box (eta, eh, thh, w); e = eta eh >= 0
    lo = np.empty(7); hi = np.empty(7); ac = np.empty(7)
    nc = 0.5 * (x0lo + x0hi); hc = 0.5 * (x1lo + x1hi); tc = 0.5 * (x2lo + x2hi); wc = 0.5 * (x3lo + x3hi)
    lo[0] = x0lo; hi[0] = x0hi; ac[0] = nc
    lo[1] = x1lo; hi[1] = x1hi; ac[1] = hc
    lo[2] = x2lo; hi[2] = x2hi; ac[2] = tc
    emax = x0hi * x1hi; emin = x0lo * x1lo; ecen = nc * hc
    lo[3] = math.sqrt(1 - emax * emax); hi[3] = math.sqrt(1 - emin * emin); ac[3] = math.sqrt(1 - ecen * ecen)
    lo[4] = 1 / hi[3]; hi[4] = 1 / lo[3]; ac[4] = 1 / ac[3]
    lo[5] = math.sin(x3lo); hi[5] = math.sin(x3hi); ac[5] = math.sin(wc)
    lo[6] = math.cos(x3hi); hi[6] = math.cos(x3lo); ac[6] = math.cos(wc)
    return ac, lo, hi
''',
    'L2E': '''
@njit(cache=True)
def atoms_box(x0lo, x0hi, x1lo, x1hi, x2lo, x2hi, x3lo, x3hi):
    # atoms nu, P, A, b, ib ; box (nu, P, A); e = S0 nu
    lo = np.empty(5); hi = np.empty(5); ac = np.empty(5)
    nc = 0.5 * (x0lo + x0hi)
    lo[0] = x0lo; hi[0] = x0hi; ac[0] = nc
    lo[1] = x1lo; hi[1] = x1hi; ac[1] = 0.5 * (x1lo + x1hi)
    lo[2] = x2lo; hi[2] = x2hi; ac[2] = 0.5 * (x2lo + x2hi)
    emax = S0 * x0hi; emin = S0 * x0lo; ecen = S0 * nc
    lo[3] = math.sqrt(1 - emax * emax); hi[3] = math.sqrt(1 - emin * emin); ac[3] = math.sqrt(1 - ecen * ecen)
    lo[4] = 1 / hi[3]; hi[4] = 1 / lo[3]; ac[4] = 1 / ac[3]
    return ac, lo, hi
''',
    'L2X': '''
@njit(cache=True)
def atoms_box(x0lo, x0hi, x1lo, x1hi, x2lo, x2hi, x3lo, x3hi):
    # atoms nu, E, (A or P), b, ib ; box (nu, E, A|P); e = S0 E nu >= 0
    lo = np.empty(5); hi = np.empty(5); ac = np.empty(5)
    nc = 0.5 * (x0lo + x0hi); Ec = 0.5 * (x1lo + x1hi)
    lo[0] = x0lo; hi[0] = x0hi; ac[0] = nc
    lo[1] = x1lo; hi[1] = x1hi; ac[1] = Ec
    lo[2] = x2lo; hi[2] = x2hi; ac[2] = 0.5 * (x2lo + x2hi)
    emax = S0 * x0hi * x1hi; emin = S0 * x0lo * x1lo; ecen = S0 * nc * Ec
    lo[3] = math.sqrt(1 - emax * emax); hi[3] = math.sqrt(1 - emin * emin); ac[3] = math.sqrt(1 - ecen * ecen)
    lo[4] = 1 / hi[3]; hi[4] = 1 / lo[3]; ac[4] = 1 / ac[3]
    return ac, lo, hi
''',
    'L2T': '''
@njit(cache=True)
def atoms_box(x0lo, x0hi, x1lo, x1hi, x2lo, x2hi, x3lo, x3hi):
    # atoms nu, E, P, A, b, ib ; box (nu, E, P, A); e = S0 E nu >= 0
    lo = np.empty(6); hi = np.empty(6); ac = np.empty(6)
    nc = 0.5 * (x0lo + x0hi); Ec = 0.5 * (x1lo + x1hi)
    lo[0] = x0lo; hi[0] = x0hi; ac[0] = nc
    lo[1] = x1lo; hi[1] = x1hi; ac[1] = Ec
    lo[2] = x2lo; hi[2] = x2hi; ac[2] = 0.5 * (x2lo + x2hi)
    lo[3] = x3lo; hi[3] = x3hi; ac[3] = 0.5 * (x3lo + x3hi)
    emax = S0 * x0hi * x1hi; emin = S0 * x0lo * x1lo; ecen = S0 * nc * Ec
    lo[4] = math.sqrt(1 - emax * emax); hi[4] = math.sqrt(1 - emin * emin); ac[4] = math.sqrt(1 - ecen * ecen)
    lo[5] = 1 / hi[4]; hi[5] = 1 / lo[4]; ac[5] = 1 / ac[4]
    return ac, lo, hi
''',
}
KIND = {'B1': ('B1', 'aff'), 'B2': ('B2', 'aff'), 'B3p': ('B3', 'mv4'), 'B3m': ('B3', 'mv4'), 'L2E': ('L2E', 'aff'),
        'L2P': ('L2X', 'aff'), 'L2Ap': ('L2X', 'aff'), 'L2Am': ('L2X', 'aff'), 'L2Tp': ('L2T', 'mv4'),
        'L2Tm': ('L2T', 'mv4'),
        # D2b corner charts (gen_bu D2*): same atom layouts with e = SD2 * (...)
        'D2E': ('L2E', 'aff'), 'D2E2': ('L2E', 'aff'), 'D2P': ('L2X', 'aff'), 'D2K': ('L2X', 'aff'), 'D2Tp': ('L2T', 'mv4'), 'D2Tm': ('L2T', 'mv4')}
SCALE_E = {'D2E': 5e-3, 'D2E2': 5e-3, 'D2P': 5e-3, 'D2K': 5e-3, 'D2Tp': 5e-3, 'D2Tm': 5e-3}

HEAD = '''"""Generated by make_engines_bu.py: certificate engine for blow-up chart {chart} ({kind}).
Enclosures as in nb_cert.py: centre values + chain-rule gradient at the centre atoms; the variation of each gradient
component over the box is bounded by Pabs(|a_c| + Da) - Pabs(|a_c|) (absolute-value polynomial); |coef| <= Pabs(|a_c| + Da).
Arcs are padded by ARC_PAD (overlapping arcs; their union covers [-pad, 2 pi + pad])."""
import math
import numpy as np
from numba import njit
from nbcoef_{chart} import coef_val, coef_abs
from nb_cert import arc_signs, dp_count

LIMIT = 11          # parity lemma: an even count <= 11 is <= 10
ARC_PAD = 1e-12
S0 = {S0!r}
{atoms_box}
'''

AFF = '''
@njit(cache=True)
def box_coeffs(x0lo, x0hi, x1lo, x1hi, x2lo, x2hi, x3lo, x3hi):
    ac, lo, hi = atoms_box(x0lo, x0hi, x1lo, x1hi, x2lo, x2hi, x3lo, x3hi)
    na = ac.size
    da = np.empty(na); aa = np.empty(na); ab = np.empty(na)
    for i in range(na):
        da[i] = max(hi[i] - ac[i], ac[i] - lo[i]) * (1 + 1e-12) + 1e-15 * (abs(ac[i]) + 1.0)
        aa[i] = abs(ac[i]); ab[i] = aa[i] + da[i]
    val = np.zeros((4, 2, 2, 7)); a0 = np.zeros((4, 2, 2, 7)); a1 = np.zeros((4, 2, 2, 7))
    coef_val(ac, val); coef_abs(aa, a0); coef_abs(ab, a1)
    cmid = np.empty((2, 2, 7)); crad = np.empty((2, 2, 7)); cmax = np.empty((2, 2, 7))
    dcen = np.empty((3, 2, 2, 7)); dwid = np.empty((3, 2, 2, 7))
    for f in range(2):
        for pt in range(2):
            for k in range(7):
                cmid[f, pt, k] = val[0, f, pt, k]
                crad[f, pt, k] = 1e-13 * a0[0, f, pt, k] + 1e-300
                cmax[f, pt, k] = a1[0, f, pt, k] * (1 + 1e-12)
                for v in range(3):
                    dcen[v, f, pt, k] = val[1 + v, f, pt, k]
                    dwid[v, f, pt, k] = (a1[1 + v, f, pt, k] - a0[1 + v, f, pt, k]) * (1 + 1e-12) + 1e-13 * a1[1 + v, f, pt, k]
    rr = np.empty(3)
    rr[0] = 0.5 * (x0hi - x0lo) * (1 + 1e-12) + 1e-15 * (abs(0.5 * (x0lo + x0hi)) + 1.0)
    rr[1] = 0.5 * (x1hi - x1lo) * (1 + 1e-12) + 1e-15 * (abs(0.5 * (x1lo + x1hi)) + 1.0)
    rr[2] = 0.5 * (x2hi - x2lo) * (1 + 1e-12) + 1e-15 * (abs(0.5 * (x2lo + x2hi)) + 1.0)
    return cmid, crad, dcen, dwid, cmax, rr


@njit(cache=True)
def signs(cmid, crad, dcen, dwid, cmax, rr, x3lo, x3hi, um, hh, J):
    return arc_signs(cmid, crad, dcen, dwid, cmax, rr, x3lo, x3hi, um, hh, J)
'''

MV4 = '''
@njit(cache=True)
def box_coeffs(x0lo, x0hi, x1lo, x1hi, x2lo, x2hi, x3lo, x3hi):
    ac, lo, hi = atoms_box(x0lo, x0hi, x1lo, x1hi, x2lo, x2hi, x3lo, x3hi)
    na = ac.size
    da = np.empty(na); aa = np.empty(na); ab = np.empty(na)
    for i in range(na):
        da[i] = max(hi[i] - ac[i], ac[i] - lo[i]) * (1 + 1e-12) + 1e-15 * (abs(ac[i]) + 1.0)
        aa[i] = abs(ac[i]); ab[i] = aa[i] + da[i]
    val = np.zeros((5, 1, 2, 7)); a0 = np.zeros((5, 1, 2, 7)); a1 = np.zeros((5, 1, 2, 7))
    coef_val(ac, val); coef_abs(aa, a0); coef_abs(ab, a1)
    cmid = np.empty((2, 7)); crad = np.empty((2, 7)); cmax = np.empty((2, 7))
    dcen = np.empty((4, 2, 7)); dwid = np.empty((4, 2, 7))
    for pt in range(2):
        for k in range(7):
            cmid[pt, k] = val[0, 0, pt, k]
            crad[pt, k] = 1e-13 * a0[0, 0, pt, k] + 1e-300
            cmax[pt, k] = a1[0, 0, pt, k] * (1 + 1e-12)
            for v in range(4):
                dcen[v, pt, k] = val[1 + v, 0, pt, k]
                dwid[v, pt, k] = (a1[1 + v, 0, pt, k] - a0[1 + v, 0, pt, k]) * (1 + 1e-12) + 1e-13 * a1[1 + v, 0, pt, k]
    rr = np.empty(4)
    rr[0] = 0.5 * (x0hi - x0lo) * (1 + 1e-12) + 1e-15 * (abs(0.5 * (x0lo + x0hi)) + 1.0)
    rr[1] = 0.5 * (x1hi - x1lo) * (1 + 1e-12) + 1e-15 * (abs(0.5 * (x1lo + x1hi)) + 1.0)
    rr[2] = 0.5 * (x2hi - x2lo) * (1 + 1e-12) + 1e-15 * (abs(0.5 * (x2lo + x2hi)) + 1.0)
    rr[3] = 0.5 * (x3hi - x3lo) * (1 + 1e-12) + 1e-15 * (abs(0.5 * (x3lo + x3hi)) + 1.0)
    return cmid, crad, dcen, dwid, cmax, rr


@njit(cache=True)
def signs(cmid, crad, dcen, dwid, cmax, rr, x3lo, x3hi, um, hh, J):
    """same bounds as nb_cert.arc_signs, for a single function with four mean-value variables."""
    M = um.size
    pos = np.ones((J + 1, M), dtype=np.bool_)
    neg = np.ones((J + 1, M), dtype=np.bool_)
    hopeless = np.ones(M, dtype=np.bool_)
    c = cmid; cr = crad; dc = dcen; dw = dwid; cm = cmax
    for i in range(M):
        u = um[i]; h = hh[i]
        vc = np.empty(J + 2); rad = np.empty(J + 2)
        for j in range(J + 2):
            s = 0.0; mag = 0.0; rd = 0.0; aab = 0.0
            g = np.zeros(4); w = np.zeros(4)
            for k in range(7):
                if j > 0 and k == 0:
                    continue
                kj = float(k) ** j
                ab_k = abs(c[0, k]) + abs(c[1, k]) + cr[0, k] + cr[1, k]
                for v in range(4):
                    ab_k += rr[v] * (abs(dc[v, 0, k]) + abs(dc[v, 1, k]) + dw[v, 0, k] + dw[v, 1, k])
                aab += kj * ab_k
                T = math.cos(k * u + j * math.pi / 2) * kj
                U = math.sin(k * u + j * math.pi / 2) * kj
                aT = abs(T); aU = abs(U)
                s += c[0, k] * T + c[1, k] * U
                mag += abs(c[0, k]) * aT + abs(c[1, k]) * aU
                rd += cr[0, k] * aT + cr[1, k] * aU
                for v in range(4):
                    g[v] += dc[v, 0, k] * T + dc[v, 1, k] * U
                    w[v] += dw[v, 0, k] * aT + dw[v, 1, k] * aU
            prad = 0.0
            for v in range(4):
                prad += rr[v] * (abs(g[v]) + w[v])
            tot = rd + prad
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
            tail = tail + 1e-15 * (abs(vc[j]) + tail)
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
'''

CERT = '''

@njit(cache=True)
def _run(x0lo, x0hi, x1lo, x1hi, x2lo, x2hi, x3lo, x3hi, J, M0, rounds, arc_cap, stall_rule):
    cmid, crad, dcen, dwid, cmax, rr = box_coeffs(x0lo, x0hi, x1lo, x1hi, x2lo, x2hi, x3lo, x3hi)
    left = np.arange(M0) * (2 * math.pi / M0)
    width = np.full(M0, 2 * math.pi / M0)
    total = 999
    prev = 10 ** 9
    stall = 0
    for rd in range(rounds + 1):
        sg, hopeless = signs(cmid, crad, dcen, dwid, cmax, rr, x3lo, x3hi, left + width / 2, width / 2 + ARC_PAD, J)
        total = dp_count(sg)
        if total <= LIMIT:
            return True, total
        if rd == rounds or hopeless or left.size > arc_cap:
            break
        if stall_rule:
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
def certify_box(x0lo, x0hi, x1lo, x1hi, x2lo, x2hi, x3lo, x3hi, J, M0, rounds):
    return _run(x0lo, x0hi, x1lo, x1hi, x2lo, x2hi, x3lo, x3hi, J, M0, rounds, 4096, True)


@njit(cache=True)
def certify_box_deep(x0lo, x0hi, x1lo, x1hi, x2lo, x2hi, x3lo, x3hi, J, M0, rounds, arc_cap):
    return _run(x0lo, x0hi, x1lo, x1hi, x2lo, x2hi, x3lo, x3hi, J, M0, rounds, arc_cap, False)
'''

for chart in sys.argv[1:]:
    ab, kind = KIND[chart]
    src = HEAD.format(chart=chart, kind=kind, S0=SCALE_E.get(chart, S0), atoms_box=ATOMS[ab]) + (AFF if kind == 'aff' else MV4) + CERT
    open(f'nbeng_{chart}.py', 'w').write(src)
    print('written', f'nbeng_{chart}.py')

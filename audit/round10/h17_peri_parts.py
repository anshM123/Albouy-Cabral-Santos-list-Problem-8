"""AUDIT-10 h17: each of the four parts of kq_run.count_peri bounds the zeros of P in its region, and the total bounds
N_mult, at 100 digits (MY OWN P).  Random boxes (e in [0.95, 1], l, w, t in [0.5, 1)); for each box the parts are
computed with the production engines exactly as count_peri calls them (aphelion KBQe_2 on |th| <= tan(phi0/2) + 1e-9,
perihelion KBPe_2 on |tp| <= tan(peri/2) + 1e-9, nb_cert2.away_count_arc2 on the two arcs; rounds 16, cap 8192,
target 0), and at 4 points of the box the 100-digit zeros of P are counted per (unpadded) region.
usage: python h17_peri_parts.py NBOX SEED"""
import os, sys, math, random
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.normpath(os.path.join(HERE, '..', '..', 'code'))
sys.path.insert(0, CODE); sys.path.insert(0, os.path.join(CODE, 'c3')); os.chdir(CODE)
import kq_run
if __name__ == '__main__':
    NBOX, SEED = int(sys.argv[1]), int(sys.argv[2]); random.seed(SEED)
    kq_run.init(peri=kq_run.PHI_P)
mp.mp.dps = 100
dn, up = kq_run.dn, kq_run.up


def zeros_f(e, l, w, t):
    e, l, w, t = [mp.mpf(v) for v in (e, l, w, t)]
    pz, qz = l * mp.sin(w), l * mp.cos(w)
    N = 32; vals = []
    for n in range(N):
        f = 2 * mp.pi * n / N
        zeta = pz * mp.cos(f) + qz * mp.sin(f); dz = -pz * mp.sin(f) + qz * mp.cos(f)
        B = -2 * zeta * dz + 2 * e * (mp.sin(f) - qz * zeta); D = 1 + e * mp.cos(f); k = 2 * e * (1 + e)
        vals.append((1 - t) * k ** 2 * mp.sin(f) ** 2 * (1 - zeta ** 2) - t * D ** 2 * B ** 2)
    c = [mp.fsum(vals[n] * mp.expj(-kk * 2 * mp.pi * n / N) for n in range(N)) / N for kk in range(-6, 7)]
    sc = max(abs(v) for v in c); lo = 0
    while abs(c[lo]) < mp.mpf(10) ** -80 * sc:
        lo += 1
    hi = 12
    while abs(c[hi]) < mp.mpf(10) ** -80 * sc:
        hi -= 1
    z = mp.polyroots(c[lo:hi + 1][::-1], maxsteps=1500, extraprec=1500) if hi > lo else []
    return [mp.arg(r) % (2 * mp.pi) for r in z if abs(mp.log(abs(r))) < mp.mpf(10) ** -25]


if __name__ == '__main__':
    viol = 0; nbox = 0; npt = 0; tight = {'A': 0, 'P': 0, '1': 0, '2': 0}
    for it in range(NBOX):
        wd = [10 ** random.uniform(-4, -2) for _ in range(4)]
        e0 = random.uniform(0.95, 1 - wd[0]); l0 = random.uniform(0, 1 - wd[1]); w0 = random.uniform(0, math.pi / 2 - wd[2])
        t0 = random.uniform(0.5, 0.9999 - wd[3])
        bx = (e0, e0 + wd[0], l0, l0 + wd[1], w0, w0 + wd[2], t0, t0 + wd[3])
        phi0 = kq_run.phi_split(bx); peri = kq_run.PHI_P
        e_0, e_1, l_0, l_1, w_0, w_1, t_0, t_1 = bx
        pz0 = dn(l_0 * dn(math.sin(w_0))); pz1 = up(l_1 * up(math.sin(w_1)))
        qz0 = dn(l_0 * dn(math.cos(w_1))); qz1 = up(l_1 * up(math.cos(w_0)))
        ep0 = max(dn(1.0 - e_1), 0.0); ep1 = up(1.0 - e_0)
        mu0 = max(dn(math.sqrt(max((1.0 - t_1) / t_1, 0.0))), 0.0); mu1 = up(math.sqrt((1.0 - t_0) / t_0))
        t_hi = math.tan(phi0 / 2) + 1e-9; tp_hi = math.tan(peri / 2) + 1e-9
        arcl = math.pi - phi0 - peri + 2e-9
        _, qA = kq_run.cnt(ep0, ep1, pz0, pz1, qz0, qz1, mu0, mu1, -t_hi, t_hi, kq_run.J, kq_run.M0, 16, 8192, 0)
        _, qP = kq_run.cntp(ep0, ep1, pz0, pz1, qz0, qz1, mu0, mu1, -tp_hi, tp_hi, kq_run.J, kq_run.M0, 16, 8192, 0)
        _, q1 = kq_run.awarc(*bx, math.pi + phi0 - 1e-9, arcl, 16, 16, kq_run.J, 0.0, 8192)
        _, q2 = kq_run.awarc(*bx, peri - 1e-9, arcl, 16, 16, kq_run.J, 0.0, 8192)
        okt, tot = kq_run.count_peri(bx, 16, 8192, phi0, peri)
        nbox += 1
        for k in range(4):
            x = [0.5 * (bx[2 * i] + bx[2 * i + 1]) if k == 0 else random.uniform(bx[2 * i], bx[2 * i + 1]) for i in range(4)]
            fz = zeros_f(*x); npt += 1
            cA = sum(1 for f in fz if abs(f - mp.pi) <= phi0)
            cP = sum(1 for f in fz if min(f, 2 * mp.pi - f) <= peri)
            c1 = sum(1 for f in fz if mp.pi + phi0 <= f <= 2 * mp.pi - peri)
            c2 = sum(1 for f in fz if peri <= f <= mp.pi - phi0)
            for nm, c, b in (('A', cA, qA), ('P', cP, qP), ('1', c1, q1), ('2', c2, q2)):
                if c > b:
                    viol += 1; print('   PART VIOLATION', nm, bx, x, c, b, flush=True)
                if c == b:
                    tight[nm] += 1
            if okt and len(fz) > tot:
                viol += 1; print('   TOTAL VIOLATION', bx, x, len(fz), tot, flush=True)
    print(f'{nbox} boxes x 4 points = {npt}: violations {viol}; parts attained (count == bound) A {tight["A"]}, '
          f'P {tight["P"]}, arc1 {tight["1"]}, arc2 {tight["2"]}')

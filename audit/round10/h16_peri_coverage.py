"""AUDIT-10 h16: coverage of the circle by the four parts of kq_run.count_peri, with the float end points the code uses
(t_hi = tan(phi0/2) + 1e-9, tp_hi = tan(peri/2) + 1e-9, arcs [pi + phi0 - 1e-9, 2 pi - peri + 1e-9] and
[peri - 1e-9, pi - phi0 + 1e-9] padded by PAD inside away_count_arc2, every polynomial root interval and every arc
widened by ARC_PAD by the engines).  For each (phi0, peri) the overlaps at the four seams are computed exactly (mpmath)
and must be positive; also in face mode (aphelion part only on |th| >= T_X)."""
import math, random
import mpmath as mp
mp.mp.dps = 40
PAD = ARC_PAD = 1e-12
random.seed(16)


def parts(phi0, peri):
    t_hi = math.tan(phi0 / 2) + 1e-9
    tp_hi = math.tan(peri / 2) + 1e-9
    arcl = math.pi - phi0 - peri + 2e-9
    A = 2 * mp.atan(mp.mpf(t_hi) + ARC_PAD)                         # aphelion: |f - pi| <= A
    Pp = 2 * mp.atan(mp.mpf(tp_hi) + ARC_PAD)                       # perihelion: |f| <= Pp
    a1 = math.pi + phi0 - 1e-9; a2 = peri - 1e-9
    arc1 = (mp.mpf(a1) - PAD - ARC_PAD, mp.mpf(a1) + mp.mpf(arcl) + PAD + ARC_PAD)
    arc2 = (mp.mpf(a2) - PAD - ARC_PAD, mp.mpf(a2) + mp.mpf(arcl) + PAD + ARC_PAD)
    return A, Pp, arc1, arc2, phi0 + peri >= math.pi


worst = mp.inf; n = 0
cases = [(0.9, 1.2), (1.25, 1.2)] + [(random.uniform(0.3, 1.6), random.uniform(0.3, 1.6)) for _ in range(2000)]
for phi0, peri in cases:
    A, Pp, arc1, arc2, drop = parts(phi0, peri)
    if drop:
        ov = [A + Pp - mp.pi]                                        # aphelion and perihelion parts must meet
    else:
        ov = [mp.pi - A - arc2[0] if False else (arc2[1] - (mp.pi - A)),     # arc2 reaches into the aphelion part
              Pp - arc2[0],                                          # perihelion part reaches into arc2
              (mp.pi + A) - arc1[0],                                 # aphelion part reaches into arc1
              arc1[1] - (2 * mp.pi - Pp)]                            # arc1 reaches into the perihelion part (at 2 pi)
    worst = min(worst, min(ov)); n += 1
    if min(ov) <= 0:
        print('  GAP', phi0, peri, [mp.nstr(v, 5) for v in ov])
print(f'{n} (phi0, peri) pairs incl. (0.9, 1.2) and (1.25, 1.2): smallest seam overlap {mp.nstr(worst, 5)} rad')
for phi0 in (0.9, 1.25):
    A, Pp, arc1, arc2, drop = parts(phi0, 1.2)
    print(f'  phi0 = {phi0}, peri = 1.2: overlaps {mp.nstr(arc2[1] - (mp.pi - A), 4)}, {mp.nstr(Pp - arc2[0], 4)}, '
          f'{mp.nstr((mp.pi + A) - arc1[0], 4)}, {mp.nstr(arc1[1] - (2 * mp.pi - Pp), 4)}')
# face mode: the aphelion part counts only |th| >= T_X = 0.06; the inner part |th| < 0.06 belongs to the face charts
T_X = 3.0 * 0.02
inner = 2 * mp.atan(mp.mpf(T_X) - ARC_PAD)
print(f'face mode: aphelion part starts at |f - pi| >= {mp.nstr(inner, 12)} = 2 atan(T_X - ARC_PAD); the face charts cover '
      f'|th| <= 3 SF = {3 * 0.02!r} (gap/near), so the seam overlaps by ARC_PAD in th')

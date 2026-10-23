"""AUDIT-4 d4: the winding-evaluation step of wind_box (lines 138-160 of wind_cert.py, replicated verbatim) on a crafted
closed polygon of node values whose arcs satisfy the labels (sanity checks pass) and whose TRUE winding is -1:
  arc0 label 3 (Re<0):  (-eps,-1) -> (-eps, 1)  through the left half-plane   (true step -pi + 2eps, clockwise)
  arc1 label 2 (Im>0):  (-eps, 1) -> (1, eps)                                   (true step ~ -pi/2)
  arc2 label 1 (Re>0):  (1, eps) -> (eps, -1)                                   (true step ~ -pi/2)
  arc3 label 4 (Im<0):  (eps, -1) -> (-eps, -1)                                 (true step ~ 0)
With eps below the resolution of atan2 at pi/2, the principal difference of arc0 is exactly fl(pi) = math.pi, which the
strict normalisation loop does not move to -pi: the computed winding is 0 instead of -1 (passes the 0.5 test).
Also shows the per-label branch evaluation that fixes it."""
import math
def wind_code(nodes, labs):
    tot = 0.0
    re0, im0 = nodes[0]; a_prev = math.atan2(im0, re0)
    M = len(labs)
    for i in range(M):
        re1, im1 = nodes[(i + 1) % M]
        lb = labs[i]
        for (rv, iv) in ((re0, im0), (re1, im1)):
            if (lb == 1 and not rv > 0) or (lb == 2 and not iv > 0) or (lb == 3 and not rv < 0) or (lb == 4 and not iv < 0):
                return None
        a1 = math.atan2(im1, re1); d = a1 - a_prev
        while d > math.pi: d -= 2 * math.pi
        while d < -math.pi: d += 2 * math.pi
        tot += d; a_prev = a1; re0, im0 = re1, im1
    w = int(round(tot / (2 * math.pi)))
    return (w if abs(tot - 2 * math.pi * w) <= 0.5 else None), tot
def branch_arg(v, lb):          # continuous branch of arg on the open half-plane of label lb
    a = math.atan2(v[1], v[0])
    if lb == 3 and a < 0: a += 2 * math.pi          # Re<0: (pi/2, 3pi/2)
    if lb == 4 and a > 0: a -= 2 * math.pi          # Im<0: (-pi, 0)  (atan2 already there)
    return a
def wind_fixed(nodes, labs):
    tot = 0.0; M = len(labs)
    for i in range(M):
        tot += branch_arg(nodes[(i + 1) % M], labs[i]) - branch_arg(nodes[i], labs[i])
    return round(tot / (2 * math.pi)), tot
eps = 1e-30
nodes = [(-eps, -1.0), (-eps, 1.0), (1.0, eps), (eps, -1.0)]
labs = [3, 2, 1, 4]
print("wind_box step logic: (winding, total) =", wind_code(nodes, labs), " -> W = 6 + winding")
print("per-label branch logic:            ", wind_fixed(nodes, labs))
print("true winding: -1 (clockwise loop around 0)")

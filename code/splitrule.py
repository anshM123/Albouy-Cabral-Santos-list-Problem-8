"""Split rules shared by kq_run.py and the replay tools (check_all.py, bb_replay.py).  The split of a box is a
deterministic function of the box, the header scale and the header split rule, so the replays rebuild the same trees.

  None / 'max'  (default; files without 'split_rule' in the header): bisect the coordinate with the largest normalised
                width (b[2i+1] - b[2i]) / scale[i] (first maximum).
  'face'        (kq_run.py --facesplit): KQ boxes (e, l, w, t) near the face e = 1, mu = 0, where eps = 1 - e and
                mu = sqrt((1 - t)/t) (C3_CORE.md section 6).
                (1) A box that meets {eps < S2, mu < S2}, S2 = FACE_S^2 (1 - 5e-9), and has e = 1 - S2 (or else
                    t = 1/(1 + S2^2)) strictly inside its range is cut there; after at most two such cuts the child on the face side lies in
                    the face neighbourhood {eps <= FACE_S^2, mu <= FACE_S^2} and is pruned ('p', covered by the face
                    charts and the KQF run).
                (2) Otherwise the widths are measured in the blow-up scale of the face (weights eps:2, mu:2): the eps
                    and mu widths relative to s2 = max(eps_hi, mu_hi, FACE_S^2), the l and w widths normalised by the
                    scale; the coordinate with the largest such width is bisected at its midpoint (t at the midpoint
                    of t).
  'face/K'      (kq_run.py --splitrule=face/K, K > 0): as 'face', with the eps and mu widths of (2) divided by K (so
                l and w are split about K times finer than eps and mu).
  'shell/K'     as 'face/K' but without the cut (1) (run KQS, whose boxes all lie in the face neighbourhood: the cut
                would only produce slivers there).
Any split rule gives a valid certificate; the rule only decides which boxes are certified.
"""
import math

FACE_S = 0.02                    # kq_run.FACE_S, c3_run.SF, check_all.FACE_S


def mu_of(t):
    return math.sqrt(max((1.0 - t) / t, 0.0))


def is_face(rule):
    return rule is not None and (rule == 'face' or rule.startswith('face/') or rule.startswith('shell/'))


def face_k(rule):
    return 1.0 if rule == 'face' else float(rule.split('/', 1)[1])


def rel_widths(bx, scale, rule=None):
    """the widths that the split rule compares (also used by kq_run for its minimum-width and deep-pass tests)"""
    if not is_face(rule):
        return [(bx[2 * i + 1] - bx[2 * i]) / scale[i] for i in range(4)]
    e0, e1, l0, l1, w0, w1, t0, t1 = bx
    mu_hi = mu_of(t0); mu_lo = mu_of(t1)
    s2 = max(1.0 - e0, mu_hi, FACE_S * FACE_S) * face_k(rule)
    return [(e1 - e0) / s2, (l1 - l0) / scale[1], (w1 - w0) / scale[2], (mu_hi - mu_lo) / s2]


def split_box(bx, scale, rule=None):
    """the two children of bx under the given rule"""
    i = None
    if is_face(rule) and not rule.startswith('shell/'):
        S2 = FACE_S * FACE_S * (1 - 5e-9)
        e_star = 1.0 - S2
        t_star = 1.0 / (1.0 + S2 * S2)
        if (1.0 - bx[1]) < S2 and mu_of(bx[7]) < S2:
            if bx[0] < e_star < bx[1]:
                i, mid = 0, e_star
            elif bx[6] < t_star < bx[7]:
                i, mid = 3, t_star
    if i is None:
        w = rel_widths(bx, scale, rule)
        i = 0
        for k in range(1, 4):
            if w[k] > w[i]:
                i = k
        mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx)
    a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)

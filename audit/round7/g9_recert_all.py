"""AUDIT-7 g9 (final pass): exact re-certification of sampled leaves of every bb_run / kq_run family with the
production engines, my own dispatch (not bb_replay's):
  regular leaves 'a'..'l'  -> the engine's certify_box (same J, M0, rounds as the header): must return ok, i.e. a
                              whole-circle total <= LIMIT = 11 (<= 10 by parity); KA (c3_away.away_count_t as in
                              bb_run): bound <= the recorded bound (<= 2);
  deep leaves 'A'..'L'     -> certify_box_deep (header J, M0, deep_rounds, deep_cap): ok as above.
  The number of leaves whose current total differs from the recorded one is reported separately ('not exact'): runs
  made before the parity change (LIMIT 10 -> 11, 03:51) were refined down to <= 10, the final code stops at <= 11.
  winding leaves '1'..'6'  -> wind_cert.wind_box on the engine's box_coeffs with the stored L and the header's wind
                              parameters: ok and the same W;
  KQ count leaves          -> kq_run.certify (current code); KQ winding leaves -> nb_cert chart 1 (f-chart), then nb_cert2.
  Engines a, f: if the first-order nb_cert certifier fails, nb_cert2 (bb_run --order2) is tried and counted separately.
N per leaf kind and family (N < 0: all leaves of that kind).
usage: python g9_recert_all.py N SEED NAME [NAME ...] [--kinds reg,deep,wind]"""
import os, sys, json, gzip, math, random, importlib
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '..', 'numba_cache'))
sys.dont_write_bytecode = True
CODE = os.path.normpath(os.path.join(HERE, '..', '..', 'code'))
argv = sys.argv[1:]
KINDS = ('reg', 'deep', 'wind')
if '--kinds' in argv:
    i = argv.index('--kinds'); KINDS = tuple(argv[i + 1].split(',')); argv = argv[:i] + argv[i + 2:]
N, SEED = int(argv[0]), int(argv[1]); names = argv[2:]
sys.path.insert(0, CODE); sys.path.insert(0, os.path.join(CODE, 'c3'))
os.chdir(CODE)
import wind_cert
import splitrule
random.seed(SEED)
ENG = {'bulk': 'a', 'F0': 'f', 'KQ1': 'KQ', 'KQ2': 'KQ', 'KQK1': 'KQ', 'KQK2': 'KQ', 'KQB': 'KQ'}


def load(name):
    for p in ('runs/%s.json' % name, 'runs/%s.json.gz' % name):
        if os.path.exists(p):
            return json.load(gzip.open(p, 'rt') if p.endswith('.gz') else open(p))
    return None


SPLIT_RULE = None                                    # header split_rule of the file being replayed (round 9)


def split(bx, S):
    if SPLIT_RULE:
        return splitrule.split_box(bx, S, SPLIT_RULE)
    w = [(bx[2 * i + 1] - bx[2 * i]) / S[i] for i in range(4)]
    i = max(range(4), key=lambda j: w[j]); mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
    a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
    return tuple(a), tuple(c)


for name in names:
    D = load(name)
    if D is None:
        print(name, 'not present'); continue
    H = D['header']; S = H['scale']; eng = H['engine']
    SPLIT_RULE = H.get('split_rule')
    groups = {'reg': [], 'deep': [], 'wind': []}; seen = {'reg': 0, 'deep': 0, 'wind': 0}

    def keep(kind, item):                     # reservoir sample of size N (all items if N < 0)
        if kind not in KINDS:
            return
        seen[kind] += 1
        g = groups[kind]
        if N < 0 or len(g) < N:
            g.append(item)
        else:
            j = random.randrange(seen[kind])
            if j < N:
                g[j] = item
    for r in D['results']:
        stack = [tuple(r['box'])]; wl = r.get('wL', []); wi = 0
        for ch in r['tree']:
            bx = stack.pop()
            if ch == '+':
                stack.extend(split(bx, S))
            elif 'a' <= ch <= 'l':
                keep('reg', (bx, ord(ch) - ord('a'), None))
            elif 'A' <= ch <= 'L':
                keep('deep', (bx, ord(ch) - ord('A'), None))
            elif '1' <= ch <= '6':
                keep('wind', (bx, int(ch), wl[wi])); wi += 1
    J, M0 = H.get('J', 6), H.get('M0', 64)
    reg2 = deep2 = wind2 = None                 # second-order fallbacks (nb_cert2; bb_run --order2 since 09:26)
    if eng in ('a', 'f'):
        import nb_cert, nb_cert2
        chn = 0 if eng == 'a' else 1
        reg = lambda bx: nb_cert.certify_box(chn, *bx, J, M0, H['rounds'])
        deep = lambda bx: nb_cert.certify_box_deep(chn, *bx, J, M0, H['deep_rounds'], H['deep_cap'])
        bcf = lambda bx: nb_cert.box_coeffs(chn, *bx[:6])
        reg2 = lambda bx: nb_cert2.certify_box2(chn, *bx, J, M0, H['rounds'])
        deep2 = lambda bx: nb_cert2.certify_box_deep2(chn, *bx, J, M0, H['deep_rounds'], H['deep_cap'])
        wind2 = lambda bx, L: nb_cert2.wind_box2(*nb_cert2.box_data2(chn, *bx[:6]), bx[6], bx[7], L, H.get('wind_m0', 64),
                                                 H.get('wind_rounds', 3), H.get('wind_cap', 1024))
    elif eng == 'KQ':
        import kq_run, nb_cert
        kq_run.init()
        if hasattr(kq_run, 'DEEP_ALWAYS'):
            kq_run.DEEP_ALWAYS = True                # round 8: two-tier counts; allow the deep pass
        reg = kq_run.certify; deep = None
        bcf = lambda bx: nb_cert.box_coeffs(1, *bx[:6])
        # AUDIT-12: KQ winding leaves are certified with nb_cert2's enclosure (kq_run.wind: box_coeffs2 + wind_box2,
        # header 'enclosure'), which is never wider than nb_cert's; re-check first with nb_cert, then with nb_cert2
        import nb_cert2
        wind2 = lambda bx, L: nb_cert2.wind_box2(*nb_cert2.box_coeffs2(*bx[:6]), bx[6], bx[7], L, H.get('wind_m0', 64),
                                                 H.get('wind_rounds', 3), H.get('wind_cap', 1024))
    elif eng == 'KQF':                        # round 9: outer count of the face construction
        import kq_run
        kq_run.init(face=True, faces=H.get('FACE_S', 0.02), tx=H.get('T_X', 0.06))   # round 11: KQF2 counts from 3 SF2
        kq_run.DEEP_ALWAYS = True
        reg = kq_run.certify_outer; deep = None; bcf = None
    elif eng == 'KA':
        import c3_away
        PHI_A = H['phi_a']
        reg = lambda bx: c3_away.away_count_t(bx[0], bx[1], bx[2], bx[3], bx[4], bx[5], 1.0 / (1.0 + bx[7] * bx[7]),
                                              1.0 / (1.0 + bx[6] * bx[6]), PHI_A, 64, 16, J, 2.0, 8192)
        deep = None; bcf = None
    else:
        m = importlib.import_module('nbeng_' + eng)
        reg = lambda bx: m.certify_box(*bx, J, M0, H['rounds'])
        deep = lambda bx: m.certify_box_deep(*bx, J, M0, H['deep_rounds'], H['deep_cap'])
        bcf = (lambda bx: m.box_coeffs(*bx[:6])) if eng in ('C1', 'C1p', 'C1pp') else (lambda bx: m.box_coeffs(*bx))
    out = []; nweak = [0]; nexact = [0]; n2 = [0]
    for kind, items in groups.items():
        if kind not in KINDS:
            continue
        smp = items
        bad = 0
        for bx, b, L in smp:
            if kind == 'wind':
                ok, W = wind_cert.wind_box(*bcf(bx), bx[6], bx[7], L, H.get('wind_m0', 64), H.get('wind_rounds', 3),
                                           H.get('wind_cap', 1024))
                good = ok and W == b
                if not good and wind2 is not None:
                    ok, W = wind2(bx, L); good = ok and W == b; n2[0] += good
            else:
                ok, tot = (reg if kind == 'reg' else deep)(bx)
                good = ok and (int(tot) <= b if eng == 'KA' else int(tot) <= 11)
                if not good and reg2 is not None:
                    ok, tot = (reg2 if kind == 'reg' else deep2)(bx)
                    good = ok and int(tot) <= 11; n2[0] += good
                if good and int(tot) != b:
                    nexact[0] += 1
            if not good:
                bad += 1
                extra = ''
                if kind != 'wind':
                    extra = 'now (%s, %s)' % (bool(ok), int(tot))
                    if deep is not None and kind == 'reg':   # regular certifier failed now: try the deep one
                        okd, totd = deep(bx)
                        weak = okd and int(totd) <= 11
                        nweak[0] += 0 if weak else 1
                        if weak:
                            bad -= 1                      # certified by the final code (deep), total <= 11
                        extra += '; deep now (%s, %s) -> %s' % (bool(okd), int(totd), 'certified (<= 11)' if weak else 'NOT RE-CERTIFIED')
                    else:
                        nweak[0] += 1
                print('   MISMATCH' if 'certified (<= 11)' not in extra else '   (regular certifier failed, deep certifies)', name, kind, bx, 'recorded', b, L, extra, flush=True)
        out.append(f'{kind} {len(smp)}/{seen[kind]} (mismatches {bad})')
    print(f"{name} [{eng}]: " + '; '.join(out) + f'; certified but total != recorded (not exact): {nexact[0]}; certified only by the second-order engine: {n2[0]}; not re-certifiable even by the deep certifier: {nweak[0]}', flush=True)

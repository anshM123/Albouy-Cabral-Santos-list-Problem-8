"""AUDIT-15 j6: re-certify EVERY tau-winding ('7') and disc-winding ('8') leaf of the given run files.
KQ / KQF files: tau_wind.check_tw(bx, which, cre, cim, r) and tau_wind.check_disc(bx, header T_X);
c3 chart files (Ca, ...): tau_wind.check_tw_chart(chart, bx, cre, cim, r).
The leaf boxes are rebuilt by replaying each tree with the file's own header scale and split rule (c3: plain
bisection with the header scale), the contours taken from the record field 'tw' in tree order.
usage: python j6_recert_all_tw_disc.py RUN.json [...]"""
import sys, os, json, time
import numpy as np
CODE = r'<repo>\code'
sys.path.insert(0, CODE); sys.path.insert(0, os.path.join(CODE, 'c3'))
files = [os.path.abspath(f) for f in sys.argv[1:]]
os.chdir(CODE)
import splitrule
import tau_wind
assert tau_wind.DISC_W == 6, tau_wind.DISC_W
for fn in files:
    t0 = time.time()
    D = json.load(open(fn)); H = D['header']; S = H['scale']; rule = H.get('split_rule')
    eng = H.get('engine') or H.get('chart')
    is_c3 = 'chart' in H and 'engine' not in H

    def split(bx):
        if rule is not None:
            return splitrule.split_box(bx, S, rule)
        w = [(bx[2 * i + 1] - bx[2 * i]) / S[i] for i in range(4)]
        i = int(np.argmax(w)); mid = 0.5 * (bx[2 * i] + bx[2 * i + 1])
        a = list(bx); c = list(bx); a[2 * i + 1] = mid; c[2 * i] = mid
        return tuple(a), tuple(c)
    n7 = n8 = b7 = b8 = 0; probs = 0
    for r in D['results']:
        stack = [tuple(r['box'])]; tw = r.get('tw', []); k = 0
        for ch in r['tree']:
            bx = stack.pop()
            if ch == '+':
                stack.extend(split(bx)); continue
            if ch == '7':
                if k >= len(tw):
                    probs += 1; k += 1; continue
                which, cre, cim, rr = tw[k]; k += 1; n7 += 1
                if is_c3:
                    ok = which == 0 and tau_wind.check_tw_chart(eng, bx, cre, cim, rr)
                else:
                    ok = tau_wind.check_tw(bx, which, cre, cim, rr)
                if not ok:
                    b7 += 1; print('   TW RECERT FAILURE', os.path.basename(fn), bx, tw[k - 1], flush=True)
            elif ch == '8':
                n8 += 1
                if is_c3 or eng != 'KQF':
                    b8 += 1; print('   DISC LEAF IN A NON-KQF FILE', os.path.basename(fn), bx, flush=True); continue
                if not tau_wind.check_disc(bx, H['T_X']):
                    b8 += 1; print('   DISC RECERT FAILURE', os.path.basename(fn), bx, flush=True)
        if k != len(tw):
            probs += 1
    print(f"{os.path.basename(fn)} [{eng}]: tau-winding leaves {n7} re-certified, failures {b7}; disc-winding leaves {n8} "
          f"re-certified (T_X {H.get('T_X')!r}), failures {b8}; tw-list problems {probs}; {time.time() - t0:.0f} s", flush=True)

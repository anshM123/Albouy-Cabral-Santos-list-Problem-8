"""Summarise a certify5 run: completion, cost per region, and where the failed leaves are."""
import json
import math
import sys
import numpy as np

d = json.load(open(sys.argv[1]))
res = d["results"]
print("initial boxes:", len(res), " complete:", sum(r["done"] for r in res))
print("leaves certified:", sum(r["ok"] for r in res), " failed leaves:", sum(len(r["fails"]) for r in res))
hist = {}
for r in res:
    for k, v in r["hist"].items():
        hist[int(k)] = hist.get(int(k), 0) + v
print("certified counts histogram:", dict(sorted(hist.items())))
fails = [f for r in res for f in r["fails"]]
if fails:
    F = np.array([f[0] for f in fails])
    mid = np.stack([(F[:, 2 * i] + F[:, 2 * i + 1]) / 2 for i in range(4)], 1)
    print("failed leaves: e range [%.4f, %.4f], l range [%.4f, %.4f], w range [%.4f, %.4f], t range [%.5f, %.5f]" % (
        math.sin(mid[:, 0].min()), math.sin(mid[:, 0].max()), mid[:, 1].min(), mid[:, 1].max(),
        mid[:, 2].min(), mid[:, 2].max(), mid[:, 3].min(), mid[:, 3].max()))
    for m in mid[:12]:
        print("   e=%.4f l=%.4f w=%.4f t=%.5f" % (math.sin(m[0]), m[1], m[2], m[3]))
# cost per initial box: certified volume fraction within the time budget
vol = np.array([r["vol"] for r in res])
bvol = np.array([np.prod([b[2 * i + 1] - b[2 * i] for i in range(4)]) for b in [r["box"] for r in res]])
frac = vol / bvol
print("per-initial-box certified fraction: min %.3f median %.3f mean %.3f" % (frac.min(), np.median(frac), frac.mean()))
sec = np.array([r["sec"] for r in res])
rate = vol.sum() / sec.sum()
domvol = np.prod([d["domain"][2 * i + 1] - d["domain"][2 * i] for i in range(4)])
print("throughput: %.3g volume per core-second; whole domain needs ~%.0f core-hours at this rate" % (rate, domvol / rate / 3600))
# where is it slow? group by e and by t
for name, idx in (("e", 0), ("l", 1), ("t", 3)):
    keys = sorted(set(round(r["box"][2 * idx], 4) for r in res))
    line = []
    for k in keys:
        sel = [i for i, r in enumerate(res) if round(r["box"][2 * idx], 4) == k]
        line.append("%.2f:%.2f" % (math.sin(k) if idx == 0 else k, frac[sel].mean()))
    print("mean certified fraction by %s lower edge: %s" % (name, "  ".join(line)))

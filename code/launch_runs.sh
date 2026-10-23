#!/bin/sh
# Launch every certificate run of the Problem 8 (Q2) proof (final run set).  Outputs: runs/*.json (+ .log,
# .fails.jsonl) and c3/runs/*.json.  PY must point to a Python with numpy + numba.  Worker counts are what was used on
# a 24-thread machine; they do not affect the certificates.
PY="${PY:-python}"
HP=1.5707963267948966
mkdir -p runs c3/runs
run() { nohup "$PY" bb_run.py "$@" > "runs/$NAME.log" 2>&1 & }
# --- bulk a-chart: e in [0.199, 0.836164] (eta = asin e in [0.20034, 0.99025]) ---
# (the recorded bulk run used the mean-value engine for its first 3825 initial boxes and was then resumed with --order2,
#  the second-order engine of nb_cert2, winding rounds 8 / cap 16384; either engine gives valid certificates.  It was
#  closed with --resume and the exclusion below: its initial boxes with eta >= 0.99025 (e >= 0.836164) are covered by
#  F1, KQ0 and KQ0T instead, where the f-chart variables are faster)
NAME=bulk run --engine a --dom 0.20033740617957885 1.2535563099909826 0 1 0 $HP 0 1 --n 8 10 10 10 --tmax 7200 --workers 12 --order2 --exclude 0.9902515840381316 1.2535563099909826 0 1 0 $HP 0 1 --out runs/bulk.json
# --- e in [0.835, 0.9501]: f-chart for t <= 0.5, engine KQ for t >= 0.5 ---
NAME=F1 run --engine f --order2 --dom 0.835 0.9501 0 1 0 $HP 0 0.5 --n 2 10 10 2 --tmax 3600 --workers 3 --out runs/F1.json
nohup "$PY" kq_run.py 0.835 0.9501 0 1 0 $HP 0.5 0.9975 2 10 10 2 1200 6 runs/KQ0.json --scale 1 1 $HP 0.5 > runs/KQ0.log 2>&1 &
nohup "$PY" kq_run.py 0.835 0.9501 0 1 0 $HP 0.9975 1.0 2 10 10 1 1200 8 runs/KQ0T.json --scale 1 1 $HP 0.005 > runs/KQ0T.log 2>&1 &
# --- e < 0.2 charts; C1 without the D1a block (covered by B1, B2, B3p, B3m) ---
NAME=C1 run --engine C1 --dom 0 0.26 0 1.1071487177940904 0 $HP 0 1 --n 8 12 4 16 --exclude 0 0.0650001 0 0.1846 0 1.5708 0.375 0.625 --tmax 7200 --workers 2 --out runs/C1.json
NAME=C1p run --engine C1p --dom 0 0.2001 0.4 1 0 $HP 0 1 --n 4 4 4 4 --tmax 7200 --workers 2 --out runs/C1p.json
NAME=C1pp run --engine C1pp --dom 0 0.51 0 0.52 0 $HP 0 1 --n 4 4 4 4 --tmax 7200 --workers 1 --out runs/C1pp.json
# --- D1a cone charts (B1 without the D1b block, covered by the L2 charts) ---
NAME=B1 run --engine B1 --dom 0 0.07 0 1.75 0 1.75 0 1 --n 7 7 7 8 --exclude 0 0.0100001 0 0.2500001 0.7499999 1.2500001 0.3749999 0.6250001 --tmax 7200 --workers 1 --out runs/B1.json
NAME=B2 run --engine B2 --dom 0 0.2 0 0.3675 0 $HP 0 1 --n 4 4 4 4 --tmax 7200 --workers 1 --out runs/B2.json
NAME=B3p run --engine B3p --dom 0 0.13 0 0.5654 0 1.6154 0 $HP --n 4 4 4 4 --tmax 7200 --workers 1 --out runs/B3p.json
NAME=B3m run --engine B3m --dom 0 0.13 0 0.5654 0 1.6154 0 $HP --n 4 4 4 4 --tmax 7200 --workers 1 --out runs/B3m.json
# --- D1b level-2 cone charts ---
NAME=L2E run --engine L2E --dom 0 1.01 0 1.05 -1.05 1.05 0 1 --n 4 4 4 4 --tmax 7200 --workers 1 --out runs/L2E.json
NAME=L2P run --engine L2P --dom 0 1.01 0 1.05 -1.05 1.05 0 1 --n 4 4 4 4 --tmax 7200 --workers 1 --out runs/L2P.json
NAME=L2Ap run --engine L2Ap --dom 0 1.01 0 1.05 0 1.05 0 1 --n 4 4 4 4 --tmax 7200 --workers 1 --out runs/L2Ap.json
NAME=L2Am run --engine L2Am --dom 0 1.01 0 1.05 0 1.05 0 1 --n 4 4 4 4 --tmax 7200 --workers 1 --out runs/L2Am.json
NAME=L2Tp run --engine L2Tp --dom 0 1.01 0 1.05 0 1.05 -1.05 1.05 --n 4 4 4 4 --tmax 7200 --workers 1 --out runs/L2Tp.json
NAME=L2Tm run --engine L2Tm --dom 0 1.01 0 1.05 0 1.05 -1.05 1.05 --n 4 4 4 4 --tmax 7200 --workers 1 --out runs/L2Tm.json
# --- D2b corner (inside C1p: e -> 0, s ~ e, l -> 1, w -> 0); covers C1p's failed leaves ---
NAME=D2E run --engine D2E --dom 0 1.01 0 1.05 0 1.05 0 1 --n 4 4 4 4 --tmax 7200 --workers 1 --out runs/D2E.json
NAME=D2P run --engine D2P --dom 0 1.01 0 1.05 0 1.05 0 1 --n 4 4 4 4 --tmax 7200 --workers 1 --out runs/D2P.json
NAME=D2K run --engine D2K --dom 0 1.01 0 1.05 0 1.05 0 1 --n 4 4 4 4 --tmax 7200 --workers 1 --out runs/D2K.json
NAME=D2Tp run --engine D2Tp --dom 0 1.01 0 1.05 0 1.05 0 1.05 --n 4 4 4 4 --tmax 1800 --workers 4 --out runs/D2Tp.json
NAME=D2Tm run --engine D2Tm --dom 0 1.01 0 1.05 0 1.05 0 1.05 --n 4 4 4 4 --tmax 7200 --workers 1 --out runs/D2Tm.json
# D2E extended to tau'/e in [1.04, 1.28]: covers the part of D2Tp's degenerate box with e/tau' < 1/1.05 (check_all D2 rules)
NAME=D2E2 run --engine D2E2 --dom 0 0.14 0 0.35 0.24 0.68 0 1 --n 4 4 4 4 --tmax 3600 --workers 2 --out runs/D2E2.json
# --- e >= 0.95 (near-parabolic): f-chart for t <= 0.5; engine KQ (tan-half-angle polynomial in eps = 1 - e near aphelion
#     + f-chart away) for t >= 0.5, for the corner region e >= 0.99 (minus the band and the C3 cell), and for the band ---
NAME=F0 run --engine f --dom 0.9499 1.0 0 1 0 $HP 0 0.5 --n 4 5 5 4 --tmax 7200 --workers 2 --out runs/F0.json
nohup "$PY" kq_run.py 0.9499 1.0 0 1 0 $HP 0.5 0.9975 4 10 10 5 7200 8 runs/KQ1.json > runs/KQ1.log 2>&1 &
nohup "$PY" kq_run.py 0.9499 0.99 0 1 0 $HP 0.9975 1.0 4 10 10 1 7200 2 runs/KQ2.json > runs/KQ2.log 2>&1 &
nohup "$PY" kq_run.py 0.99 1.0 0 0.95 0 $HP 0.997456 1.0 2 19 10 2 7200 6 runs/KQK1.json > runs/KQK1.log 2>&1 &
nohup "$PY" kq_run.py 0.99 1.0 0.95 1.0 0 0.47123889803846897 0.997456 1.0 2 1 3 2 7200 1 runs/KQK2.json > runs/KQK2.log 2>&1 &
nohup "$PY" kq_run.py 0.99 1.0 0.95 1.0 0.47123889803846897 1.413716694115407 0.997456 1.0 2 5 6 2 7200 6 runs/KQB.json > runs/KQB.log 2>&1 &
# The recorded KQ2, KQK1, KQK2 and KQB runs (and the C3 chart Ca below) were completed in two stages: the base run was
# closed with TMAX = -1 --resume (each unfinished initial box stored as an open record), and the open boxes were
# continued with a t-first split scale (Ca: K-first):
#   $PY kq_run.py --boxes runs/KQ2_r1_boxes.json --scale 1 1 $HP 0.005 7200 5 runs/KQ2.r1.json
#   $PY kq_run.py --boxes runs/KQK1_r1_boxes.json --scale 1 1 $HP 0.0025 7200 5 runs/KQK1.r1.json      (same for KQK2, KQB)
#   (cd c3; $PY c3_run.py Ca --boxes runs/Ca_r1_boxes.json --scale 0.37 0.25 1.002 0.5 7200 6 runs/Ca.r1.json)
# where runs/NAME_r1_boxes.json lists the open boxes of the base run (key "boxes").  From 09:46 on, KQ1 (base closed
# after 991 complete initial boxes), KQ2, KQK1, KQK2, KQB and Ca continued in levels r1, r2, ... with a 20-minute cap
# per box (TMAX 1200): a box still open after the cap passes its remaining stack to the next level.  The levels are
# driven by chain.sh (next_level.py lists the open boxes of a finished level), e.g.
#   PY=$PY sh chain.sh kq KQ1 1 8 1200 1 1 $HP 0.5      PY=$PY sh chain.sh c3 Ca 2 6 1200 0.37 0.25 1.002 0.5
# From level r4 (KQK1, KQK2) and r5 (KQB) on, the KQ continuation levels use the split rule 'face' (splitrule.py: the
# face neighbourhood eps, mu <= 0.02^2 is cut off and pruned 'p', widths measured in the blow-up scale), and from 11:50
# on most levels use a 5-minute cap per box (TMAX 300), e.g.
#   EXTRA=--facesplit PY=$PY sh chain.sh kq KQK1 4 10 300 1 1 $HP 0.0025
# KQ0 and KQ0T continue the same way (scale 1 1 $HP 0.5 and 1 1 $HP 0.005).  From 12:12 on, KQK1 r4+ also uses the
# perihelion count (--peri, C3_CORE.md section 7), and from 12:34 on KQB r5+ uses the split rule face/30 (l and w
# split about 30 times finer than eps and mu); every level file records its options in its header:
#   EXTRA="--facesplit --peri" PY=$PY sh chain.sh kq KQK1 4 10 300 1 1 $HP 0.0025
#   EXTRA="--splitrule=face/30" PY=$PY sh chain.sh kq KQB 5 8 300 1 1 $HP 0.0025
# --- face e = 1, mu = 0 (C3_CORE.md section 6): near charts FE, FM; gap charts FIP, FIM; outer count KQF ---
(cd c3; for c in FE FM FIP FIM; do nohup "$PY" c3_run.py $c 2 4 8 8 1200 3 runs/$c.json > runs/$c.log 2>&1 & done)
nohup "$PY" kq_run.py 0.9996 1.0 0 1 0 $HP 0.9999998300000288 1.0 1 10 10 1 1200 4 runs/KQF.json --face --scale 0.0008 1 $HP 3.4e-7 > runs/KQF.log 2>&1 &
#   continuation levels: PY=$PY sh chain.sh c3 FE 0 3 1200 0.02002 1.001 1 1   (FM, FIP likewise; next_level.py drops
#   open boxes with s_lo > SF = 0.02, outside the needed range);  KQF continues with an l/w-first scale:
#   EXTRA=--face PY=$PY sh chain.sh kq KQF 0 8 300 0.024 1 $HP 1.02e-5
# --- corner K_lw = {l >= 0.88, w >= 1.29} of the face neighbourhood (C3_CORE.md section 6): face charts only for
#     s <= SF2 = 0.010005 there, outer count KQF2 from 3 SF2, and the shell SF2 < s <= 0.02 by the whole-circle run KQS ---
nohup "$PY" kq_run.py 0.99989989 1.0 0.88 1.0 1.29 $HP 0.9999999899 1.0 1 4 6 1 300 8 runs/KQF2.json --face --faces=0.010005 --peri --scale 0.003 1 $HP 3e-7 > runs/KQF2.log 2>&1 &
nohup "$PY" kq_run.py 0.9996 1.0 0.88 1.0 1.29 $HP 0.99999983 1.0 2 4 6 2 300 8 runs/KQS.json --scale 1 1 $HP 0.0025 --faces=0.010005 --c3prune --splitrule=shell/30 --peri > runs/KQS.log 2>&1 &
#   continuation levels: chain.sh with EXTRA="--face --faces=0.010005 --peri" (KQF2) and
#   EXTRA="--faces=0.010005 --c3prune --splitrule=shell/30 --peri" (KQS)
# --- C3 cell: away part ---
NAME=KA run --engine KA --dom 0.99 1.0 0.95 1.0 1.413716694115407 $HP 0 0.0505 --n 2 2 4 4 --tmax 7200 --workers 1 --minw 1e-5 --out runs/KA.json
# --- C3 core charts ---
cd c3
for c in PhiP PhiM L2aYp L2aYm L2KYp L2KYm L2aE L2am L2KE L2Km CK CE Cm Ca; do
  nohup "$PY" c3_run.py $c 8 4 4 4 7200 2 runs/$c.json > runs/$c.log 2>&1 &
done
cd ..
echo launched

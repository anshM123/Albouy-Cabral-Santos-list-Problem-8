#!/bin/sh
# Final verification pass over the recorded runs (run from code/ after all runs have finished):
#   1. check_all.py: exact replay of every split tree, coverage of every family, writes ../RESULTS.md;
#   2. bb_replay.py --recert: re-certification of random samples of leaves (all methods, incl. winding leaves), for
#      every run file of the bb_run and kq_run families (base and continuation files);
#   3. verify_mp.py: 60-digit spot checks of random points of random leaves against the original polynomial;
#   4. recert_tw_c3.py: every tau-winding leaf of the c3 charts; disc_mp_check.py: 60-digit outer counts at random
#      points of the disc-winding leaves (code '8') of KQF and KQF2.
# Outputs go to ../verify/ (one file per run file for 2 and 3).  P parallel jobs (default 6).
PY="${PY:-python}"
N="${N:-300}"
P="${P:-6}"
export PY N
mkdir -p ../verify
"$PY" -B check_all.py > ../verify/check_all.out 2>&1
FAMS="bulk C1 C1p C1pp B1 B2 B3p B3m L2E L2P L2Ap L2Am L2Tp L2Tm D2E D2E2 D2P D2K D2Tp D2Tm F0 F1 KQ0 KQ0T KQ1 KQ2 KQK1 KQK2 KQB KA KQF KQF2 KQS"
for f in $FAMS; do ls runs/$f.json runs/$f.r*.json 2>/dev/null; done | \
  xargs -P "$P" -I{} sh -c 'b=$(basename {} .json); "$PY" -B bb_replay.py {} --recert "$N" --seed 23 > ../verify/recert_$b.out 2>&1'
for f in bulk C1 C1p C1pp B1 B2 B3p B3m L2E L2P L2Ap L2Am L2Tp L2Tm F0 F1 KQ0 KQ0T KQ1 KQ2 KQK1 KQK2 KQB KA KQS; do
  ls runs/$f.json runs/$f.r*.json 2>/dev/null; done | \
  xargs -P "$P" -I{} sh -c 'b=$(basename {} .json); "$PY" -B verify_mp.py {} --n 200 --seed 23 > ../verify/mp_$b.out 2>&1'
for c in Ca CK CE Cm; do ls c3/runs/$c.json c3/runs/$c.r*.json 2>/dev/null; done | \
  xargs -P "$P" -I{} sh -c 'b=$(basename {} .json); "$PY" -B verify_mp.py {} --n 200 --seed 23 > ../verify/mp_C3_$b.out 2>&1'
"$PY" -B recert_tw_c3.py Ca FE FM > ../verify/recert_tw_c3.out 2>&1      # every c3 tau-winding leaf (audit R13-1)
for f in KQF KQF2; do ls runs/$f.json runs/$f.r*.json 2>/dev/null; done | \
  xargs -P "$P" -I{} sh -c 'b=$(basename {} .json); "$PY" -B disc_mp_check.py {} --n 100 --seed 23 > ../verify/disc_mp_$b.out 2>&1'
echo done

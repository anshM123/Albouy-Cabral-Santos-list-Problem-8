#!/bin/sh
# AUDIT-12 final pass over the families whose files are FINAL (authors, 14:05): two serial chains (2 cores),
# outputs in audit/round7/final12/.  KQF, KQF2, KQS (paused, not final) are only replayed by g2, to report the 'p'
# dependency of KQ0, KQ1, KQK2 honestly.
PY="${PY:-python}"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
OUT="$HERE/final12"; mkdir -p "$OUT"
NUMBA_CACHE_DIR="$(cd ../numba_cache && pwd)"; export NUMBA_CACHE_DIR
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
R=../../code/runs; C=../../code/c3/runs
SMALL="F1 C1 C1p C1pp B1 B2 B3p B3m L2E L2P L2Ap L2Am L2Tp L2Tm D2E D2E2 D2P D2K D2Tp D2Tm F0 KA"
KQN="KQ0 KQ1 KQ1.r1 KQK2 KQK2.r1 KQK2.r2 KQK2.r3 KQK2.r4 KQK2.r5 KQK2.r6"
kqfiles=""; for f in $KQN; do kqfiles="$kqfiles $R/$f.json"; done
C3N="CK CE Cm PhiP PhiM L2aE L2am L2aYp L2aYm L2KE L2Km L2KYp L2KYm FE FM FIP FIM"
c3files=""; for c in $C3N; do c3files="$c3files $C/$c.json"; done
ONLY="bulk,F1,C1,C1p,C1pp,B1,B2,B3p,B3m,L2E,L2P,L2Ap,L2Am,L2Tp,L2Tm,D2E,D2E2,D2P,D2K,D2Tp,D2Tm,F0,KQ0,KQ1,KA,KQK2"
for c in $C3N; do case $c in F*) ONLY="$ONLY,face $c";; *) ONLY="$ONLY,C3 $c";; esac; done
ONLY="$ONLY,KQF,KQF2,KQS"
{
  date "+== chain A start %H:%M:%S"
  echo "== g2: own replay of the final families (+ KQF, KQF2, KQS for the p-dependency)"
  "$PY" -B g2_replay_all.py --prune-samples 6 --only "$ONLY"
  date "+== %H:%M:%S g9: EVERY leaf of the small whole-circle families, F0 and KA"
  "$PY" -B g9_recert_all.py -1 1 $SMALL
  date "+== %H:%M:%S g9: bulk, KQ0, KQ1(+r1), KQK2(+r1..r6): 300 regular/deep leaves per file"
  "$PY" -B g9_recert_all.py 300 23 bulk $KQN --kinds reg,deep,tw,disc
  date "+== %H:%M:%S g9: EVERY winding leaf of KQ0, KQ1(+r1), KQK2(+r1..r6)"
  "$PY" -B g9_recert_all.py -1 1 $KQN --kinds wind
  date "+== %H:%M:%S g9: bulk winding leaves, 20000 sampled"
  "$PY" -B g9_recert_all.py 20000 23 bulk --kinds wind
  date "+== %H:%M:%S g12: project check_all.py (unchanged; RESULTS.md redirected)"
  "$PY" -B g12_check_all_guarded.py "$OUT"
  date "+== chain A done %H:%M:%S"
} > "$OUT/chainA.out" 2>&1 &
{
  date "+== chain B start %H:%M:%S"
  echo "== g1: KQ0, KQ1(+r1), KQK2(+r1..r6): replay, every winding leaf, 100 count leaves per file, 100-digit checks"
  "$PY" -B g1_kq_wind.py -1 100 30 3 23 $kqfiles
  date "+== %H:%M:%S g20: D2E2 chart vs physics, D2Tp rules (i)/(ii)"
  "$PY" -B g20_d2e2_checks.py 40 23
  date "+== %H:%M:%S g3: 100-digit spot checks, whole-circle families + KA"
  "$PY" -B g3_spot_all.py 60 3 23 bulk $SMALL $KQN
  date "+== %H:%M:%S g4: final C3 charts + face charts"
  "$PY" -B g4_c3_spot.py 30 3 23 200 $c3files
  date "+== %H:%M:%S g13: winding step margins (F0 and the final KQ files; 20000 of bulk)"
  "$PY" -B g13_wind_margins.py -1 $R/F0.json $kqfiles
  "$PY" -B g13_wind_margins.py 20000 $R/bulk.json
  date "+== chain B done %H:%M:%S"
} > "$OUT/chainB.out" 2>&1 &
wait
grep -h -E "VIOLATION|MISMATCH|PROBLEM|NOT OK|NOT COMPLETE|BOUNDS: NO|not present|problems [1-9]|violations [1-9]|mismatches [1-9]|failures [1-9]|odd N_mult [1-9]|NOT IN ANY LEAF|Traceback|Error" "$OUT"/chain*.out > "$OUT/ALERTS.txt"
echo "final12 done; alerts: $(wc -l < "$OUT/ALERTS.txt") line(s)" >> "$OUT/chainA.out"

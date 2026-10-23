#!/bin/sh
# AUDIT-12 chain C (one core): the steps of final12.sh chain B after g1 (g1 keeps running alone in chain B).
PY="${PY:-python}"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
OUT="$HERE/final12"
NUMBA_CACHE_DIR="$(cd ../numba_cache && pwd)"; export NUMBA_CACHE_DIR
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
R=../../code/runs; C=../../code/c3/runs
SMALL="F1 C1 C1p C1pp B1 B2 B3p B3m L2E L2P L2Ap L2Am L2Tp L2Tm D2E D2E2 D2P D2K D2Tp D2Tm F0 KA"
KQN="KQ0 KQ1 KQ1.r1 KQK2 KQK2.r1 KQK2.r2 KQK2.r3 KQK2.r4 KQK2.r5 KQK2.r6"
kqfiles=""; for f in $KQN; do kqfiles="$kqfiles $R/$f.json"; done
C3N="CK CE Cm PhiP PhiM L2aE L2am L2aYp L2aYm L2KE L2Km L2KYp L2KYm FE FM FIP FIM"
c3files=""; for c in $C3N; do c3files="$c3files $C/$c.json"; done
{
  date "+== chain C start %H:%M:%S"
  date "+== %H:%M:%S g20: D2E2 chart vs physics, D2Tp rules (i)/(ii)"
  "$PY" -B g20_d2e2_checks.py 40 23
  date "+== %H:%M:%S g3: 100-digit spot checks, whole-circle families + KA"
  "$PY" -B g3_spot_all.py 60 3 23 bulk $SMALL $KQN
  date "+== %H:%M:%S g4: final C3 charts + face charts"
  "$PY" -B g4_c3_spot.py 30 3 23 200 $c3files
  date "+== %H:%M:%S g13: winding step margins (F0 and the final KQ files; 20000 of bulk)"
  "$PY" -B g13_wind_margins.py -1 $R/F0.json $kqfiles
  "$PY" -B g13_wind_margins.py 20000 $R/bulk.json
  date "+== %H:%M:%S g9: bulk winding leaves, another 60000 sampled (seed 41)"
  "$PY" -B g9_recert_all.py 60000 41 bulk --kinds wind
  date "+== chain C done %H:%M:%S"
} > "$OUT/chainC.out" 2>&1
grep -h -E "VIOLATION|MISMATCH|PROBLEM|NOT OK|NOT COMPLETE|BOUNDS: NO|not present|problems [1-9]|violations [1-9]|mismatches [1-9]|failures [1-9]|odd N_mult [1-9]|NOT IN ANY LEAF|Traceback|Error|WARNING" "$OUT"/chainC.out > "$OUT/ALERTS_C.txt"
echo "chain C finished; alerts: $(wc -l < "$OUT/ALERTS_C.txt")"

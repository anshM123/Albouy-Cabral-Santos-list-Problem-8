#!/bin/sh
# AUDIT-12 snapshot pass (after the authors closes the running levels, ~15:45): replay + sampled recertification of
# the families that are not final: KQ0T, KQ2, KQB, KQK1, KQF, KQF2, KQS (engine KQ / KQF) and C3 Ca.  2 cores.
PY="${PY:-python}"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
OUT="$HERE/snap12"; mkdir -p "$OUT"
NUMBA_CACHE_DIR="$(cd ../numba_cache && pwd)"; export NUMBA_CACHE_DIR
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
R=../../code/runs; C=../../code/c3/runs
FAMS="KQ0T KQ2 KQB KQK1 KQF KQF2 KQS"
names=""; files=""
for f in $FAMS; do
  for p in $R/$f.json $(ls $R/$f.r[0-9]*.json 2>/dev/null); do b=$(basename $p .json); names="$names $b"; files="$files $p"; done
done
cafiles="$C/Ca.json $(ls $C/Ca.r[0-9]*.json 2>/dev/null | tr '\n' ' ')"
ONLY="KQ0T,KQ2,KQB,KQK1,KQF,KQF2,KQS,C3 Ca,face FE,face FM,face FIP,face FIM"
{
  date "+== snapshot chain A start %H:%M:%S"
  echo "== g2: replay of the snapshot families (+ face charts for the p-dependency)"
  "$PY" -B g2_replay_all.py --prune-samples 4 --only "$ONLY"
  date "+== %H:%M:%S g9: 200 regular/deep leaves per file"
  "$PY" -B g9_recert_all.py 200 29 $names --kinds reg,deep
  date "+== %H:%M:%S g9: 2000 winding leaves per file"
  "$PY" -B g9_recert_all.py 2000 29 $names --kinds wind
  date "+== snapshot chain A done %H:%M:%S"
} > "$OUT/chainA.out" 2>&1 &
{
  date "+== snapshot chain B start %H:%M:%S"
  echo "== g1: replay + 40 count leaves per file, 15 spot leaves x 2 points at 100 digits (KQ and KQF engines)"
  "$PY" -B g1_kq_wind.py 0 40 15 2 29 $files
  date "+== %H:%M:%S g4: Ca + continuations: replay, 20 spot leaves x 2 points, 60 leaves re-certified per file"
  "$PY" -B g4_c3_spot.py 20 2 29 60 $cafiles
  date "+== snapshot chain B done %H:%M:%S"
} > "$OUT/chainB.out" 2>&1 &
wait
grep -h -E "VIOLATION|MISMATCH|PROBLEM|NOT OK|NOT COMPLETE|BOUNDS: NO|not present|problems [1-9]|violations [1-9]|mismatches [1-9]|failures [1-9]|odd N_mult [1-9]|NOT IN ANY LEAF|Traceback|Error" "$OUT"/chain*.out > "$OUT/ALERTS.txt"
echo "snap12 done; alerts: $(wc -l < "$OUT/ALERTS.txt") line(s)" >> "$OUT/chainA.out"

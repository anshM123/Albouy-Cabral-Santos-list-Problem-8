#!/bin/sh
# AUDIT-13: j1 on the TW leaves of the running levels (records.jsonl as of now), two serial chains
PY="python"
cd "$(dirname "$0")"
NUMBA_CACHE_DIR="$(cd ../numba_cache && pwd)"; export NUMBA_CACHE_DIR
R="<repo>/code/runs"
C="<repo>/code/c3/runs"
HP=1.5707963267948966
{
  date "+== chain 1 start %H:%M:%S"
  "$PY" -B j1_tw_leaves.py "$R/KQK1.r7.json.records.jsonl" kq face/30 1 1 $HP 0.0025 30 4 4 11
  "$PY" -B j1_tw_leaves.py "$R/KQB.r7.json.records.jsonl" kq face/30 1 1 $HP 0.0025 30 4 4 12
  "$PY" -B j1_tw_leaves.py "$R/KQ2.r5.json.records.jsonl" kq face/30 1 1 $HP 0.005 30 4 4 13
  "$PY" -B j1_tw_leaves.py "$R/KQ0T.r4.json.records.jsonl" kq face/30 1 1 $HP 0.005 30 4 4 14
  date "+== chain 1 done %H:%M:%S"
} > j1_chain1.out 2>&1 &
{
  date "+== chain 2 start %H:%M:%S"
  "$PY" -B j1_tw_leaves.py "$R/KQF.r2.json.records.jsonl" kq - 0.024 1 $HP 1.02e-5 30 4 4 21
  "$PY" -B j1_tw_leaves.py "$R/KQF2.r3.json.records.jsonl" kq - 0.003 1 $HP 3e-7 30 4 4 22
  "$PY" -B j1_tw_leaves.py "$R/KQS.r2.json.records.jsonl" kq shell/30 1 1 $HP 0.0025 30 4 4 23
  "$PY" -B j1_tw_leaves.py "$C/Ca.r8.json.records.jsonl" Ca - 0.37 0.25 1.002 0.5 30 4 4 24
  date "+== chain 2 done %H:%M:%S"
} > j1_chain2.out 2>&1 &
wait
echo "j1 done $(date +%H:%M:%S)"

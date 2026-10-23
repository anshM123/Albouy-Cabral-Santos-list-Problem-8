#!/bin/sh
# AUDIT-7 final pass over the FINAL run files (code/runs/*.json[.gz], code/c3/runs/*.json[.gz], incl. continuations).
# Run from audit/round7/ after every run has finished:   sh final_pass.sh
# Uses at most 2 cores (two serial chains), writes only into audit/round7/final/ (and audit/numba_cache).
# Expected wall time on this machine: about 1-2 h (dominated by bulk: replay, all ~1.5M winding leaves re-certified,
# and by the ~9k deep leaves of D2Tp).  (round6/f3 is not used: its exact-bound test has the R7-1 problem on D2E/P/K.)
PY="${PY:-python}"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
OUT="$HERE/final"; mkdir -p "$OUT"
NUMBA_CACHE_DIR="$(cd ../numba_cache && pwd)"; export NUMBA_CACHE_DIR
R=../../code/runs; C=../../code/c3/runs
BB="bulk F1 C1 C1p C1pp B1 B2 B3p B3m L2E L2P L2Ap L2Am L2Tp L2Tm D2E D2E2 D2P D2K D2Tp D2Tm F0 KQ0 KQ0T KQ1 KQ2 KQK1 KQK2 KQB KQF KA"
# continuation files (NAME.r1.json, ...: base runs closed, open boxes continued, possibly with another split scale)
CONT="$(cd $R && ls *.r[0-9]*.json 2>/dev/null | sed 's/\.json$//' | tr '\n' ' ')"
C3="Ca CK CE Cm PhiP PhiM L2aE L2am L2aYp L2aYm L2KE L2Km L2KYp L2KYm FE FM FIP FIM"
c3files=""; for c in $C3; do c3files="$c3files $C/$c.json"; done
c3files="$c3files $(ls $C/*.r[0-9]*.json 2>/dev/null | tr '\n' ' ')"
kqfiles="$R/KQ0.json $R/KQ0T.json $R/KQ1.json $R/KQ2.json $R/KQK1.json $R/KQK2.json $R/KQB.json $R/KQF.json $(ls $R/KQ*.r[0-9]*.json 2>/dev/null | tr '\n' ' ')"

# ---- chain 1: structure, re-certification, project checker ----
{
  echo "== g2: own replay of every family (grid, codes, per-family bounds, owed boxes, continuations, pruning)"
  "$PY" -B g2_replay_all.py --prune-samples 6
  echo "== g12: project check_all.py (unchanged; RESULTS.md redirected into final/)"
  "$PY" -B g12_check_all_guarded.py "$OUT"
  echo "== g9: EVERY leaf of every family except bulk and KQ*, final code, parity-aware (R7-1)"
  "$PY" -B g9_recert_all.py -1 1 F1 C1 C1p C1pp B1 B2 B3p B3m L2E L2P L2Ap L2Am L2Tp L2Tm D2E D2E2 D2P D2K D2Tp D2Tm F0 KA
  echo "== g9: bulk, KQ*, KQF and continuations: 300 regular/deep leaves per family, plus EVERY winding leaf"
  "$PY" -B g9_recert_all.py 300 23 bulk KQ0 KQ0T KQ1 KQ2 KQK1 KQK2 KQB KQF $CONT --kinds reg,deep
  "$PY" -B g9_recert_all.py -1 1 bulk KQ0 KQ0T KQ1 KQ2 KQK1 KQK2 KQB $CONT --kinds wind
} > "$OUT/chain1.out" 2>&1 &

# ---- chain 2: high-precision checks with my own formulas ----
{
  echo "== g1: KQ families + continuations (replay, every winding leaf re-certified, 300 count leaves, 100-digit checks)"
  "$PY" -B g1_kq_wind.py -1 300 40 3 23 $kqfiles
  echo "== g3: 100-digit spot checks, own chart maps, all whole-circle families + KA + continuations"
  "$PY" -B g3_spot_all.py 60 3 23 $BB $CONT
  echo "== g20: D2E2 chart vs physics, D2Tp rules (i) and (ii) with 100-digit checks in D2E / D2E2 leaves"
  "$PY" -B g20_d2e2_checks.py 40 23
  echo "== g4: C3 charts + continuations (own replay incl. zone codes up to F, 100-digit spots, 200 leaves re-certified)"
  "$PY" -B g4_c3_spot.py 30 3 23 200 $c3files
  echo "== g13: winding step margins (all winding leaves of F0 and KQ*, 20000 sampled of bulk)"
  "$PY" -B g13_wind_margins.py -1 $R/F0.json $kqfiles
  "$PY" -B g13_wind_margins.py 20000 $R/bulk.json
} > "$OUT/chain2.out" 2>&1 &
wait
sha256sum ../../code/runs/*.json* ../../code/c3/runs/*.json* 2>/dev/null | grep -v records > "$OUT/run_file_hashes.txt"
grep -h -E "VIOLATION|MISMATCH|PROBLEM|NOT OK|NOT COMPLETE|BOUNDS: NO|not present|problems [1-9]|violations [1-9]|mismatches [1-9]|failures [1-9]|odd N_mult [1-9]|NOT IN ANY LEAF|Traceback|Error" "$OUT"/chain*.out > "$OUT/ALERTS.txt"
echo "final pass done; alerts: $(wc -l < "$OUT/ALERTS.txt") line(s) in final/ALERTS.txt"

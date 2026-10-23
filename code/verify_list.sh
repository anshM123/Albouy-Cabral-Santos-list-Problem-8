#!/bin/sh
# Re-certification (bb_replay.py --recert) and 60-digit spot checks (verify_mp.py) for the given run families
# (base and continuation files), outputs in ../verify/ as in final_verify.sh.  usage: sh verify_list.sh FAMILY ...
PY="${PY:-python}"
N="${N:-300}"
P="${P:-6}"
export PY N
mkdir -p ../verify
for f in "$@"; do ls runs/$f.json runs/$f.r*.json 2>/dev/null; done | \
  xargs -P "$P" -I{} sh -c 'b=$(basename {} .json); "$PY" -B bb_replay.py {} --recert "$N" --seed 23 > ../verify/recert_$b.out 2>&1; "$PY" -B verify_mp.py {} --n 200 --seed 23 > ../verify/mp_$b.out 2>&1'
echo done

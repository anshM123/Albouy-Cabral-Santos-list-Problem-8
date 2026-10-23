#!/bin/sh
# core 2: after the final-family g1 (last file KQK2.r6) has printed its 100-digit line, run snapshot part C (g4 on Ca)
cd "$(dirname "$0")"
until [ "$(sed -n '/KQK2.r6.json: records/,$p' final12/chainB.out 2>/dev/null | grep -c '100-digit')" -ge 1 ]; do sleep 10; done
sleep 5
sh snap12c.sh
echo "queue12b finished $(date +%H:%M:%S)"

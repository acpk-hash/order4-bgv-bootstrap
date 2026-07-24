#!/bin/bash
# T2: p=257 order-4 measurement, m=2^16, d=256, 128 slots, i=14.
# OFF = generic bounded-support PS (HELIB_EXPLICIT_AUX=16)
# ON  = order-4 (HELIB_EXPLICIT_AUX=16 HELIB_AUX_ORDER4_EVAL=1)
# Sequential same-session OFF then ON. h passed as runtime arg. t=-1 newbts=1 newks=1 thick=0 repeat=3.
cd /home/user/experiments/baselines/bgv-bootstrapping-with-homomorphic-NTT
export HELIB_ZZX_CACHE_DIR=$PWD/cache/p257_order4
LOG=$PWD/results/p257_order4/driver.log
BIN=./build/fatboot
run_one() {
  local name=$1 hh=$2 ord4=$3
  echo "[$(date +%F_%T)] START $name h=$hh load=$(cut -d\" \" -f1-3 /proc/loadavg)" >> $LOG
  if [ "$ord4" = "1" ]; then
    HELIB_EXPLICIT_AUX=16 HELIB_AUX_ORDER4_EVAL=1 taskset -c 8-25 $BIN i=14 h=$hh t=-1 newbts=1 newks=1 thick=0 repeat=3 > results/p257_order4/$name.log 2>&1
  else
    HELIB_EXPLICIT_AUX=16 taskset -c 8-25 $BIN i=14 h=$hh t=-1 newbts=1 newks=1 thick=0 repeat=3 > results/p257_order4/$name.log 2>&1
  fi
  local rc=$?
  echo "[$(date +%F_%T)] END $name rc=$rc load=$(cut -d\" \" -f1-3 /proc/loadavg)" >> $LOG
  return $rc
}
H=${1:-7}
echo "[$(date +%F_%T)] P257 SEQUENCE BEGIN h=$H pid=$$" >> $LOG
run_one h${H}_OFF $H 0
run_one h${H}_ON  $H 1
echo "[$(date +%F_%T)] P257 DONE h=$H" >> $LOG

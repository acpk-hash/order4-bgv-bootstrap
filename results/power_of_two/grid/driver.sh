#!/bin/bash
cd /home/user/experiments/baselines/bgv-bootstrapping-with-homomorphic-NTT
export HELIB_ZZX_CACHE_DIR=$PWD/cache/po2_grid
LOG=$PWD/results/po2_grid/driver.log
run_one() {
  local name=$1 idx=$2 hh=$3 aux4=$4
  local mg
  mg=$(pgrep -c -f "magma.exe.*case_wrapper")
  echo "[$(date +%F_%T)] START $name load=$(cut -d" " -f1-3 /proc/loadavg) magma_running=$mg" >> $LOG
  if [ "$aux4" = "1" ]; then
    HELIB_EXPLICIT_AUX=256 HELIB_AUX_ORDER4_EVAL=1 taskset -c 8-25 ./build/fatboot i=$idx h=$hh t=-1 newbts=1 newks=1 thick=0 repeat=3 > results/po2_grid/$name.log 2>&1
    rc=$?
  else
    HELIB_EXPLICIT_AUX=256 taskset -c 8-25 ./build/fatboot i=$idx h=$hh t=-1 newbts=1 newks=1 thick=0 repeat=3 > results/po2_grid/$name.log 2>&1
    rc=$?
  fi
  echo "[$(date +%F_%T)] END $name rc=$rc load=$(cut -d" " -f1-3 /proc/loadavg)" >> $LOG
}
echo "[$(date +%F_%T)] GRID SEQUENCE BEGIN pid=$$ taskset=8-25" >> $LOG
for cfg in "i0_h26 0 26" "i3_h26 3 26" "i4_h26 4 26" "i0_h12 0 12" "i0_h18 0 18" "i0_h40 0 40"; do
  set -- $cfg
  run_one ${1}_OFF $2 $3 0
  run_one ${1}_ON  $2 $3 1
done
echo "[$(date +%F_%T)] GRID_DONE" >> $LOG

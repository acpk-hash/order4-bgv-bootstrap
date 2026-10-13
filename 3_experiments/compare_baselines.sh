#!/usr/bin/env bash
# One-click comparison against the original baselines, with control groups.
#
#   3_experiments/compare_baselines.sh            both experiments (~2.5 h, 32 GB, 2 cores)
#   3_experiments/compare_baselines.sh m50731     Ma et al. Large-p, p = 65537, m = 50731   (~1 h, 17 GB)
#   3_experiments/compare_baselines.sh po2        Ma et al. NTT,     p = 65537, m = 2^16    (~1.5 h, 32 GB)
#   3_experiments/compare_baselines.sh report     only rebuild the report from existing logs
#
# Each experiment runs a fixed set of groups, one bootstrapping run each
# (REPEAT trials, default 3), all decrypting and checking every slot:
#
#   G0  the ORIGINAL baseline, compiled from 2_baselines/ with none of our code
#   G1  our executable, evaluator OFF, the baseline's own radix
#         -> equivalence control: G0/G1 should be 1.00, i.e. our changes do
#            not alter the baseline when switched off
#   G2  our executable, evaluator OFF, order-four radix A
#         -> radix control: separates the effect of choosing A from the
#            effect of the evaluator
#   G3  our executable, evaluator ON, radix A      -> treatment
#   (m50731) G4 = G3 + parallel CoeffToSlot      -> the paper's "combined" row
#   (po2)    G4-G6 repeat G1-G3 at the recommended 80-bit chain (Table 4)
#
# Output: out/compare/<experiment>/<group>.log, and
#         out/compare/report.md, out/compare/summary.csv, out/compare/figures/
# produced by 3_experiments/compare_report.py.
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${BIN:-$ROOT/build/bin}"
OUT="${OUT:-$ROOT/out}"
REPEAT="${REPEAT:-3}"
export HELIB_ZZX_CACHE_DIR="${HELIB_ZZX_CACHE_DIR:-$ROOT/cache/saved_ZZX}"
mkdir -p "$HELIB_ZZX_CACHE_DIR"
TIME=(/usr/bin/time -v); [[ -x /usr/bin/time ]] || TIME=()

# group <experiment> <id> <binary> [VAR=value ...] -- <driver arguments>
group() {
  local exp=$1 id=$2 bin=$3; shift 3
  local envs=(HELIB_ZZX_CACHE_DIR="$HELIB_ZZX_CACHE_DIR")
  while [[ $1 != -- ]]; do envs+=("$1"); shift; done; shift
  [[ -x $bin ]] || { echo "missing $bin: run 3_experiments/build_all.sh" >&2; exit 2; }
  local d="$OUT/compare/$exp" log="$OUT/compare/$exp/$id.log"
  mkdir -p "$d"
  if [[ -s $log ]] && grep -q 'everything ok' "$log" && [[ ${FORCE:-0} != 1 ]]; then
    echo "[$(date +%T)] $exp/$id: kept existing log (FORCE=1 to rerun)"; return
  fi
  echo "[$(date +%T)] $exp/$id: ${envs[*]:1} $(basename "$bin") $*"
  ( cd "$d" && env "${envs[@]}" "${TIME[@]}" "$bin" "$@" ) > "$log" 2>&1
  if grep -q 'everything ok' "$log"; then
    grep 'time for linear1' "$log" | tail -1 | sed 's/^/           mean: /'
  else
    echo "           FAILED: no 'everything ok' marker (see $log)"
  fi
}

m50731() {
  local a=(i=4 h=12 t=-1 newbts=1 newks=1 thick=0 repeat="$REPEAT")
  group m50731 G0_original_ma_largep "$BIN/baseline_ma_largep" -- i=4 h=12 t=-1 newbts=1 newks=1 repeat="$REPEAT"
  group m50731 G1_ours_off_radix35   "$BIN/fatboot_order4" -- "${a[@]}"
  group m50731 G2_ours_off_radix256  "$BIN/fatboot_order4" HELIB_EXPLICIT_AUX=256 -- "${a[@]}"
  group m50731 G3_ours_on_radix256   "$BIN/fatboot_order4" HELIB_EXPLICIT_AUX=256 HELIB_AUX_ORDER4_EVAL=1 -- "${a[@]}"
  group m50731 G4_ours_on_parallel   "$BIN/fatboot_order4" HELIB_EXPLICIT_AUX=256 HELIB_AUX_ORDER4_EVAL=1 \
        HELIB_AUX_PARALLEL_COEFF2SLOT=1 -- "${a[@]}"
}

po2() {
  # upstream parameter set i=0: p = 65537, m = 2^16, the baseline's own 1332-bit chain
  local own=(i=0 h=26 t=-1 newbts=1 newks=1 thick=0 cache=1 repeat="$REPEAT")
  # the paper's recommended 80-bit chain (1066 bits built)
  local rec=(i=0 bits=782 mainhwt=120 h=26 t=-1 newbts=1 newks=1 thick=0 cache=1 repeat="$REPEAT")
  group po2 G0_original_ma_ntt         "$BIN/baseline_ma_ntt" -- "${own[@]}"
  group po2 G1_ours_off_radix27        "$BIN/fatboot_po2" -- "${own[@]}"
  group po2 G2_ours_off_radix256       "$BIN/fatboot_po2" HELIB_EXPLICIT_AUX=256 -- "${own[@]}"
  group po2 G3_ours_on_radix256        "$BIN/fatboot_po2" HELIB_EXPLICIT_AUX=256 HELIB_AUX_ORDER4_EVAL=1 -- "${own[@]}"
  group po2 G4_ours_off_radix27_rec    "$BIN/fatboot_po2" -- "${rec[@]}"
  group po2 G5_ours_off_radix256_rec   "$BIN/fatboot_po2" HELIB_EXPLICIT_AUX=256 -- "${rec[@]}"
  group po2 G6_ours_on_radix256_rec    "$BIN/fatboot_po2" HELIB_EXPLICIT_AUX=256 HELIB_AUX_ORDER4_EVAL=1 -- "${rec[@]}"
}

report() { python3 "$ROOT/3_experiments/compare_report.py" --dir "$OUT/compare"; }

case ${1:-all} in
  m50731) m50731; report ;;
  po2)    po2; report ;;
  all)    m50731; po2; report ;;
  report) report ;;
  *) sed -n 2,26p "$0"; exit 1 ;;
esac

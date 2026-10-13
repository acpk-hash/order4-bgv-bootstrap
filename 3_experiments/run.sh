#!/usr/bin/env bash
# Run the ciphertext experiments of Section 6 of the paper.
#
#   3_experiments/run.sh list                       targets, cost and memory
#   3_experiments/run.sh smoke                      one po2 configuration, OFF+ON, one trial (~10 min, 16 GB)
#   3_experiments/run.sh admissibility              the two witnessed admissibility failures (<1 min)
#   3_experiments/run.sh butterfly                  Section 5 encrypted butterfly prototype (~3 min, 5 GB)
#   3_experiments/run.sh po2     [p ...]            Table 5
#   3_experiments/run.sh general [set ...]          Table 6
#   3_experiments/run.sh m50731  [row ...]          Table 8
#   3_experiments/run.sh all                        everything above except smoke (~9.5 h)
#
# Logs are written to $OUT (default out/) as <experiment>/<run>.log|.out, in
# the same layout as the logs of the reported runs in 5_results/paper_runs/,
# so 3_experiments/make_tables.py --logs out tabulates them next to the paper.
#
# Environment: REPEAT (trials per run, default 3), OUT, BIN (default build/bin).
# Runs are sequential and single-threaded, except the m50731 'combined' row.
# Timings are sensitive to other load: use an otherwise idle machine.
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${BIN:-$ROOT/build/bin}"
OUT="${OUT:-$ROOT/out}"
REPEAT="${REPEAT:-3}"
export HELIB_ZZX_CACHE_DIR="${HELIB_ZZX_CACHE_DIR:-$ROOT/cache/saved_ZZX}"
mkdir -p "$HELIB_ZZX_CACHE_DIR"

TIME=(/usr/bin/time -v)
[[ -x /usr/bin/time ]] || TIME=()

# run <logfile> <binary> [VAR=value ...] -- <driver arguments>
run() {
  local log=$1 bin=$2; shift 2
  local envs=(HELIB_ZZX_CACHE_DIR="$HELIB_ZZX_CACHE_DIR")
  while [[ $1 != -- ]]; do envs+=("$1"); shift; done; shift
  [[ -x $bin ]] || { echo "missing $bin: run 3_experiments/build_all.sh first" >&2; return 2; }
  mkdir -p "$(dirname "$log")"
  echo "[$(date +%T)] ${log#"$OUT"/}: ${envs[*]:1} $(basename "$bin") $*"
  ( cd "$(dirname "$log")" && env "${envs[@]}" "${TIME[@]}" "$bin" "$@" ) > "$log" 2>&1
  local rc=$?
  if grep -q 'everything ok' "$log"; then
    grep 'time for linear1' "$log" | tail -1 | sed 's/^/           mean: /'
  else
    echo "           rc=$rc, no 'everything ok' marker (see $log)"
  fi
}

# ------------------------------------------------------------------ po2
# m = 2^16, recommended 80-bit chain. OFF is Ma et al.'s driver at its default
# radix A = 2B+1 = 27, exactly as measured for the paper (README, Sec. 8);
# ON uses the order-four radix A of the prime.
#      p      idx    A
PO2="  1297     9    36
       1601     5    40
       2521    13    71
       3137     6    56
       4513     8    95
       7057    10    84
       13457   11   116
       14401    7   120
       15377   12   124
       65537    0   256"

po2_one() {  # po2_one <p> [dir]
  local p=$1 dir=${2:-$OUT/po2} line idx a
  line=$(awk -v p="$p" '$1==p' <<< "$PO2")
  [[ -z $line ]] && { echo "p=$p is not in the po2 table; its order-four radix is unknown" >&2; return 1; }
  read -r _ idx a <<< "$line"
  local args=(i="$idx" bits=782 mainhwt=120 h=26 t=-1 newbts=1 newks=1 thick=0 cache=1 repeat="$REPEAT")
  run "$dir/R_p${p}_OFF.out" "$BIN/fatboot_po2" -- "${args[@]}"
  run "$dir/R_p${p}_ON.out"  "$BIN/fatboot_po2" HELIB_EXPLICIT_AUX="$a" HELIB_AUX_ORDER4_EVAL=1 -- "${args[@]}"
}

# ------------------------------------------------------------------ general
# Radix matched: OFF and ON use the same order-four radix A.
#      set   idx   h   bits    A
GEN="  80A     9  12   1101   95
       80B    12  12   1228  120
       80C    10  12   1318  116
       80D    11  12   1498   71
       128A   18  26   1055  116
       128B   14  26    856   54
       128C   15  26    892  124
       128D   16  26    918   90
       128E   17  26    839  120"

gen_one() {
  local s=$1 line idx h bits a
  line=$(awk -v s="$s" '$1==s' <<< "$GEN")
  [[ -z $line ]] && { echo "unknown set $s (one of: 80A-80D, 128A-128E)" >&2; return 1; }
  read -r _ idx h bits a <<< "$line"
  local args=(i="$idx" h="$h" bits="$bits" t=-1 newbts=1 newks=1 thick=0 repeat="$REPEAT")
  run "$OUT/general/${s}_off.log" "$BIN/fatboot_order4" HELIB_EXPLICIT_AUX="$a" -- "${args[@]}"
  run "$OUT/general/${s}_on.log"  "$BIN/fatboot_order4" HELIB_EXPLICIT_AUX="$a" HELIB_AUX_ORDER4_EVAL=1 -- "${args[@]}"
}

# ------------------------------------------------------------------ m50731
m50731() {  # m50731 [ma sparse order4 combined]
  local d=$OUT/m50731 args=(i=4 h=12 t=-1 newbts=1 newks=1 thick=0 repeat="$REPEAT")
  local rows=("$@"); [[ ${#rows[@]} -eq 0 ]] && rows=(ma sparse order4 combined)
  for r in "${rows[@]}"; do case $r in
    ma)       run "$d/ma_largep.log" "$BIN/baseline_ma_largep" -- \
                  i=4 h=12 t=-1 newbts=1 newks=1 repeat="$REPEAT" ;;
    sparse)   run "$d/sparse_aux256.log" "$BIN/fatboot_order4" HELIB_EXPLICIT_AUX=256 -- "${args[@]}" ;;
    order4)   run "$d/order4.log" "$BIN/fatboot_order4" HELIB_EXPLICIT_AUX=256 \
                  HELIB_AUX_ORDER4_EVAL=1 -- "${args[@]}" ;;
    # runs the two CoeffToSlot transforms concurrently: give it two cores
    combined) run "$d/combined.log" "$BIN/fatboot_order4" HELIB_EXPLICIT_AUX=256 \
                  HELIB_AUX_ORDER4_EVAL=1 HELIB_AUX_PARALLEL_COEFF2SLOT=1 -- "${args[@]}" ;;
    *)        echo "unknown m50731 row $r (one of: ma sparse order4 combined)" >&2 ;;
  esac; done
}

# ------------------------------------------------------------------ butterfly
# Section 5: the Step2 (CoeffToSlot) dim-0 block of m = 50731 evaluated on
# ciphertexts as a mixed-radix/Rader butterfly, checked slot by slot against
# the matrix HElib applies, followed by HElib's BSGS MatMul1D for reference.
butterfly() {
  local d=$OUT/butterfly bin=$BIN/enc_butterfly
  [[ -x $bin ]] || { echo "missing $bin: run 3_experiments/build_all.sh first" >&2; return 2; }
  mkdir -p "$d"
  echo "[$(date +%T)] butterfly/run_bits600.log: enc_butterfly 600 1 1"
  ( cd "$ROOT/1_implementation/butterfly" && "${TIME[@]}" "$bin" 600 1 1 ) > "$d/run_bits600.log" 2>&1
  grep -E 'ENCRYPTED CHAIN|VERIFY' "$d/run_bits600.log" | sed 's/^/           /'
}

# ------------------------------------------------------------------ admissibility
admissibility() {
  local d=$OUT/admissibility
  # Theorem 3 (ii): (2B+1)^2 < p fails at p = 257, h = 12 -> abort while building tables
  run "$d/p257_h12_abort.log" "$BIN/fatboot_order4" HELIB_EXPLICIT_AUX=35 HELIB_AUX_ORDER4_EVAL=1 -- \
      i=2 h=12 t=-1 newbts=1 newks=1 thick=0 repeat=1
  if grep -q 'p too small to hold the overflow part' "$d/p257_h12_abort.log"; then
    echo "           aborted as expected: p too small to hold the overflow part"
  else
    echo "           UNEXPECTED: abort message not found"
  fi
  # Theorem 3 (i): an order-four radix exists iff p = 1 mod 4 (exhaustive check)
  python3 - <<'EOF' | tee "$d/order_four_radix.txt"
for p in (8191, 257, 65537, 14401):
    roots = [a for a in range(p) if a * a % p == p - 1]
    print(f"           p={p:6d}  p mod 4 = {p % 4}  A with A^2 = -1 mod p: {roots or 'none'}")
EOF
}

list() {
  cat <<'EOF'
Measured on the paper's host, single-threaded, REPEAT=3 unless noted.

target            runs  wall clock   peak RAM          paper
smoke                2  ~10 min      15 GB             one row of Table 5, REPEAT=1
admissibility        1  <1 min       <2 GB             Sec. 6.3, Theorem 3
butterfly            1  ~3 min       5 GB              Sec. 5 (prototype, m = 50731)
po2     [p ..]      20  ~2.7 h       15-59 GB;         Table 5
                                     118 GB for p=2521
general [set ..]    18  ~5.6 h       8-61 GB           Table 6
m50731  [row ..]     4  ~1 h         17 GB, 2 cores    Table 8
                                     for 'combined'
all                 43  ~9.5 h       118 GB (p=2521)

Cheap subsets:  po2 1601 3137 14401    (15 GB, ~8 min per prime)
                general 80B 80A 128D   (8-13 GB, ~20 min per set)
                m50731 ma order4       (17 GB, ~30 min)
EOF
}

cmd=${1:-list}; shift || true
case $cmd in
  list)          list ;;
  smoke)         REPEAT=1 po2_one 14401 "$OUT/smoke/po2" ;;
  admissibility) admissibility ;;
  butterfly)     butterfly ;;
  po2)           for p in ${*:-1297 1601 2521 3137 4513 7057 13457 14401 15377 65537}; do po2_one "$p"; done ;;
  general)       for s in ${*:-80A 80B 80C 80D 128A 128B 128C 128D 128E}; do gen_one "$s"; done ;;
  m50731)        m50731 "$@" ;;
  all)           admissibility; butterfly; m50731; po2; general ;;
  *)             echo "unknown target '$cmd'" >&2; list; exit 1 ;;
esac

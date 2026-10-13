#!/usr/bin/env bash
# Run the ciphertext experiments of Section 7. Logs go to ae/out/<table>/,
# in the same format as the shipped logs under results/paper/, so that
# ae/make_tables.py can tabulate either.
#
#   ae/run.sh list                      what each target runs and costs
#   ae/run.sh smoke                     ~15 min, 20 GB: one prime on m=2^16, OFF+ON, 1 trial
#   ae/run.sh tab7                      Table 7 (m=50731): Ma et al., ours, combined
#   ae/run.sh tab5 [p ...]              Table 5 (m=2^16), all ten primes or the ones given
#   ae/run.sh tab6 [set ...]            Table 6 (general rings), all nine sets or the ones given
#   ae/run.sh admissibility             Section 7.3, the two witnessed failures
#
# REPEAT (default 3) sets the number of trials per run. Runs are sequential and
# single-threaded (except the "combined" row of Table 7, see README).
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="$ROOT/build/bin"
OUT="${OUT:-$ROOT/ae/out}"
REPEAT="${REPEAT:-3}"
export HELIB_ZZX_CACHE_DIR="$ROOT/cache/saved_ZZX"
mkdir -p "$HELIB_ZZX_CACHE_DIR"

TIME=(/usr/bin/time -v)
[[ -x /usr/bin/time ]] || TIME=()

# run <logfile> <binary> <env...> -- <args...>
run() {
  local log=$1 bin=$2; shift 2
  local envs=(HELIB_ZZX_CACHE_DIR="$HELIB_ZZX_CACHE_DIR")
  while [[ $1 != -- ]]; do envs+=("$1"); shift; done; shift
  mkdir -p "$(dirname "$log")"
  echo "[$(date +%T)] $(basename "$log"): ${envs[*]} $(basename "$bin") $*"
  ( cd "$(dirname "$log")" && env "${envs[@]}" "${TIME[@]}" "$bin" "$@" ) > "$log" 2>&1
  local rc=$?
  if grep -q 'everything ok' "$log"; then
    grep 'time for linear1' "$log" | tail -1 | sed 's/^/          /'
  else
    echo "          rc=$rc, no 'everything ok' marker (see $log)"
  fi
}

# ---------------------------------------------------------------- Table 5
#      p      idx  A
PO2="1297     9   36
     1601     5   40
     2521    13   71
     3137     6   56
     4513     8   95
     7057    10   84
     13457   11  116
     14401    7  120
     15377   12  124
     65537    0  256"

po2_one() {  # po2_one <p> [dir]
  local p=$1 dir=${2:-$OUT/tab5_po2} line idx a
  line=$(awk -v p="$p" '$1==p' <<< "$PO2")
  [[ -z $line ]] && { echo "p=$p is not in Table 5 (its order-four radix is unknown)" >&2; return 1; }
  read -r _ idx a <<< "$line"
  local args=(i="$idx" bits=782 mainhwt=120 h=26 t=-1 newbts=1 newks=1 thick=0 cache=1 repeat="$REPEAT")
  # OFF is the driver's default radix A = 2B+1 = 27, exactly as measured for
  # the paper (README, "Known deviations"); ON uses the order-four radix.
  run "$dir/R_p${p}_OFF.out" "$BIN/fatboot_po2" -- "${args[@]}"
  run "$dir/R_p${p}_ON.out"  "$BIN/fatboot_po2" HELIB_EXPLICIT_AUX="$a" HELIB_AUX_ORDER4_EVAL=1 -- "${args[@]}"
}

# ---------------------------------------------------------------- Table 6
#      set   idx  h   bits  A
GEN="80A     9  12  1101   95
     80B    12  12  1228  120
     80C    10  12  1318  116
     80D    11  12  1498   71
     128A   18  26  1055  116
     128B   14  26   856   54
     128C   15  26   892  124
     128D   16  26   918   90
     128E   17  26   839  120"

gen_one() {
  local s=$1 line idx h bits a
  line=$(awk -v s="$s" '$1==s' <<< "$GEN")
  [[ -z $line ]] && { echo "unknown set $s" >&2; return 1; }
  read -r _ idx h bits a <<< "$line"
  local args=(i="$idx" h="$h" bits="$bits" t=-1 newbts=1 newks=1 thick=0 repeat="$REPEAT")
  run "$OUT/tab6_general/${s}_off.log" "$BIN/fatboot_order4" HELIB_EXPLICIT_AUX="$a" -- "${args[@]}"
  run "$OUT/tab6_general/${s}_on.log"  "$BIN/fatboot_order4" HELIB_EXPLICIT_AUX="$a" HELIB_AUX_ORDER4_EVAL=1 -- "${args[@]}"
}

# ---------------------------------------------------------------- Table 7
tab7() {
  local d=$OUT/tab7_m50731 args=(i=4 h=12 t=-1 newbts=1 newks=1 thick=0 repeat="$REPEAT")
  run "$d/ma_largep.log" "$BIN/fatboot_ma_largep" -- i=4 h=12 t=-1 newbts=1 newks=1 repeat="$REPEAT"
  run "$d/sparse_aux256.log" "$BIN/fatboot_order4" HELIB_EXPLICIT_AUX=256 -- "${args[@]}"
  run "$d/order4.log" "$BIN/fatboot_order4" HELIB_EXPLICIT_AUX=256 HELIB_AUX_ORDER4_EVAL=1 -- "${args[@]}"
  run "$d/combined.log" "$BIN/fatboot_order4" HELIB_EXPLICIT_AUX=256 HELIB_AUX_ORDER4_EVAL=1 \
      HELIB_AUX_PARALLEL_COEFF2SLOT=1 -- "${args[@]}"
}

# ---------------------------------------------------------------- Sec 7.3
admissibility() {
  local d=$OUT/admissibility
  # (ii) (2B+1)^2 < p fails at p=257, h=12: aborts while building the tables
  run "$d/p257_h12_abort.log" "$BIN/fatboot_order4" HELIB_EXPLICIT_AUX=35 HELIB_AUX_ORDER4_EVAL=1 -- \
      i=2 h=12 t=-1 newbts=1 newks=1 thick=0 repeat=1
  grep -H 'p too small' "$d/p257_h12_abort.log" || echo "          (abort message not found)"
  # (i) p = 3 mod 4 has no A with A^2 = -1; checked exhaustively
  python3 - <<'EOF'
for p in (8191, 65537, 14401):
    roots = [a for a in range(p) if a * a % p == p - 1]
    print(f"          p={p:6d} p mod 4={p % 4}  order-four radices: {roots[:2] or 'none'}")
EOF
}

list() {
  cat <<'EOF2'
Measured on the paper's host, single-threaded, REPEAT=3 (wall clock per run):

target          what                                   runs   wall total   peak RAM
smoke           Table 5, p=14401, OFF+ON, REPEAT=1        2   ~6 min       15 GB
admissibility   p=257 abort + p mod 4 check               1   <1 min       <2 GB
tab7            Table 7, m=50731 (Ma, sparse, ON, comb.)  4   ~1 h         17 GB
tab5 [p ..]     Table 5, ten primes on m=2^16            20   ~2.7 h       15-59 GB, 118 GB for p=2521
tab6 [set ..]   Table 6, nine general rings              18   ~5.6 h       8-61 GB

Cheapest Table 5 primes: 1601, 3137, 14401 (15 GB, ~8 min per pair).
Cheapest Table 6 sets: 80B, 80A, 128D (8-13 GB, ~20 min per pair).
EOF2
}

cmd=${1:-list}; shift || true
case $cmd in
  list)          list ;;
  smoke)         REPEAT=1 po2_one 14401 "$OUT/smoke/tab5_po2" ;;
  tab5)          for p in "${@:-1297 1601 2521 3137 4513 7057 13457 14401 15377 65537}"; do for q in $p; do po2_one "$q"; done; done ;;
  tab6)          for s in "${@:-80A 80B 80C 80D 128A 128B 128C 128D 128E}"; do for t in $s; do gen_one "$t"; done; done ;;
  tab7)          tab7 ;;
  admissibility) admissibility ;;
  *)             echo "unknown target $cmd" >&2; list; exit 1 ;;
esac

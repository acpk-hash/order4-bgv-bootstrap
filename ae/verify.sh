#!/usr/bin/env bash
# Automated checks of the paper's claims, with a PASS/FAIL/SKIP line per test
# and a report in ae/out/verify_report.txt. Exit status 0 iff nothing failed.
#
#   ae/verify.sh            quick: no ciphertext benchmark          (~2 min,   2 GB)
#   ae/verify.sh smoke      quick + one po2 run + butterfly         (~15 min, 16 GB)
#   ae/verify.sh full       quick + every ciphertext experiment     (~9.5 h, 128 GB)
#
#   --security-rerun SET..  also recompute these security sets       (~20 min each, SageMath)
#   --lean                  also build the Lean proofs               (needs elan; ~10 min with cache)
#
# Test IDs match the README (Sec. 4, "Claims and how they are checked").
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
OUT="${OUT:-$ROOT/ae/out}"
export OUT
mkdir -p "$OUT"
REPORT="$OUT/verify_report.txt"
LEVEL=quick; SEC_SETS=(); LEAN=0
while [[ $# -gt 0 ]]; do
  case $1 in
    quick|smoke|full) LEVEL=$1 ;;
    --lean) LEAN=1 ;;
    --security-rerun) shift; while [[ $# -gt 0 && $1 != --* ]]; do SEC_SETS+=("$1"); shift; done; continue ;;
    -h|--help) sed -n 2,14p "$0"; exit 0 ;;
    *) echo "unknown argument $1" >&2; exit 2 ;;
  esac; shift
done

PASS=0; FAIL=0; SKIP=0
: > "$REPORT"
{ echo "order4-bgv artifact verification, level=$LEVEL, $(date -u +%FT%TZ)"
  echo "host: $(uname -srm), $(nproc) CPUs, $(awk '/MemTotal/{printf "%.0f GB", $2/1048576}' /proc/meminfo)"
  echo; } | tee -a "$REPORT"

# test <id> <description> <command...>: run, keep its output in $OUT/tests/<id>.log
test_() {
  local id=$1 desc=$2; shift 2
  mkdir -p "$OUT/tests"
  local t0=$SECONDS
  "$@" > "$OUT/tests/$id.log" 2>&1
  local rc=$? dt=$((SECONDS - t0))
  if [[ $rc -eq 0 ]]; then res=PASS; PASS=$((PASS + 1))
  elif [[ $rc -eq 77 ]]; then res=SKIP; SKIP=$((SKIP + 1))
  else res=FAIL; FAIL=$((FAIL + 1)); fi
  printf '%-4s %-5s %-62s %5ss\n' "$res" "$id" "$desc" "$dt" | tee -a "$REPORT"
  if [[ $res == FAIL ]]; then tail -5 "$OUT/tests/$id.log" | sed 's/^/           | /' | tee -a "$REPORT"; fi
}

have_bins() { [[ -x build/bin/fatboot_order4 && -x build/bin/fatboot_po2 && -x build/bin/fatboot_ma_largep ]]; }

t_env() {
  local cxx=${CXX:-g++}
  "$cxx" ae/check_ntl.cpp -o "$OUT/check_ntl" -lntl -lgmp -lpthread || return 1
  "$OUT/check_ntl" || return 1
  if have_bins; then ls -l build/bin; else echo "executables missing: run ae/build_all.sh"; return 1; fi
}

t_admissibility() {
  have_bins || return 77
  ae/run.sh admissibility || return 1
  grep -q 'p too small to hold the overflow part' "$OUT/admissibility/p257_h12_abort.log" || return 1
  grep -q 'p=  8191  p mod 4 = 3  A with A^2 = -1 mod p: none' "$OUT/admissibility/order_four_radix.txt"
}

t_smoke() {
  have_bins || return 77
  ae/run.sh smoke || return 1
  local on="$OUT/smoke/po2/R_p14401_ON.out" off="$OUT/smoke/po2/R_p14401_OFF.out"
  grep -q 'INFO: total bits = 1066,' "$off" || { echo "chain is not 1066 bits"; return 1; }
  grep -q 'HELIB_AUX_ORDER4_EVAL enabled: deg(P)=727' "$on" || { echo "order-four evaluator not engaged"; return 1; }
  python3 ae/make_tables.py --logs "$OUT/smoke" --check rerun po2
}

t_butterfly() {
  [[ -x build/bin/enc_butterfly ]] || return 77
  ae/run.sh butterfly || return 1
  local log="$OUT/butterfly/run_bits600.log"
  grep -q 'ENCRYPTED CHAIN: .* 41 rotations' "$log" || { echo "rotation count is not 41"; return 1; }
  grep -q 'VERIFY vs dump matvec (up to fixed Rader permutation): PASS  (2784 slots checked, 0 mismatches)' "$log"
}

t_full() {
  have_bins || return 77
  ae/run.sh all || return 1
  python3 ae/make_tables.py --logs "$OUT" --check rerun
}

t_lean() {
  command -v lake >/dev/null || { echo "lake not found"; return 77; }
  ( cd lean && lake exe cache get && lake build ) || return 1
  ! grep -rnw --include='*.lean' -e 'sorry' -e 'admit' lean/BootstrapAlgebra | grep -v -- '--\|/-\|-/\|`sorry`'
}

test_ E1 "environment: NTL configuration, executables built"          t_env
test_ P1 "plaintext: P_A = X/2 + X^3 Q(X^4) on every (p,A,B) used"     python3 ae/tests/test_plaintext.py
test_ L1 "shipped logs regenerate Tables 4-6 within 0.5%"              python3 ae/make_tables.py --check shipped
test_ S1 "security at built chain: transcripts equal paper (0.02 bit)" python3 ae/tests/test_security.py
test_ A1 "admissibility failures: p=257 abort, p=8191 no radix"       t_admissibility
if [[ ${#SEC_SETS[@]} -gt 0 ]]; then
  test_ S2 "security recomputed for: ${SEC_SETS[*]}"                   python3 ae/tests/test_security.py --rerun "${SEC_SETS[@]}"
fi
if [[ $LEAN -eq 1 ]]; then
  test_ F1 "Lean 4 formalisation builds, no sorry/admit"              t_lean
fi
if [[ $LEVEL == smoke || $LEVEL == full ]]; then
  test_ C1 "ciphertext smoke: po2 p=14401 OFF/ON, decrypts, >= 1.5x"   t_smoke
  test_ B1 "Sec. 5 butterfly: 41 rotations, 2784 slots, 0 mismatches"  t_butterfly
fi
if [[ $LEVEL == full ]]; then
  test_ C2 "ciphertext: all rows of Tables 4-6 rerun within 15%"       t_full
fi

{ echo; echo "summary: $PASS passed, $FAIL failed, $SKIP skipped   (logs: ${OUT#"$ROOT"/}/tests/)"; } | tee -a "$REPORT"
[[ $FAIL -eq 0 ]]

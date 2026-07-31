#!/bin/bash
# run_bench.sh - reproduce the reported ciphertext bootstrapping measurements.
#
#   ./run_bench.sh list
#   ./run_bench.sh po2-80 off
#   ./run_bench.sh po2-80 on
#   ./run_bench.sh po2-80 on 13457      # a different plaintext prime
#
# Why this is a wrapper and not a standalone driver.
#
# The order-four evaluator needs an auxiliary modulus that is not a power of the
# plaintext prime. Stock HElib chooses that modulus in RecryptData::setAE, which
# only considers powers of p, so on a power-of-two ring at p = 65537 it fails
# outright with "setAE: cannot find suitable e" before a context can even be
# built. The fork carried in baselines/ extends that choice, which is what the
# t = -1 flag and the HELIB_EXPLICIT_AUX variable below select. A driver written
# against stock HElib therefore cannot reproduce any of these numbers; it does
# not merely differ, it does not run.
#
# So the benchmark is the fork's own executable driven with the exact arguments
# the reported measurements used, and parameter.h is the record of what those
# arguments are for each set. Everything that varies between OFF and ON is in
# the environment, so the two differ only in the evaluator.
set -u

BIN=${FATBOOT:-../baselines/bgv-bootstrapping-with-homomorphic-NTT/build/fatboot}
TRIALS=${TRIALS:-3}

#   set        idx  p      A     h_enc  bits   built  target
SETS="
po2-80        0    65537  256   26     782    1066   80
po2-80-p1297  9     1297   36   26     782    1071   80
po2-80-p1601  5     1601   40   26     782    1071   80
po2-80-p2521 13     2521   71   26     782    1071   80
po2-80-p3137  6     3137   56   26     782    1071   80
po2-80-p4513  8     4513   95   26     782    1066   80
po2-80-p7057 10     7057   84   26     782    1066   80
po2-80-p13457 11   13457  116   26     782    1066   80
po2-80-p14401  7   14401  120   26     782    1066   80
po2-80-p15377 12   15377  124   26     782    1066   80
"

if [ "${1:-}" = "list" ] || [ $# -lt 2 ]; then
  echo "usage: $0 <set> <off|on> [prime]"
  echo
  printf "%-14s %4s %7s %5s %6s %7s %7s %7s\n" \
         set idx p A h_enc bits built target
  echo "$SETS" | while read -r n i p a h b bu t; do
    [ -z "${n:-}" ] && continue
    printf "%-14s %4s %7s %5s %6s %7s %7s %7s\n" "$n" "$i" "$p" "$a" "$h" "$b" "$bu" "$t"
  done
  echo
  echo "The m = 2^17 sets of the recommendation, po2-100 and po2-128, need a"
  echo "parameter-table entry for that ring which the baseline source does not"
  echo "carry; see README.md for the entry we added and validated."
  echo
  echo "General cyclotomic sets need their own mvec, gens and ords, which depend"
  echo "on m = q1*q2 and on ord_m(p). They live in the baseline's own table and"
  echo "are selected by index there, not by this script."
  exit 0
fi

NAME=$1; MODE=$2; WANT_P=${3:-}
LINE=$(echo "$SETS" | awk -v n="$NAME" '$1==n {print; exit}')
[ -z "$LINE" ] && { echo "unknown set '$NAME'; try: $0 list" >&2; exit 1; }
set -- $LINE
SETNAME=$1; IDX=$2; P=$3; A=$4; HENC=$5; BITS=$6; BUILT=$7

if [ -n "$WANT_P" ]; then
  L2=$(echo "$SETS" | awk -v p="$WANT_P" '$3==p {print; exit}')
  [ -z "$L2" ] && {
    echo "p=$WANT_P is not in the table. Its order-four radix would be unknown," >&2
    echo "and passing another prime's radix does not fail: the evaluator falls" >&2
    echo "back to the generic path and the run reports a speedup near 1.0, which" >&2
    echo "reads as a negative result rather than a misconfiguration." >&2
    exit 1; }
  set -- $L2
  SETNAME=$1; IDX=$2; P=$3; A=$4; HENC=$5; BITS=$6; BUILT=$7
fi

[ -x "$BIN" ] || { echo "fatboot not found at $BIN; set FATBOOT=<path>" >&2; exit 2; }

case "$MODE" in
  off) ENVV="HELIB_AUX_ORDER4_EVAL=0" ;;
  on)  ENVV="HELIB_EXPLICIT_AUX=$A HELIB_AUX_ORDER4_EVAL=1" ;;
  *)   echo "mode must be off or on" >&2; exit 1 ;;
esac

echo "set=$SETNAME p=$P A=$A idx=$IDX mode=$MODE h_main=120 h_enc=$HENC"
echo "requesting bits=$BITS, which built $BUILT bits on the reference machine;"
echo "quote the value the run prints, not this one, since the offset between the"
echo "request and the chain HElib builds is not portable across builds."
echo

exec env $ENVV /usr/bin/time -v "$BIN" \
     i="$IDX" bits="$BITS" mainhwt=120 h="$HENC" \
     t=-1 newbts=1 newks=1 thick=0 cache=1 repeat="$TRIALS"

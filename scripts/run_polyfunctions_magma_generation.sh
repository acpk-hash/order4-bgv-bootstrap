#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BASELINE="$ROOT/../baselines/Bootstrapping_Polyfunctions"
MAGMA_BIN="${MAGMA_BIN:-$ROOT/tools/bin/magma}"
WORKDIR="${POLYFUNCTIONS_MAGMA_WORKDIR:-$ROOT/results/polyfunctions_magma_generation}"
LOG="${POLYFUNCTIONS_MAGMA_LOG:-$ROOT/results/baseline_matrix/polyfunctions_magma_generation.log}"

if [[ ! -x "$MAGMA_BIN" ]]; then
  echo "Magma binary not found: $MAGMA_BIN" >&2
  echo "Run scripts/install_magma_local.sh first, or set MAGMA_BIN=/path/to/magma." >&2
  exit 2
fi

if [[ ! -d "$BASELINE/Scripts" ]]; then
  echo "Bootstrapping_Polyfunctions baseline not found: $BASELINE" >&2
  exit 2
fi

mkdir -p "$WORKDIR/HElib/HElib/src" "$WORKDIR/HElib/Polynomials/polynomials" "$(dirname "$LOG")"
cp -a "$BASELINE/Scripts" "$WORKDIR/"

RUNNER="$WORKDIR/run_find_digit_poly.m"
cat > "$RUNNER" <<'EOF'
SetSeed(1);
load "Scripts/Find_digit_poly.m";
quit;
EOF

(
  cd "$WORKDIR"
  "$MAGMA_BIN" "$RUNNER"
) 2>&1 | tee "$LOG"

echo "Generated files under: $WORKDIR/HElib"
echo "Log: $LOG"

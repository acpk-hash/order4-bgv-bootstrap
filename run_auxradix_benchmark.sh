#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="$ROOT/src/BGV-Boot-auxradix-opt/build/fatboot"
LOG_DIR="$ROOT/results"
mkdir -p "$LOG_DIR" "$ROOT/cache/saved_ZZX"

if [[ ! -x "$BIN" ]]; then
  echo "Missing $BIN. Run $ROOT/build_upload.sh first." >&2
  exit 1
fi

REPEAT="${REPEAT:-3}"
TIMEOUT="${TIMEOUT:-1200s}"
COMMON_ARGS=(i=4 h=12 t=-1 newbts=1 newks=1 thick=0 repeat="$REPEAT")

export HELIB_ZZX_CACHE_DIR="$ROOT/cache/saved_ZZX"

echo "Running isolated upload baseline aux=default repeat=$REPEAT"
timeout "$TIMEOUT" "$BIN" "${COMMON_ARGS[@]}" \
  > "$LOG_DIR/upload_aux_default_repeat${REPEAT}.log" 2>&1

echo "Running isolated upload optimized aux=256 repeat=$REPEAT"
HELIB_EXPLICIT_AUX=256 timeout "$TIMEOUT" "$BIN" "${COMMON_ARGS[@]}" \
  > "$LOG_DIR/upload_aux_256_repeat${REPEAT}.log" 2>&1

echo "Logs:"
echo "  $LOG_DIR/upload_aux_default_repeat${REPEAT}.log"
echo "  $LOG_DIR/upload_aux_256_repeat${REPEAT}.log"

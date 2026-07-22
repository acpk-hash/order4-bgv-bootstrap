#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HELIB_SRC="$ROOT/src/HElib_auxradix_opt"
HELIB_BUILD="$ROOT/build/HElib_auxradix_opt"
HELIB_INSTALL="$ROOT/local/helib_auxradix_opt"
FATBOOT_SRC="$ROOT/src/BGV-Boot-auxradix-opt"
FATBOOT_BUILD="$FATBOOT_SRC/build"

mkdir -p "$ROOT/cache/saved_ZZX" "$ROOT/results"

cmake -S "$HELIB_SRC" -B "$HELIB_BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$HELIB_INSTALL" \
  -DENABLE_TEST=OFF \
  -DENABLE_THREADS=ON \
  -DBUILD_SHARED=OFF

cmake --build "$HELIB_BUILD" --target install --parallel "${JOBS:-4}"

cmake -S "$FATBOOT_SRC" -B "$FATBOOT_BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -Dhelib_DIR="$HELIB_INSTALL/share/cmake/helib"

cmake --build "$FATBOOT_BUILD" --target fatboot --parallel "${JOBS:-4}"

echo "Built isolated upload implementation:"
echo "  HElib install: $HELIB_INSTALL"
echo "  fatboot:       $FATBOOT_BUILD/fatboot"

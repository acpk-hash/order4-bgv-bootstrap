#!/usr/bin/env bash
# Build every executable of the artifact from the sources in this repository.
# No network access is needed. Output goes to build/ and local/ (git-ignored);
# 2_baselines/ is never modified: each tree is copied to build/src/ first.
#
#   JOBS=8 3_experiments/build_all.sh              # everything (~10 min on 16 cores)
#   3_experiments/build_all.sh ours                # only some of: ours po2 baselines
#
# Executables in build/bin/:
#
#   ours (treatment and same-binary controls; OFF/ON selected at run time)
#     fatboot_order4       general rings and m = 50731   HElib = 1_implementation/helib_order4
#                                                          driver = 1_implementation/driver_order4
#     fatboot_po2          m = 2^16                       HElib 3e337a6 + Ma et al. NTT patch
#                                                          + 1_implementation/patches/helib-ntt-order4.patch;
#                                                          driver = Ma et al. NTT driver
#                                                          + 1_implementation/patches/ma-ntt-fatboot-order4.patch
#     enc_butterfly        Section 5 prototype            HElib as fatboot_order4
#
#   original baselines (control group G0; none of our code)
#     baseline_ma_largep   Ma et al. Large-p, unmodified  HElib 3e337a6 + bgvboot-final.patch
#     baseline_ma_ntt      Ma et al. NTT, unmodified      HElib 3e337a6 + github_20241003.patch
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
JOBS="${JOBS:-$(nproc)}"
IMPL="$ROOT/1_implementation"
BASE="$ROOT/2_baselines"
HELIB_UP="$BASE/HElib-3e337a6"
MA_LP="$BASE/BGV-Boot-for-Large-p-83b54d5"
MA_NTT="$BASE/bgv-bootstrapping-with-homomorphic-NTT-800bae0"
SRC="$ROOT/build/src"
BIN="$ROOT/build/bin"
CACHE="$ROOT/cache"
mkdir -p "$SRC" "$BIN" "$CACHE/saved_ZZX" "$CACHE/baseline_ma_largep" "$CACHE/baseline_ma_ntt"

log() { printf '\n=== %s\n' "$*"; }

fresh_copy() {  # fresh_copy <tree> <name>: a clean, writable copy in build/src
  local dst="$SRC/$2"
  [[ -e $dst ]] && chmod -R u+w "$dst" && find "$dst" -mindepth 1 -delete
  mkdir -p "$dst"
  cp -a "$1/." "$dst/"
  echo "$dst"
}

apply() {  # apply <dir> <patch>
  ( cd "$1" && git apply --whitespace=nowarn "$2" )
}

build_helib() {  # build_helib <src> <name>: install into local/<name>
  local src=$1 name=$2
  local bdir="$ROOT/build/$name" prefix="$ROOT/local/$name"
  cmake -S "$src" -B "$bdir" -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$prefix" -DENABLE_TEST=OFF \
        -DENABLE_THREADS=ON -DBUILD_SHARED=OFF > "$bdir.cmake.log"
  cmake --build "$bdir" --target install --parallel "$JOBS" > "$bdir.build.log"
}

build_driver() {  # build_driver <src> <helib name> <output name>
  local src=$1 helib="$ROOT/local/$2/share/cmake/helib" out=$3
  # upstream drivers hard-code helib_DIR; their README says to edit it
  sed -i "s|^set(helib_DIR .*)|set(helib_DIR $helib)|" "$src/CMakeLists.txt"
  cmake -S "$src" -B "$ROOT/build/$out" -DCMAKE_BUILD_TYPE=Release -Dhelib_DIR="$helib" > /dev/null
  cmake --build "$ROOT/build/$out" --target fatboot --parallel "$JOBS" > "$ROOT/build/$out.build.log"
  cp "$ROOT/build/$out/fatboot" "$BIN/$out"
}

build_ours() {
  log "ours: fatboot_order4 and enc_butterfly (HElib = 1_implementation/helib_order4)"
  build_helib "$IMPL/helib_order4" helib_order4
  local d; d=$(fresh_copy "$IMPL/driver_order4" driver_order4)
  build_driver "$d" helib_order4 fatboot_order4
  local L="$ROOT/local/helib_order4"
  "${CXX:-g++}" -O2 -std=c++17 -march=native -DHELIB_BOOT_THREADS -DHELIB_THREADS \
      "$IMPL/butterfly/enc_butterfly.cpp" -I"$L/include" -L"$L/lib" \
      -lhelib -lntl -lgmp -lpthread -o "$BIN/enc_butterfly"
}

build_po2() {
  log "ours: fatboot_po2 (HElib 3e337a6 + Ma et al. NTT patch + our port)"
  local h; h=$(fresh_copy "$HELIB_UP" helib_ntt_order4)
  apply "$h" "$MA_NTT/github_20241003.patch"
  apply "$h" "$IMPL/patches/helib-ntt-order4.patch"
  build_helib "$h" helib_ntt_order4
  local d; d=$(fresh_copy "$MA_NTT" driver_po2)
  apply "$d" "$IMPL/patches/ma-ntt-fatboot-order4.patch"
  build_driver "$d" helib_ntt_order4 fatboot_po2
}

build_baselines() {
  log "baseline: Ma et al. Large-p, unmodified (bgvboot-final.patch on HElib 3e337a6)"
  local h; h=$(fresh_copy "$HELIB_UP" helib_ma_largep)
  apply "$h" "$MA_LP/bgvboot-final.patch"
  sed -i "s|your_path_to_main_dir/saved_ZZX|$CACHE/baseline_ma_largep|" "$h/src/extractDigits.cpp"
  build_helib "$h" helib_ma_largep
  local d; d=$(fresh_copy "$MA_LP" baseline_ma_largep_src)
  build_driver "$d" helib_ma_largep baseline_ma_largep

  log "baseline: Ma et al. NTT, unmodified (github_20241003.patch on HElib 3e337a6)"
  h=$(fresh_copy "$HELIB_UP" helib_ma_ntt)
  apply "$h" "$MA_NTT/github_20241003.patch"
  # the patch leaves an unquoted placeholder that the upstream README says to replace
  sed -i "s|path_of_fatboot.cpp/saved_ZZX/%ld_%ld_ZZX.txt|\"$CACHE/baseline_ma_ntt/%ld_%ld_ZZX.txt\"|" \
      "$h/src/extractDigits.cpp"
  build_helib "$h" helib_ma_ntt
  d=$(fresh_copy "$MA_NTT" baseline_ma_ntt_src)
  build_driver "$d" helib_ma_ntt baseline_ma_ntt
}

targets=("$@"); [[ ${#targets[@]} -eq 0 ]] && targets=(ours po2 baselines)
for t in "${targets[@]}"; do "build_$t"; done
log "done"; ls -l "$BIN"

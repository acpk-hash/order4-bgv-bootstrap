#!/usr/bin/env bash
# Build the three bootstrapping executables used for the paper's ciphertext
# measurements. Everything is written below $ROOT/build, $ROOT/local and
# $ROOT/third_party/src (all git-ignored).
#
#   build/bin/fatboot_order4     ours, general cyclotomic rings and m = 50731
#                                (Tables 6 and 7). HElib = src/HElib_auxradix_opt,
#                                driver = src/BGV-Boot-auxradix-opt.
#   build/bin/fatboot_po2        ours, power-of-two ring m = 2^16 (Table 5).
#                                HElib 3e337a6 + Ma et al. homomorphic-NTT patch
#                                + third_party/patches/helib-ntt-order4.patch;
#                                driver = Ma et al. NTT fatboot + our patch.
#   build/bin/enc_butterfly      ours, Section 5 encrypted butterfly prototype
#                                (m = 50731), same HElib as fatboot_order4.
#   build/bin/fatboot_ma_largep  baseline of Ma et al. (Table 7, row "Ma et al."),
#                                unmodified upstream.
#
# Upstream code without an open-source licence (Ma et al. homomorphic-NTT) is
# fetched at a pinned commit, never redistributed; only our diff is shipped.
#
#   JOBS=8 ae/build_all.sh            # all three
#   ae/build_all.sh order4            # only one of: order4 po2 ma_largep
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
JOBS="${JOBS:-$(nproc)}"
TP="$ROOT/third_party"
SRC="$TP/src"
BIN="$ROOT/build/bin"
mkdir -p "$SRC" "$BIN" "$ROOT/cache/saved_ZZX"

HELIB_URL=https://github.com/homenc/HElib.git
HELIB_REV=3e337a66a91a92d49de6a9505340826b0eb71081
MA_NTT_URL=https://github.com/msh086/bgv-bootstrapping-with-homomorphic-NTT.git
MA_NTT_REV=800bae0
MA_LP_URL=https://github.com/msh086/BGV-Boot-for-Large-p.git
MA_LP_REV=83b54d5

log() { printf '\n=== %s\n' "$*"; }

fetch() {  # fetch <url> <rev> <dir>
  local url=$1 rev=$2 dir=$3
  local try
  for try in 1 2 3 4 5 6 7 8; do
    [[ -d $dir/.git ]] && break
    git -c http.version=HTTP/1.1 clone --quiet "$url" "$dir" && break
    rm -rf "$dir"; echo "clone of $url failed (attempt $try), retrying" >&2; sleep 10
  done
  git -C "$dir" checkout --quiet --force "$rev"
  git -C "$dir" clean -fdxq
}

build_helib() {  # build_helib <src> <prefix>
  local src=$1 prefix=$2
  local bdir="$ROOT/build/$(basename "$prefix")"
  cmake -S "$src" -B "$bdir" -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$prefix" -DENABLE_TEST=OFF \
        -DENABLE_THREADS=ON -DBUILD_SHARED=OFF > "$bdir.cmake.log"
  cmake --build "$bdir" --target install --parallel "$JOBS" > "$bdir.build.log"
}

build_order4() {
  log "ours (general rings, m=50731): HElib_auxradix_opt + BGV-Boot-auxradix-opt"
  build_helib "$ROOT/src/HElib_auxradix_opt" "$ROOT/local/helib_auxradix_opt"
  cmake -S "$ROOT/src/BGV-Boot-auxradix-opt" -B "$ROOT/build/fatboot_order4" \
        -DCMAKE_BUILD_TYPE=Release \
        -Dhelib_DIR="$ROOT/local/helib_auxradix_opt/share/cmake/helib" > /dev/null
  cmake --build "$ROOT/build/fatboot_order4" --target fatboot --parallel "$JOBS" > /dev/null
  cp "$ROOT/build/fatboot_order4/fatboot" "$BIN/fatboot_order4"
  # Section 5: encrypted butterfly prototype, same HElib
  local L="$ROOT/local/helib_auxradix_opt"
  "${CXX:-g++}" -O2 -std=c++17 -march=native -DHELIB_BOOT_THREADS -DHELIB_THREADS \
      "$ROOT/encrypted_linear_transform/enc_butterfly.cpp" -I"$L/include" -L"$L/lib" \
      -lhelib -lntl -lgmp -lpthread -o "$BIN/enc_butterfly"
}

build_po2() {
  log "ours (m=2^16): HElib $HELIB_REV + Ma NTT patch + order-four port"
  fetch "$MA_NTT_URL" "$MA_NTT_REV" "$SRC/ma_ntt"
  fetch "$HELIB_URL" "$HELIB_REV" "$SRC/helib_ntt"
  git -C "$SRC/helib_ntt" apply --whitespace=nowarn "$SRC/ma_ntt/github_20241003.patch"
  git -C "$SRC/helib_ntt" apply "$TP/patches/helib-ntt-order4.patch"
  build_helib "$SRC/helib_ntt" "$ROOT/local/helib_ntt_order4"
  git -C "$SRC/ma_ntt" apply "$TP/patches/ma-ntt-fatboot-order4.patch"
  sed -i "s|^set(helib_DIR .*)|set(helib_DIR $ROOT/local/helib_ntt_order4/share/cmake/helib)|" \
      "$SRC/ma_ntt/CMakeLists.txt"
  cmake -S "$SRC/ma_ntt" -B "$ROOT/build/fatboot_po2" -DCMAKE_BUILD_TYPE=Release > /dev/null
  cmake --build "$ROOT/build/fatboot_po2" --target fatboot --parallel "$JOBS" > /dev/null
  cp "$ROOT/build/fatboot_po2/fatboot" "$BIN/fatboot_po2"
}

build_ma_largep() {
  log "baseline Ma et al. (Large-p), unmodified upstream $MA_LP_REV"
  fetch "$MA_LP_URL" "$MA_LP_REV" "$SRC/ma_largep"
  fetch "$HELIB_URL" "$HELIB_REV" "$SRC/helib_largep"
  git -C "$SRC/helib_largep" apply --whitespace=nowarn "$SRC/ma_largep/bgvboot-final.patch"
  # upstream README: replace the placeholder by the polynomial-cache directory
  mkdir -p "$ROOT/cache/ma_largep_ZZX"
  sed -i "s|your_path_to_main_dir/saved_ZZX|$ROOT/cache/ma_largep_ZZX|" \
      "$SRC/helib_largep/src/extractDigits.cpp"
  build_helib "$SRC/helib_largep" "$ROOT/local/helib_ma_largep"
  sed -i "s|^set(helib_DIR .*)|set(helib_DIR $ROOT/local/helib_ma_largep/share/cmake/helib)|" \
      "$SRC/ma_largep/CMakeLists.txt"
  cmake -S "$SRC/ma_largep" -B "$ROOT/build/fatboot_ma_largep" -DCMAKE_BUILD_TYPE=Release > /dev/null
  cmake --build "$ROOT/build/fatboot_ma_largep" --target fatboot --parallel "$JOBS" > /dev/null
  cp "$ROOT/build/fatboot_ma_largep/fatboot" "$BIN/fatboot_ma_largep"
}

targets=("$@"); [[ ${#targets[@]} -eq 0 ]] && targets=(order4 po2 ma_largep)
for t in "${targets[@]}"; do "build_$t"; done
log "done"; ls -l "$BIN"

#!/usr/bin/env bash
# build.sh - Build the order-four BGV bootstrapping artifact
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HELIB_SRC="$ROOT/src/HElib_auxradix_opt"
HELIB_INSTALL="$ROOT/local/helib"
FATBOOT_SRC="$ROOT/src/BGV-Boot-auxradix-opt"

echo "=== Building Order-Four BGV Bootstrap Artifact ==="
echo "Root: $ROOT"

# Step 1: Build HElib with auxradix optimization
if [[ ! -f "$HELIB_INSTALL/lib/libhelib.a" ]]; then
    echo "[1/2] Building HElib (auxradix optimized)..."
    mkdir -p "$HELIB_SRC/build"
    cd "$HELIB_SRC/build"
    cmake .. \
        -DCMAKE_INSTALL_PREFIX="$HELIB_INSTALL" \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DBUILD_SHARED=OFF \
        -DENABLE_THREADS=ON \
        -DPACKAGE_BUILD=OFF \
        -DENABLE_TEST=OFF
    make -j$(nproc)
    make install
    echo "  HElib installed to $HELIB_INSTALL"
else
    echo "[1/2] HElib already built at $HELIB_INSTALL"
fi

# Step 2: Build fatboot driver
echo "[2/2] Building fatboot driver..."
mkdir -p "$FATBOOT_SRC/build"
cd "$FATBOOT_SRC/build"
cmake .. -DCMAKE_PREFIX_PATH="$HELIB_INSTALL"
make -j$(nproc)
echo "  fatboot binary: $FATBOOT_SRC/build/fatboot"

echo ""
echo "=== Build complete ==="
echo "Run benchmarks with: ./run_benchmark.sh"

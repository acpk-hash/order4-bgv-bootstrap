#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INSTALL_ROOT="${MAGMA_INSTALL_ROOT:-$ROOT/tools/magma}"
BIN_DIR="${MAGMA_BIN_DIR:-$ROOT/tools/bin}"
DOWNLOAD_DIR="${MAGMA_DOWNLOAD_DIR:-$ROOT/tools/magma_downloads}"

VERSION_PATH="${MAGMA_VERSION_PATH:-2/29/1}"
EXEC_NAME="${MAGMA_EXEC_NAME:-magma.avx2.exe.gz}"
EXEC_GZ="${MAGMA_EXEC_GZ:-$DOWNLOAD_DIR/$EXEC_NAME}"
SHARED_TGZ="${MAGMA_SHARED_TGZ:-$DOWNLOAD_DIR/shared_complete.tar.gz}"
PASSFILE="${MAGMA_PASSFILE:-$DOWNLOAD_DIR/magmapassfile}"

EXEC_URL="${MAGMA_EXEC_URL:-https://magma.maths.usyd.edu.au/magma/download/x86_64-linux/$VERSION_PATH/$EXEC_NAME}"
SHARED_URL="${MAGMA_SHARED_URL:-https://magma.maths.usyd.edu.au/magma/download/shared/shared_complete.tar.gz}"

mkdir -p "$INSTALL_ROOT" "$BIN_DIR" "$DOWNLOAD_DIR"

download_if_needed() {
  local url="$1"
  local out="$2"
  if [[ -f "$out" ]]; then
    return
  fi
  if [[ -z "${MAGMA_DOWNLOAD_USER:-}" || -z "${MAGMA_DOWNLOAD_PASS:-}" ]]; then
    echo "Missing $out" >&2
    echo "Put the official file there, or set MAGMA_DOWNLOAD_USER and MAGMA_DOWNLOAD_PASS." >&2
    exit 2
  fi
  curl -fL -u "$MAGMA_DOWNLOAD_USER:$MAGMA_DOWNLOAD_PASS" -o "$out" "$url"
}

download_if_needed "$EXEC_URL" "$EXEC_GZ"
download_if_needed "$SHARED_URL" "$SHARED_TGZ"

if [[ ! -f "$PASSFILE" ]]; then
  echo "Missing Magma passfile: $PASSFILE" >&2
  echo "Copy the host-bound magmapassfile from the Magma license portal to this path." >&2
  exit 2
fi

tmp_exec="$INSTALL_ROOT/magma.exe.tmp"
gzip -dc "$EXEC_GZ" > "$tmp_exec"
chmod 755 "$tmp_exec"
mv "$tmp_exec" "$INSTALL_ROOT/magma.exe"

tar -xzf "$SHARED_TGZ" -C "$INSTALL_ROOT"
cp "$PASSFILE" "$INSTALL_ROOT/magmapassfile"
chmod 600 "$INSTALL_ROOT/magmapassfile"

cat > "$BIN_DIR/magma" <<EOF
#!/usr/bin/env bash
export MAGMA_HOME="$INSTALL_ROOT"
exec "$INSTALL_ROOT/magma.exe" "\$@"
EOF
chmod 755 "$BIN_DIR/magma"

echo "Installed Magma wrapper: $BIN_DIR/magma"
"$BIN_DIR/magma" -v || true

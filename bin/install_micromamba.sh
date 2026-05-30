#!/usr/bin/env bash
# Install micromamba and expose it as the `mamba` frontend Snakemake expects.
# micromamba's solver is far faster than classic conda, which makes building the
# ~10 per-rule environments quick and reliable. After running this, build/run with
# `--conda-frontend mamba` (see README).
#
#   bash bin/install_micromamba.sh
#   export PATH="$HOME/.local/bin:$PATH"   # add to your shell rc
#
# Notes:
# * Pinned to micromamba 1.5.8: micromamba 2.x aborts on Snakemake's pre-created
#   --prefix dir ("Non-conda folder exists at prefix"); 1.5.x works with the wrapper below.
# * The `mamba` wrapper removes any incomplete prefix (a dir lacking conda-meta) that
#   Snakemake pre-creates, so micromamba will populate it.
set -euo pipefail

VERSION="${MICROMAMBA_VERSION:-1.5.8}"
BIN="$HOME/.local/bin"
mkdir -p "$BIN"

echo ">> downloading micromamba ${VERSION}"
curl -Ls "https://micro.mamba.pm/api/micromamba/linux-64/${VERSION}" | tar -xj -C "$HOME/.local" bin/micromamba
chmod +x "$BIN/micromamba"

echo ">> writing mamba shim -> micromamba (Snakemake-compatible)"
cat > "$BIN/mamba" <<'EOF'
#!/usr/bin/env bash
real="$HOME/.local/bin/micromamba"
prev=""
for a in "$@"; do
  if [[ "$prev" == "--prefix" || "$prev" == "-p" ]]; then
    [[ -d "$a" && ! -d "$a/conda-meta" ]] && rm -rf "$a"
  fi
  prev="$a"
done
exec "$real" "$@"
EOF
chmod +x "$BIN/mamba"

echo ">> done. micromamba $($BIN/micromamba --version), mamba shim $($BIN/mamba --version)"
echo ">> ensure PATH includes $BIN, then run snakemake with: --conda-frontend mamba"

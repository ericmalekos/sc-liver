#!/usr/bin/env bash
# Build the runtime images and convert to Singularity .sif for use on SLURM.
# Run from the repo root: bash containers/build.sh [base|gpu|all]
set -euo pipefail
cd "$(dirname "$0")/.."

REG="${REG:-ghcr.io/ericmalekos/sc-liver}"
TAG="${TAG:-0.1.0}"
WHICH="${1:-base}"

build_one () {
    local name="$1" dockerfile="$2"
    echo ">> docker build ${REG}-${name}:${TAG}"
    docker build -f "containers/${dockerfile}" -t "${REG}-${name}:${TAG}" .
    echo ">> singularity build containers/${name}.sif"
    singularity build --force "containers/${name}.sif" "docker-daemon://${REG}-${name}:${TAG}"
}

case "$WHICH" in
    base) build_one base Dockerfile.base ;;
    gpu)  build_one gpu  Dockerfile.gpu ;;
    all)  build_one base Dockerfile.base; build_one gpu Dockerfile.gpu ;;
    *) echo "usage: build.sh [base|gpu|all]"; exit 1 ;;
esac
echo "Done. Point rules at containers/<name>.sif and run with --use-singularity."

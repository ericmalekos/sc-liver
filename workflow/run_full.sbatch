#!/bin/bash
# Submit the full real-data pipeline as a single detached SLURM job (snakemake runs
# locally within one large allocation -> robust to login-node memory limits and session
# disconnects). For per-rule cluster submission instead, use `--profile profiles/slurm`.
#   sbatch workflow/run_full.sbatch
#SBATCH --job-name=sc-fibrosis
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm/sc-fibrosis.%j.out
#SBATCH --error=logs/slurm/sc-fibrosis.%j.err
set -euo pipefail
# run from the directory sbatch was submitted from (the repo root); works regardless of path
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
mkdir -p logs/slurm
source /private/home/emalekos/miniconda3/etc/profile.d/conda.sh
conda activate sc-liver
export PATH="$HOME/.local/bin:$PATH"
export MAMBA_ROOT_PREFIX="$HOME/.local/micromamba"
# CRITICAL: ignore ~/.local user-site packages so each rule uses its conda env's PINNED
# versions (otherwise the user's numpy 2.4 / anndata 0.12 shadow the envs and crash e.g.
# scrublet). Inherited by all snakemake-spawned rule jobs.
export PYTHONNOUSERSITE=1
# keep caches/tmp on the large project filesystem (HOME has a tight quota)
export XDG_CACHE_HOME="$PWD/.cache"
export TMPDIR="$PWD/.tmp"
mkdir -p "$XDG_CACHE_HOME" "$TMPDIR"

snakemake --unlock --configfile config/config.yaml || true
snakemake --use-conda --conda-frontend mamba --configfile config/config.yaml \
  --cores 12 --keep-going --rerun-incomplete --restart-times 1
echo "=== pipeline finished; key outputs ==="
ls -la results/08_score/candidate_scores.tsv results/10_report/report.html 2>/dev/null || true

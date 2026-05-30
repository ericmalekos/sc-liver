#!/usr/bin/env bash
# Robust SLURM job-status script for Snakemake's `cluster-status`.
# Snakemake passes the job id (the value printed by `sbatch --parsable`).
# Must print exactly one of: "success", "failed", "running".
set -euo pipefail
jobid="$1"

# Prefer sacct (survives after the job leaves the queue); fall back to scontrol.
status=$(sacct -j "$jobid" --format=State --noheader --parsable2 2>/dev/null | head -n1 || true)
if [[ -z "$status" ]]; then
    status=$(scontrol show job "$jobid" 2>/dev/null | grep -oP 'JobState=\K\S+' || true)
fi

case "$status" in
    COMPLETED)                                  echo "success" ;;
    RUNNING|PENDING|REQUEUED|RESIZING|SUSPENDED|COMPLETING|CONFIGURING) echo "running" ;;
    BOOT_FAIL|CANCELLED*|DEADLINE|FAILED|NODE_FAIL|OUT_OF_MEMORY|PREEMPTED|TIMEOUT) echo "failed" ;;
    "")                                         echo "running" ;;   # not yet visible to sacct
    *)                                          echo "running" ;;
esac

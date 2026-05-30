"""Provenance — capture the resolved config, git SHA, and tool versions for the run."""
import json
import os
import subprocess
import sys
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import ensure_parent, get_logger  # noqa: E402

import yaml  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
config = dict(snakemake.params.config)  # noqa: F821


def sh(cmd):
    try:
        return subprocess.check_output(cmd, shell=True, text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return None


meta = {
    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    "git_sha": sh("git rev-parse HEAD"),
    "git_dirty": bool(sh("git status --porcelain")),
    "python": sys.version.split()[0],
    "snakemake": sh("snakemake --version"),
    "conda": sh("conda --version"),
    "host": sh("hostname"),
    "project": config.get("project", {}),
    "datasets": {k: v.get("geo_accession") for k, v in config.get("datasets", {}).items()},
    "integration_method": config.get("integration", {}).get("method"),
    "de_engine": config.get("de", {}).get("engine"),
    "score_weights": config.get("score", {}).get("weights"),
}
ensure_parent(snakemake.output.meta)  # noqa: F821
with open(snakemake.output.meta, "w") as fh:  # noqa: F821
    json.dump(meta, fh, indent=2)
with open(snakemake.output.cfg, "w") as fh:  # noqa: F821
    yaml.safe_dump(config, fh, sort_keys=False)
log.info(f"provenance written (git {meta['git_sha']}, snakemake {meta['snakemake']})")

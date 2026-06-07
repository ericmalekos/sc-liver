"""Shared helpers for the pipeline's Python rule-scripts.

Used via Snakemake's `script:` directive, which injects a global `snakemake`
object. Keep this import-safe (no side effects at import time).
"""

from __future__ import annotations

import logging
import os
import random
from pathlib import Path

import numpy as np


def get_logger(snakemake=None, name: str = "sc") -> logging.Logger:
    """Logger that writes to the rule's log file (if any) and stderr."""
    handlers = [logging.StreamHandler()]
    if snakemake is not None and getattr(snakemake, "log", None):
        logpath = str(snakemake.log[0])
        Path(logpath).parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(logpath))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=handlers,
        force=True,
    )
    rule = getattr(getattr(snakemake, "rule", None), "__str__", lambda: name)()
    return logging.getLogger(rule or name)


def set_all_seeds(seed: int) -> None:
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch

        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except Exception:
        pass


def ensure_parent(path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def read_h5ad(path):
    import anndata as ad

    return ad.read_h5ad(str(path))


def write_h5ad(adata, path) -> None:
    ensure_parent(path)
    adata.write_h5ad(str(path))


def use_agg():
    import matplotlib

    matplotlib.use("Agg")

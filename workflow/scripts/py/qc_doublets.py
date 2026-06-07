"""Stage 02 — doublet scoring with Scrublet (Python path).

Robust to tiny inputs: if Scrublet can't fit a threshold (small fixture data), it falls
back to zero scores so the QC filter simply keeps all cells.
"""

import os

# Pin BLAS/OpenMP to single-threaded BEFORE numpy is imported. On the 256-core compute
# nodes, OpenBLAS otherwise spawns ~256 threads per process; running many Scrublet jobs in
# parallel then explodes the thread count and segfaults (the original empty-log crashes).
for _v in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMBA_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_v, "1")

import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import get_logger, read_h5ad, write_h5ad  # noqa: E402

import numpy as np  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
seed = int(snakemake.params.seed)  # noqa: F821
rate = float(snakemake.params.expected_rate)  # noqa: F821
method = snakemake.params.method  # noqa: F821
log.info(f"doublet detection: method={method}")
adata = read_h5ad(snakemake.input.h5ad)  # noqa: F821

scores = np.zeros(adata.n_obs, dtype=float)
predicted = np.zeros(adata.n_obs, dtype=bool)
if method == "none":
    log.info("doublet detection disabled (already-filtered/processed inputs)")
    adata.obs["doublet_score"] = scores
    adata.obs["predicted_doublet"] = predicted
    write_h5ad(adata, snakemake.output.h5ad)  # noqa: F821
    sys.exit(0)
try:
    import scrublet as scr

    counts = adata.layers["counts"]
    scrub = scr.Scrublet(counts, expected_doublet_rate=rate, random_state=seed)
    # use_approx_neighbors=False -> sklearn EXACT kNN instead of annoy, whose pip wheel
    # raises SIGILL (illegal instruction) on the compute nodes. Fast enough at this scale.
    scores, pred = scrub.scrub_doublets(
        min_counts=2, min_cells=3, use_approx_neighbors=False, verbose=False
    )
    if pred is not None:
        predicted = np.asarray(pred, dtype=bool)
    log.info(f"Scrublet flagged {int(predicted.sum())}/{adata.n_obs} doublets")
except Exception as e:  # tiny / degenerate data
    log.warning(f"Scrublet fallback (kept all cells): {e}")

adata.obs["doublet_score"] = scores
adata.obs["predicted_doublet"] = predicted
write_h5ad(adata, snakemake.output.h5ad)  # noqa: F821

"""Stage 02 — decontX ambient-RNA removal via a lightweight Python<->R file hand-off.

Avoids the zellkonverter/basilisk bridge (which bootstraps a huge nested env and stalls):
Python reads the .h5ad and writes counts as a MatrixMarket file; a small R script
(decontx_run.R) runs decontX and writes corrected counts back; Python loads them and writes
the .ambient.h5ad. Falls back to passthrough on any failure so the pipeline never stalls.
"""
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import get_logger, read_h5ad, write_h5ad  # noqa: E402

import numpy as np  # noqa: E402
import scipy.io  # noqa: E402
import scipy.sparse as sp  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
method = snakemake.params.method  # noqa: F821
adata = read_h5ad(snakemake.input.h5ad)  # noqa: F821
log.info(f"ambient removal: {method} (file hand-off; no basilisk)")

counts = adata.layers.get("counts", adata.X)
genes_x_cells = sp.csr_matrix(counts).T.tocsc()  # decontX wants genes x cells

ok = False
with tempfile.TemporaryDirectory() as tmp:
    mtx_in = os.path.join(tmp, "counts.mtx")
    out_prefix = os.path.join(tmp, "out")
    scipy.io.mmwrite(mtx_in, genes_x_cells, field="integer")
    rscript = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "R", "decontx_run.R")
    res = subprocess.run(["Rscript", rscript, mtx_in, out_prefix], capture_output=True, text=True)
    log.info((res.stdout or "")[-1500:])
    dec_path = out_prefix + "_decontaminated.mtx"
    if res.returncode == 0 and os.path.exists(dec_path):
        dec = sp.csr_matrix(scipy.io.mmread(dec_path).T)  # back to cells x genes
        cont = np.loadtxt(out_prefix + "_contamination.txt")
        adata.X = dec
        adata.layers["counts"] = dec.copy()
        adata.obs["ambient_fraction"] = np.asarray(cont, dtype=float).ravel()
        ok = True
    else:
        log.warning(f"decontX failed; passing counts through. stderr: {(res.stderr or '')[-800:]}")
        adata.obs["ambient_fraction"] = 0.0

adata.uns["ambient_method"] = method if ok else f"{method}_failed_passthrough"
write_h5ad(adata, snakemake.output.h5ad)  # noqa: F821
log.info(f"ambient done (decontX ok={ok}): {adata.n_obs} cells")

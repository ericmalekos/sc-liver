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


def quick_clusters(counts_cells_x_genes, n_top=2000, n_pcs=30, k=15, seed=0):
    """Coarse cell clusters (scipy/numpy only -- this env has no scanpy) to pass to decontX
    as priors. With cluster priors decontX removes CROSS-population ambient without stripping
    a population's genuine marker genes (e.g. SPP1 in scar macrophages); its default self-
    clustering otherwise over-corrects high-expressed markers."""
    from scipy.cluster.vq import kmeans2
    X = sp.csr_matrix(counts_cells_x_genes).astype(float)          # cells x genes
    if X.shape[0] < 6:
        return np.ones(X.shape[0], dtype=int)
    lib = np.asarray(X.sum(1)).ravel(); lib[lib == 0] = 1.0
    Xn = X.multiply(1e4 / lib[:, None]).tocsr(); Xn.data = np.log1p(Xn.data)
    mean = np.asarray(Xn.mean(0)).ravel()
    var = np.asarray(Xn.multiply(Xn).mean(0)).ravel() - mean ** 2
    top = np.argsort(var)[::-1][:min(n_top, Xn.shape[1])]
    D = np.asarray(Xn[:, top].todense()); D -= D.mean(0)
    npc = int(min(n_pcs, D.shape[1] - 1, D.shape[0] - 1))
    U, S, _ = np.linalg.svd(D, full_matrices=False)
    scores = U[:, :npc] * S[:npc]
    kk = int(max(2, min(k, D.shape[0] // 50)))
    _, labels = kmeans2(scores, kk, seed=seed, minit="++", missing="warn")
    _, labels = np.unique(labels, return_inverse=True)   # contiguous
    return labels.astype(int) + 1                        # 1-based for R

genes_x_cells = sp.csr_matrix(counts).T.tocsc()  # decontX wants genes x cells

ok = False
with tempfile.TemporaryDirectory() as tmp:
    mtx_in = os.path.join(tmp, "counts.mtx")
    out_prefix = os.path.join(tmp, "out")
    scipy.io.mmwrite(mtx_in, genes_x_cells, field="integer")
    # cluster priors so decontX protects genuine cell-type markers
    z_file = os.path.join(tmp, "z.txt")
    try:
        z = quick_clusters(counts)
        np.savetxt(z_file, z, fmt="%d")
        log.info(f"decontX cluster priors: {len(np.unique(z))} clusters over {len(z)} cells")
    except Exception as e:
        z_file = ""
        log.warning(f"cluster-prior computation failed ({e}); decontX will self-cluster")
    rscript = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "R", "decontx_run.R")
    cmd = ["Rscript", rscript, mtx_in, out_prefix] + ([z_file] if z_file else [])
    res = subprocess.run(cmd, capture_output=True, text=True)
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

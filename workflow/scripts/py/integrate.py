"""Stage 03 — concatenate QC'd samples, normalize, find HVGs, and batch-correct
OVER DONOR/SAMPLE (never over disease condition) so fibrosis biology is preserved (Q3).

Default method is Harmony (CPU). scVI/scANVI is used only when gpu.enabled. Full gene
set is retained (counts layer kept) so downstream pseudobulk DE uses all genes; PCA is
computed on HVGs only.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import get_logger, set_all_seeds, use_agg, write_h5ad  # noqa: E402

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import scanpy as sc  # noqa: E402

use_agg()
import matplotlib.pyplot as plt  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
p = snakemake.params  # noqa: F821
set_all_seeds(int(p.seed))
batch_key = p.batch_key
n_hvg = int(p.n_hvg)

adatas = [sc.read_h5ad(f) for f in snakemake.input.h5ads]  # noqa: F821
# OUTER join (union of genes), NOT inner: with FACS-sorted samples (CD45+/CD45-), an inner
# join keeps only genes detected in EVERY sample, which silently deletes the cell-type-
# specific markers that fibrosis biology lives in (TREM2/SPP1 in macrophages, COL3A1/PDGFRB
# in mesenchyme) before pseudobulk DE ever sees them. Missing genes fill with 0 (true zeros
# for cells of the wrong sorted fraction); DESeq2's own per-compartment filter handles the
# rest. PCA/integration still run on HVGs only, so the embedding is unchanged.
adata = ad.concat(adatas, join="outer", fill_value=0, index_unique=None, merge="same")
adata.obs_names_make_unique()
adata.layers["counts"] = adata.layers.get("counts", adata.X).copy()
log.info(f"Concatenated {len(adatas)} samples -> {adata.n_obs} cells x {adata.n_vars} genes; "
         f"batches({batch_key})={adata.obs[batch_key].nunique()}")

# normalize / log1p on the full gene set
adata.X = adata.layers["counts"].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["lognorm"] = adata.X.copy()

n_hvg = min(n_hvg, adata.n_vars)
sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat_v3",
                            layer="counts", batch_key=batch_key, subset=False)

n_pcs = int(min(int(p.n_latent), adata.var["highly_variable"].sum() - 1, adata.n_obs - 1))
n_pcs = max(n_pcs, 2)
sc.pp.pca(adata, n_comps=n_pcs, use_highly_variable=True, random_state=int(p.seed))

method = p.method
gpu = bool(p.gpu)
emb = "X_emb"
if method in ("scvi", "scanvi") and gpu:
    import scvi

    sub = adata[:, adata.var["highly_variable"]].copy()
    sub.layers["counts"] = sub.layers["counts"].copy()
    scvi.settings.seed = int(p.seed)
    scvi.model.SCVI.setup_anndata(sub, layer="counts", batch_key=batch_key)
    model = scvi.model.SCVI(sub, n_latent=int(p.n_latent))
    model.train()
    if method == "scanvi" and "cell_type" in adata.obs:
        model = scvi.model.SCANVI.from_scvi_model(model, unlabeled_category="Unknown",
                                                  labels_key="cell_type")
        model.train()
    adata.obsm[emb] = model.get_latent_representation()
    log.info(f"{method} latent: {adata.obsm[emb].shape}")
else:
    if method in ("scvi", "scanvi"):
        log.warning(f"{method} requested but gpu.enabled is false -> using Harmony (CPU).")
    sc.external.pp.harmony_integrate(adata, key=batch_key, basis="X_pca",
                                     adjusted_basis="X_pca_harmony", random_state=int(p.seed))
    adata.obsm[emb] = adata.obsm["X_pca_harmony"]
    method = "harmony"

adata.uns["integration"] = {"method": method, "batch_key": batch_key, "n_hvg": int(n_hvg)}
sc.pp.neighbors(adata, use_rep=emb, n_neighbors=int(p.n_neighbors), random_state=int(p.seed))
sc.tl.umap(adata, random_state=int(p.seed))

# overview UMAP: batch vs condition (visual check that batches mix but condition persists)
color = [c for c in [batch_key, "condition", "fibrosis_axis"] if c in adata.obs]
fig = sc.pl.umap(adata, color=color, ncols=3, show=False, return_fig=True)
os.makedirs(os.path.dirname(snakemake.output.umap), exist_ok=True)  # noqa: F821
fig.savefig(snakemake.output.umap, dpi=110, bbox_inches="tight")  # noqa: F821
plt.close(fig)

write_h5ad(adata, snakemake.output.h5ad)  # noqa: F821
log.info(f"Integrated with {method}; wrote {adata.n_obs} cells")

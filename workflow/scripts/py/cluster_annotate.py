"""Stage 04a — Leiden clustering + marker-based cluster labelling (Q4).

Scores each curated cell-type signature per cell, assigns every Leiden cluster the
best-scoring cell type, and writes a marker dotplot. Final compartment assignment and
validation happen in celltypist_annotate.py.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import get_logger, set_all_seeds, use_agg, write_h5ad  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import yaml  # noqa: E402

use_agg()
import matplotlib.pyplot as plt  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
seed = int(snakemake.params.seed)  # noqa: F821
res = float(snakemake.params.resolution)  # noqa: F821
set_all_seeds(seed)

adata = sc.read_h5ad(snakemake.input.h5ad)  # noqa: F821
with open(snakemake.input.markers) as fh:  # noqa: F821
    markers = yaml.safe_load(fh)
adata.X = adata.layers.get("lognorm", adata.X)

sc.tl.leiden(adata, resolution=res, random_state=seed, key_added="leiden")
log.info(f"Leiden: {adata.obs['leiden'].nunique()} clusters")

# score each cell-type signature
ct_markers = markers["cell_types"]
score_cols = {}
for ct, genes in ct_markers.items():
    present = [g for g in genes if g in adata.var_names]
    col = f"score_{ct}"
    if present:
        sc.tl.score_genes(adata, present, score_name=col, random_state=seed)
    else:
        adata.obs[col] = 0.0
    score_cols[ct] = col

# assign each cluster the argmax mean signature
per_cluster = adata.obs.groupby("leiden")[list(score_cols.values())].mean()
cluster_label = {}
for cl, row in per_cluster.iterrows():
    best = row.idxmax()
    cluster_label[cl] = best.replace("score_", "")
adata.obs["cell_type_marker"] = adata.obs["leiden"].map(cluster_label).astype("category")
log.info(f"Cluster->cell_type: {cluster_label}")

# marker dotplot (a few canonical markers per major type)
panel = {}
for ct, genes in ct_markers.items():
    present = [g for g in genes if g in adata.var_names][:4]
    if present:
        panel[ct] = present
try:
    fig = sc.pl.dotplot(adata, panel, groupby="leiden", show=False, return_fig=True)
    os.makedirs(os.path.dirname(snakemake.output.dotplot), exist_ok=True)  # noqa: F821
    fig.savefig(snakemake.output.dotplot, dpi=110, bbox_inches="tight")  # noqa: F821
    plt.close("all")
except Exception as e:
    log.warning(f"dotplot failed, writing placeholder: {e}")
    plt.figure(); plt.text(0.5, 0.5, "dotplot unavailable", ha="center")
    plt.savefig(snakemake.output.dotplot, dpi=90)  # noqa: F821

write_h5ad(adata, snakemake.output.h5ad)  # noqa: F821

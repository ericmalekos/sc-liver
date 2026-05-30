"""Stage 04b — finalize cell types, roll up to compartments, and validate (Q4).

* Final cell type per cluster from the marker labels (+ optional CellTypist cross-check).
* Maps cell types to the analysis compartments (stellate/myofibroblast, macrophage/monocyte,
  endothelial, ...).
* Builds compartment_validation.tsv: per cluster, the cell type / compartment, the
  cell-type signature scores, the disease-state scores (scar-macrophage, scar-endothelium,
  activated-HSC), AND the stromal discriminators (HSC vs portal-fibroblast vs myofibroblast
  vs VSMC) — the explicit answer to the Q4 disambiguation.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import ensure_parent, get_logger, use_agg, write_h5ad  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import yaml  # noqa: E402

use_agg()
import matplotlib.pyplot as plt  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
model_name = snakemake.params.celltypist_model  # noqa: F821
adata = sc.read_h5ad(snakemake.input.h5ad)  # noqa: F821
with open(snakemake.input.markers) as fh:  # noqa: F821
    markers = yaml.safe_load(fh)
adata.X = adata.layers.get("lognorm", adata.X)

# cell_type -> compartment reverse map
ct2comp = {}
for comp, cts in markers["compartments"].items():
    for ct in cts:
        ct2comp[ct] = comp

adata.obs["cell_type"] = adata.obs["cell_type_marker"].astype(str)

# optional CellTypist cross-check (majority vote over clusters)
if model_name:
    try:
        import celltypist
        from celltypist import models

        model = models.Model.load(model=model_name)
        pred = celltypist.annotate(adata, model=model, majority_voting=True)
        adata.obs["celltypist"] = pred.predicted_labels["majority_voting"].values
        log.info("CellTypist annotation added as cross-check column 'celltypist'.")
    except Exception as e:
        log.warning(f"CellTypist skipped: {e}")

adata.obs["compartment"] = (
    adata.obs["cell_type"].map(ct2comp).fillna("myeloid_other").astype("category")
)
log.info(f"Compartments: {adata.obs['compartment'].value_counts().to_dict()}")

# score disease states + stromal discriminators for the validation table
extra = {}
for grp in ("disease_states", "stromal_discriminators"):
    for name, genes in markers[grp].items():
        present = [g for g in genes if g in adata.var_names]
        col = f"state_{name}"
        if present:
            sc.tl.score_genes(adata, present, score_name=col)
        else:
            adata.obs[col] = 0.0
        extra[name] = col

# per-cluster validation table
score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
state_cols = [c for c in adata.obs.columns if c.startswith("state_")]
g = adata.obs.groupby("leiden")
val = g[score_cols + state_cols].mean()
val.insert(0, "n_cells", g.size())
val.insert(1, "cell_type", g["cell_type"].agg(lambda s: s.mode().iat[0]))
val.insert(2, "compartment", g["compartment"].agg(lambda s: s.mode().iat[0]))
# fibrosis enrichment per cluster (fraction from high-fibrosis samples), supports Q4
if "fibrosis_axis" in adata.obs:
    val["frac_fibrotic_F3plus"] = g["fibrosis_axis"].agg(lambda s: float((s.astype(int) >= 3).mean()))
val = val.round(3)
ensure_parent(snakemake.output.validation)  # noqa: F821
val.to_csv(snakemake.output.validation, sep="\t")  # noqa: F821

# UMAP coloured by cell type / compartment / a couple disease markers
color = [c for c in ["cell_type", "compartment", "fibrosis_axis"] if c in adata.obs]
for g_ in ["COL1A1", "TREM2", "PLVAP"]:
    if g_ in adata.var_names:
        color.append(g_)
fig = sc.pl.umap(adata, color=color, ncols=3, show=False, return_fig=True)
os.makedirs(os.path.dirname(snakemake.output.umap), exist_ok=True)  # noqa: F821
fig.savefig(snakemake.output.umap, dpi=110, bbox_inches="tight")  # noqa: F821
plt.close("all")

write_h5ad(adata, snakemake.output.h5ad)  # noqa: F821
log.info(f"Annotated {adata.n_obs} cells; validation table -> {snakemake.output.validation}")  # noqa: F821

"""Stage 02 — adaptive MAD filtering (screening Q2).

Removes poor-quality cells and doublets WITHOUT discarding biologically meaningful
stressed/diseased cells:
  * adaptive n-MAD bounds on log1p(total_counts) and log1p(n_genes) (Heumos 2023);
  * a liver/modality-aware upper bound on mitochondrial % (hepatocytes are mito-high;
    snRNA is mito-low) — only an UPPER bound, never a lower one;
  * drops Scrublet-flagged doublets;
  * does NOT filter on stress/IEG signatures (keep_stressed) so activated HSCs,
    scar-associated macrophages, etc. survive.
Writes a metrics JSON and a before/after figure.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import ensure_parent, get_logger, read_h5ad, use_agg, write_h5ad  # noqa: E402

import numpy as np  # noqa: E402
import scanpy as sc  # noqa: E402

use_agg()
import matplotlib.pyplot as plt  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
p = snakemake.params  # noqa: F821
modality = p.modality
mito_max = float(p.mito_pct_max[modality]) if isinstance(p.mito_pct_max, dict) else float(p.mito_pct_max)
n_mads = float(p.n_mads)

adata = read_h5ad(snakemake.input.h5ad)  # noqa: F821
adata.layers.setdefault("counts", adata.X.copy())
if "total_counts" not in adata.obs:
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(("MT-", "MT."))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
n_before = adata.n_obs


def mad_bounds(x, nmads, upper_only=False):
    med = np.median(x)
    mad = np.median(np.abs(x - med)) or 1e-9
    lo = -np.inf if upper_only else med - nmads * mad
    return lo, med + nmads * mad


lc_lo, lc_hi = mad_bounds(np.log1p(adata.obs["total_counts"].values), n_mads)
lg_lo, lg_hi = mad_bounds(np.log1p(adata.obs["n_genes_by_counts"].values), n_mads)
_, mt_hi = mad_bounds(adata.obs["pct_counts_mt"].values, n_mads, upper_only=True)
mt_ceiling = min(mito_max, mt_hi) if np.isfinite(mt_hi) else mito_max

logc = np.log1p(adata.obs["total_counts"].values)
logg = np.log1p(adata.obs["n_genes_by_counts"].values)
keep = (
    (logc >= lc_lo) & (logc <= lc_hi)
    & (logg >= lg_lo) & (logg <= lg_hi)
    & (adata.obs["pct_counts_mt"].values <= mt_ceiling)
    & (adata.obs["total_counts"].values >= float(p.min_counts))
    & (adata.obs["n_genes_by_counts"].values >= float(p.min_genes))
)
n_mad_removed = int((~keep).sum())
dbl = adata.obs.get("predicted_doublet")
if dbl is not None:
    keep = keep & ~np.asarray(dbl, dtype=bool)

adata_f = adata[keep].copy()
sc.pp.filter_genes(adata_f, min_cells=int(p.min_cells))
n_after = adata_f.n_obs

metrics = dict(
    sample=str(adata.obs["sample_id"].iloc[0]) if "sample_id" in adata.obs else "NA",
    dataset=str(adata.obs["dataset"].iloc[0]) if "dataset" in adata.obs else "NA",
    modality=modality,
    n_cells_before=int(n_before),
    n_cells_after=int(n_after),
    n_genes_after=int(adata_f.n_vars),
    pct_kept=round(100 * n_after / max(n_before, 1), 1),
    doublets_removed=int(np.asarray(dbl, dtype=bool).sum()) if dbl is not None else 0,
    thresholds=dict(
        log1p_counts=[round(float(lc_lo), 3), round(float(lc_hi), 3)],
        log1p_genes=[round(float(lg_lo), 3), round(float(lg_hi), 3)],
        mito_pct_ceiling=round(float(mt_ceiling), 2),
        n_mads=n_mads,
        keep_stressed=bool(p.keep_stressed),
    ),
)
ensure_parent(snakemake.output.metrics)  # noqa: F821
with open(snakemake.output.metrics, "w") as fh:  # noqa: F821
    json.dump(metrics, fh, indent=2)

# before/after QC figure
ensure_parent(snakemake.output.fig)  # noqa: F821
fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))
for ax, col, title in zip(
    axes,
    ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
    ["total counts", "n genes", "% mito"],
):
    ax.hist(adata.obs[col], bins=40, alpha=0.5, label="before")
    ax.hist(adata_f.obs[col], bins=40, alpha=0.7, label="after")
    ax.set_title(title)
    ax.set_yscale("log")
axes[0].legend()
fig.suptitle(f"{metrics['sample']} QC ({n_before}->{n_after} cells)")
fig.tight_layout()
fig.savefig(snakemake.output.fig, dpi=110)  # noqa: F821

write_h5ad(adata_f, snakemake.output.h5ad)  # noqa: F821
log.info(f"QC {metrics['sample']}: {n_before}->{n_after} cells "
         f"(MAD removed {n_mad_removed}, mito ceiling {mt_ceiling:.1f}%)")

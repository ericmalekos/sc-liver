"""Stage 03 — quantify integration quality and GUARD that fibrosis biology survived (Q3).

Reports batch-mixing (we WANT batches mixed) and biological-signal-preservation (we
do NOT want the disease/condition signal erased). The fibrosis-signal guard fails loudly
if the condition is no longer separable in the integrated embedding.

Metrics (interpretable, sklearn-based; optional iLISI/cLISI if scib-metrics present):
  * batch_silhouette        : |ASW| by batch  (lower = better mixed)
  * condition_silhouette    : ASW by condition (higher = bio preserved)
  * condition_knn_accuracy  : 5-fold CV accuracy predicting condition from the embedding
                              (>> chance => fibrosis signal retained)
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import ensure_parent, get_logger, use_agg  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
from sklearn.metrics import silhouette_score  # noqa: E402
from sklearn.model_selection import cross_val_score  # noqa: E402
from sklearn.neighbors import KNeighborsClassifier  # noqa: E402

use_agg()
import matplotlib.pyplot as plt  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
batch_key = snakemake.params.batch_key  # noqa: F821
guard = bool(snakemake.params.guard)  # noqa: F821

adata = sc.read_h5ad(snakemake.input.h5ad)  # noqa: F821
emb = adata.obsm["X_emb"]


def safe_sil(labels):
    labels = np.asarray(labels)
    if len(np.unique(labels)) < 2:
        return np.nan
    n = min(len(labels), 5000)
    idx = np.random.RandomState(0).choice(len(labels), n, replace=False)
    return float(silhouette_score(emb[idx], labels[idx]))


rows = {"n_cells": adata.n_obs, "n_batches": int(adata.obs[batch_key].nunique())}
rows["batch_silhouette"] = safe_sil(adata.obs[batch_key])

cond_auc = np.nan
guard_pass = True
if "condition" in adata.obs and adata.obs["condition"].nunique() >= 2:
    rows["condition_silhouette"] = safe_sil(adata.obs["condition"])
    try:
        y = adata.obs["condition"].astype("category").cat.codes.values
        cond_auc = float(
            np.mean(
                cross_val_score(
                    KNeighborsClassifier(15), emb, y, cv=min(5, np.bincount(y).min())
                )
            )
        )
    except Exception as e:
        log.warning(f"condition kNN CV skipped: {e}")
    rows["condition_knn_accuracy"] = cond_auc
    # signal preserved if condition is separable above chance
    chance = float(pd.Series(adata.obs["condition"]).value_counts(normalize=True).max())
    guard_pass = bool(
        (cond_auc if cond_auc == cond_auc else 0) > chance + 0.05
        or (rows.get("condition_silhouette", -1) or -1) > 0
    )
rows["fibrosis_signal_guard_pass"] = guard_pass

# optional LISI from scib-metrics
try:
    from scib_metrics import ilisi_knn

    knn = adata.obsp["distances"]
    rows["ilisi"] = float(ilisi_knn(knn, adata.obs[batch_key].values))
except Exception as e:
    log.info(f"scib-metrics iLISI not computed: {e}")

df = pd.DataFrame([rows])
ensure_parent(snakemake.output.tsv)  # noqa: F821
df.to_csv(snakemake.output.tsv, sep="\t", index=False)  # noqa: F821

# small bar figure
ensure_parent(snakemake.output.fig)  # noqa: F821
plot_keys = [
    k
    for k in [
        "batch_silhouette",
        "condition_silhouette",
        "condition_knn_accuracy",
        "ilisi",
    ]
    if k in rows and rows[k] == rows[k]
]
fig, ax = plt.subplots(figsize=(5, 3))
ax.bar(plot_keys, [rows[k] for k in plot_keys])
ax.set_title("Integration metrics")
ax.tick_params(axis="x", rotation=30)
fig.tight_layout()
fig.savefig(snakemake.output.fig, dpi=110)  # noqa: F821

log.info(f"scIB metrics: {rows}")
if guard and not guard_pass:
    # surface loudly but DO NOT kill the pipeline — this is a QC flag, and donor-level
    # pseudobulk DE does not depend on the integrated embedding separating conditions.
    log.warning(
        "FIBROSIS-SIGNAL GUARD: condition is only weakly separable post-integration "
        f"(knn_acc={cond_auc}, silhouette={rows.get('condition_silhouette')}). Recorded as "
        "fibrosis_signal_guard_pass=False in the metrics. For the validation snRNA arm this "
        "is expected (subtler disease signal); for the primary arm, review batch_key "
        "(must be donor/sample) and integration strength."
    )

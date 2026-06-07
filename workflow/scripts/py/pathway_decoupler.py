"""Stage 06b — pathway & TF activity + per-cell signature scores (mechanism, Q7).

PROGENy pathway activity and CollecTRI TF activity via decoupler (fetched from OmniPath),
aggregated per compartment x condition; plus AUCell scores for the curated fibrosis gene
sets (local GMT, always offline-capable). Robust to OmniPath being unreachable.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import ensure_parent, get_logger, use_agg  # noqa: E402

import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402

use_agg()
import matplotlib.pyplot as plt  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
p = snakemake.params  # noqa: F821
adata = sc.read_h5ad(snakemake.input.h5ad)  # noqa: F821
adata.X = adata.layers.get("lognorm", adata.X)


def read_gmt_long(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                for g in parts[2:]:
                    rows.append((parts[0], g, 1.0))
    return pd.DataFrame(rows, columns=["source", "target", "weight"])


import decoupler as dc  # noqa: E402

group = [c for c in ["compartment", "condition"] if c in adata.obs]
records = []


def aggregate(estimate_df, level):
    """mean activity per compartment x condition for an obsm estimate DataFrame."""
    est = estimate_df.copy()
    est.index = adata.obs_names
    meta = adata.obs[group].copy()
    merged = meta.join(est)
    grouped = merged.groupby(group, observed=True).mean(numeric_only=True)
    long = grouped.reset_index().melt(
        id_vars=group, var_name="source", value_name="mean_activity"
    )
    long["level"] = level
    return long


# PROGENy pathway activity
try:
    prog = dc.get_progeny(organism="human", top=int(p.progeny_top))
    dc.run_mlm(
        mat=adata,
        net=prog,
        source="source",
        target="target",
        weight="weight",
        verbose=False,
        use_raw=False,
    )
    records.append(aggregate(adata.obsm["mlm_estimate"], "pathway_progeny"))
    log.info("PROGENy pathway activity computed")
except Exception as e:
    log.warning(f"PROGENy skipped (OmniPath unreachable?): {e}")

# CollecTRI TF activity
if bool(p.collectri):
    try:
        net = dc.get_collectri(organism="human", split_complexes=False)
        dc.run_ulm(
            mat=adata,
            net=net,
            source="source",
            target="target",
            weight="weight",
            verbose=False,
            use_raw=False,
        )
        records.append(aggregate(adata.obsm["ulm_estimate"], "tf_collectri"))
        log.info("CollecTRI TF activity computed")
    except Exception as e:
        log.warning(f"CollecTRI skipped: {e}")

# AUCell scores for curated fibrosis sets (offline)
geneset_long = read_gmt_long(os.path.join(p.genesets_dir, "fibrosis_core.gmt"))
try:
    dc.run_aucell(
        mat=adata,
        net=geneset_long,
        source="source",
        target="target",
        verbose=False,
        use_raw=False,
    )
    scores = aggregate(adata.obsm["aucell_estimate"], "geneset_aucell")
except Exception as e:
    log.warning(f"AUCell failed, using mean-expression score: {e}")
    rows = []
    for s in geneset_long["source"].unique():
        genes = [
            g
            for g in geneset_long.loc[geneset_long.source == s, "target"]
            if g in adata.var_names
        ]
        if genes:
            sc.tl.score_genes(adata, genes, score_name=f"gs_{s}")
            rows.append((s, f"gs_{s}"))
    if rows:
        est = adata.obs[[c for _, c in rows]].copy()
        est.columns = [s for s, _ in rows]
        scores = aggregate(est, "geneset_aucell")
    else:
        scores = pd.DataFrame(columns=group + ["source", "mean_activity", "level"])

# write activity (pathway + TF) and geneset score tables
ensure_parent(snakemake.output.activity)  # noqa: F821
activity = (
    pd.concat(records, ignore_index=True)
    if records
    else pd.DataFrame(columns=group + ["source", "mean_activity", "level"])
)
activity.round(4).to_csv(snakemake.output.activity, sep="\t", index=False)  # noqa: F821
scores.round(4).to_csv(snakemake.output.scores, sep="\t", index=False)  # noqa: F821

# heatmap of pathway activity by compartment (condition-averaged) if available
ensure_parent(snakemake.output.fig)  # noqa: F821
try:
    src = (
        activity[activity.level == "pathway_progeny"]
        if not activity.empty
        else activity
    )
    if src.empty:
        src = scores
    piv = src.pivot_table(
        index="compartment", columns="source", values="mean_activity", aggfunc="mean"
    )
    fig, ax = plt.subplots(
        figsize=(min(12, 1 + 0.5 * piv.shape[1]), 1 + 0.4 * piv.shape[0])
    )
    im = ax.imshow(piv.values, aspect="auto", cmap="RdBu_r")
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels(piv.columns, rotation=90, fontsize=6)
    ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels(piv.index, fontsize=8)
    fig.colorbar(im, ax=ax, shrink=0.6)
    ax.set_title("Pathway/geneset activity by compartment")
    fig.tight_layout()
    fig.savefig(snakemake.output.fig, dpi=110)  # noqa: F821
except Exception as e:
    log.warning(f"activity figure fallback: {e}")
    plt.figure()
    plt.text(0.5, 0.5, "no pathway activity", ha="center")
    plt.savefig(snakemake.output.fig, dpi=90)  # noqa: F821

log.info(f"pathway activity rows={len(activity)}, geneset rows={len(scores)}")

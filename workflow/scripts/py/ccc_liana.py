"""Stage 07a — consensus ligand-receptor inference with LIANA (Q7).

Runs LIANA's rank-aggregate consensus over compartments, SEPARATELY per condition so the
next step can compute condition-differential communication. Applies an expression-proportion
filter (a first overinterpretation guard). Degrades to a header-only table if LIANA fails.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import ensure_parent, get_logger, set_all_seeds  # noqa: E402

import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
p = snakemake.params  # noqa: F821
set_all_seeds(int(p.seed))
expr_prop = float(p.expr_prop)
min_cells = int(p.min_cells)
out = snakemake.output.liana  # noqa: F821
ensure_parent(out)

adata = sc.read_h5ad(snakemake.input.h5ad)  # noqa: F821
adata.X = adata.layers.get("lognorm", adata.X)
groupby = "compartment" if "compartment" in adata.obs else "cell_type"


def run_one(ad_sub, label):
    import liana as li

    # keep groups with enough cells
    vc = ad_sub.obs[groupby].value_counts()
    keep = vc[vc >= min_cells].index
    ad_sub = ad_sub[ad_sub.obs[groupby].isin(keep)].copy()
    if ad_sub.obs[groupby].nunique() < 2:
        return None
    li.mt.rank_aggregate(
        ad_sub, groupby=groupby, expr_prop=expr_prop, use_raw=False, verbose=False
    )
    res = ad_sub.uns["liana_res"].copy()
    res["condition"] = label
    return res


frames = []
try:
    if "condition" in adata.obs and adata.obs["condition"].nunique() >= 2:
        for cond in adata.obs["condition"].unique():
            r = run_one(adata[adata.obs["condition"] == cond].copy(), str(cond))
            if r is not None:
                frames.append(r)
    if not frames:  # pooled fallback
        r = run_one(adata, "all")
        if r is not None:
            frames.append(r)
except Exception as e:
    log.warning(f"LIANA failed: {e}")

cols = [
    "source",
    "target",
    "ligand_complex",
    "receptor_complex",
    "specificity_rank",
    "magnitude_rank",
    "condition",
]
if frames:
    res = pd.concat(frames, ignore_index=True)
    keep_cols = [c for c in cols if c in res.columns]
    res[keep_cols].to_csv(out, sep="\t", index=False)
    log.info(f"LIANA: {len(res)} interactions across {len(frames)} condition(s)")
else:
    pd.DataFrame(columns=cols).to_csv(out, sep="\t", index=False)
    log.warning("LIANA produced no results; wrote empty table")

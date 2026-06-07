"""Stage 05a — aggregate a compartment's single cells to DONOR-LEVEL pseudobulk (Q5).

The donor is the replication unit (not the cell), which is the key to avoiding
pseudoreplication-inflated false discoveries (Squair 2021). Donors with fewer than
`min_cells_per_donor` cells in this compartment are dropped. Emits a genes x donors
count matrix and a donor-level coldata table for the R DE step.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import ensure_parent, get_logger  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
comp = snakemake.params.compartment  # noqa: F821
min_cells = int(snakemake.params.min_cells_per_donor)  # noqa: F821
agg = snakemake.params.aggregate  # noqa: F821

adata = sc.read_h5ad(snakemake.input.h5ad)  # noqa: F821
sub = adata[adata.obs["compartment"].astype(str) == comp]
log.info(f"[{comp}] {sub.n_obs} cells across {sub.obs['donor_id'].nunique()} donors")

counts = sub.layers["counts"]
counts = counts.tocsr() if hasattr(counts, "tocsr") else counts
genes = list(sub.var_names)

cols, coldata_rows = {}, []
for donor, idx in sub.obs.groupby("donor_id").indices.items():
    if len(idx) < min_cells:
        continue
    block = counts[idx]
    vec = (
        np.asarray(block.sum(axis=0)).ravel() if hasattr(block, "sum") else block.sum(0)
    )
    if agg == "mean":
        vec = vec / len(idx)
    cols[donor] = np.rint(vec).astype(int)
    o = sub.obs.iloc[idx[0]]
    row = dict(
        sample=donor,
        donor_id=donor,
        condition=str(o.get("condition", "NA")),
        fibrosis_axis=int(o.get("fibrosis_axis", 0)),
        fibrosis_bin=("high" if int(o.get("fibrosis_axis", 0)) >= 2 else "low"),
        sort_gate=str(o.get("sort_gate", "NA")),
        n_cells=int(len(idx)),
    )
    if "sex" in sub.obs:
        row["sex"] = str(o.get("sex", "NA"))
    coldata_rows.append(row)

ensure_parent(snakemake.output.counts)  # noqa: F821
if cols:
    counts_df = pd.DataFrame(cols, index=genes)
    coldata = pd.DataFrame(coldata_rows).set_index("sample").loc[counts_df.columns]
else:
    counts_df = pd.DataFrame(index=genes)
    coldata = pd.DataFrame(
        columns=[
            "donor_id",
            "condition",
            "fibrosis_axis",
            "fibrosis_bin",
            "sort_gate",
            "n_cells",
        ]
    )
    log.warning(
        f"[{comp}] no donors passed min_cells={min_cells}; writing empty pseudobulk"
    )

counts_df.index.name = "gene"
counts_df.to_csv(snakemake.output.counts, sep="\t")  # noqa: F821
coldata.index.name = "sample"
coldata.to_csv(snakemake.output.coldata, sep="\t")  # noqa: F821
log.info(
    f"[{comp}] pseudobulk: {counts_df.shape[1]} donors x {counts_df.shape[0]} genes; "
    f"groups={coldata['condition'].value_counts().to_dict() if not coldata.empty else {}}"
)

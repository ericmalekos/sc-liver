"""Stage 01 — load one sample's 10x MEX triplet into an AnnData, attach harmonized
metadata to every cell, and stash raw counts in a dedicated layer.

Reads the MEX triplet directly (scipy/pandas) rather than via scanpy.read_10x_mtx,
which makes it robust to 10x v2 (genes.tsv, 2 cols) vs v3 (features.tsv, 3 cols) and
to scanpy's format auto-detection quirks. Gzipped or plain files both work.
"""

import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import get_logger, write_h5ad  # noqa: E402

import anndata as ad  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.io  # noqa: E402
import scipy.sparse as sp  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
data_path = str(snakemake.params.data_path)  # noqa: F821
meta = dict(snakemake.params.meta)  # noqa: F821
out = snakemake.output.h5ad  # noqa: F821


def find(*names):
    for n in names:
        p = os.path.join(data_path, n)
        if os.path.isfile(p):
            return p
    raise FileNotFoundError(f"none of {names} in {data_path}")


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


mtx_p = find("matrix.mtx.gz", "matrix.mtx")
feat_p = find("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")
bc_p = find("barcodes.tsv.gz", "barcodes.tsv")

log.info(f"Ingesting {meta.get('sample_id')} from {data_path}")
with gzip.open(mtx_p, "rb") if mtx_p.endswith(".gz") else open(mtx_p, "rb") as fh:
    M = scipy.io.mmread(fh).tocsr()  # features x cells (10x convention)
X = sp.csr_matrix(M.T)  # cells x features

feat = pd.read_csv(feat_p, sep="\t", header=None)
gene_ids = feat[0].astype(str).values
gene_sym = feat[1].astype(str).values if feat.shape[1] > 1 else gene_ids
barcodes = pd.read_csv(bc_p, sep="\t", header=None)[0].astype(str).values

var = pd.DataFrame({"gene_ids": gene_ids}, index=pd.Index(gene_sym, name=None))
adata = ad.AnnData(X=X, obs=pd.DataFrame(index=barcodes), var=var)
adata.var_names_make_unique()

for key in [
    "sample_id",
    "dataset",
    "condition",
    "fibrosis_stage_raw",
    "fibrosis_axis",
    "sort_gate",
    "modality",
    "donor_id",
]:
    if key in meta and meta[key] is not None:
        adata.obs[key] = int(meta[key]) if key == "fibrosis_axis" else str(meta[key])

adata.obs_names = [f"{meta['sample_id']}::{bc}" for bc in adata.obs_names]
adata.layers["counts"] = adata.X.copy()
adata.var["mt"] = (
    pd.Index(adata.var_names).str.upper().str.startswith(("MT-", "MT.", "MTRNR"))
)

log.info(
    f"Loaded {adata.n_obs} cells x {adata.n_vars} genes "
    f"(features file: {os.path.basename(feat_p)}, {feat.shape[1]} cols)"
)
write_h5ad(adata, out)

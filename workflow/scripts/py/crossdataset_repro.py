"""Stage 09 — cross-dataset reproducibility (Q1 validation arm).

Compares the primary (GSE136103, cirrhotic vs healthy, scRNA) and validation (GSE244832,
F2+ vs low, snRNA) donor-aware DE by DIRECTION/RANK rather than raw expression, which is
robust to the scRNA/snRNA modality gap. Emits a per-(gene, compartment) reproducibility
score consumed by the biomarker scoring.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import ensure_parent, get_logger, use_agg  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

use_agg()
import matplotlib.pyplot as plt  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
primary_ds = snakemake.params.primary  # noqa: F821
valid_ds = list(snakemake.params.validation)  # noqa: F821
log.info(f"crossdataset reproducibility: primary={primary_ds}, validation={valid_ds}")


def parse(files):
    d = {}
    for f in files:
        comp = os.path.basename(os.path.dirname(f))
        df = pd.read_csv(f, sep="\t")
        if not df.empty and "gene" in df:
            d[comp] = df
    return d


prim = parse(list(snakemake.input.primary))  # noqa: F821
vald = parse(list(snakemake.input.validation)) if "validation" in snakemake.input.keys() else {}  # noqa: F821

rows, summary = [], []
for comp, pdf in prim.items():
    pdf = pdf.dropna(subset=["gene"]).drop_duplicates("gene").set_index("gene")
    vdf = vald.get(comp)
    if vdf is not None and not vdf.empty:
        vdf = vdf.dropna(subset=["gene"]).drop_duplicates("gene").set_index("gene")
        shared = pdf.index.intersection(vdf.index)
        rho = np.nan
        if len(shared) >= 5:
            rho, _ = spearmanr(pdf.loc[shared, "stat"], vdf.loc[shared, "stat"], nan_policy="omit")
        summary.append({"compartment": comp, "n_shared": len(shared), "spearman_stat": rho})
    else:
        shared = pdf.index[:0]
        summary.append({"compartment": comp, "n_shared": 0, "spearman_stat": np.nan})

    for g in pdf.index:
        plfc = pdf.loc[g, "log2FoldChange"]
        ppadj = pdf.loc[g].get("padj", np.nan)
        in_val = vdf is not None and g in vdf.index
        vlfc = vdf.loc[g, "log2FoldChange"] if in_val else np.nan
        vpadj = vdf.loc[g].get("padj", np.nan) if in_val else np.nan
        concord = bool(in_val and np.sign(plfc) == np.sign(vlfc))
        if not in_val:
            score = 0.30
        elif concord and (vpadj is not None and vpadj == vpadj and vpadj < 0.25):
            score = 1.00
        elif concord:
            score = 0.60
        else:
            score = 0.00
        rows.append(dict(gene=g, compartment=comp, primary_lfc=plfc, primary_padj=ppadj,
                         valid_lfc=vlfc, valid_padj=vpadj, in_validation=in_val,
                         sign_concordant=concord, repro_score=score))

out_df = pd.DataFrame(rows)
ensure_parent(snakemake.output.repro)  # noqa: F821
out_df.to_csv(snakemake.output.repro, sep="\t", index=False)  # noqa: F821

# concordance scatter (never let plotting kill the rule)
ensure_parent(snakemake.output.fig)  # noqa: F821
try:
    fig, ax = plt.subplots(figsize=(5, 5))
    sub = out_df[out_df["in_validation"]]
    if not sub.empty:
        ax.scatter(sub["primary_lfc"], sub["valid_lfc"], s=8, alpha=0.5)
        ax.axhline(0, c="grey", lw=0.5)
        ax.axvline(0, c="grey", lw=0.5)
        ax.set_xlabel("primary log2FC")
        ax.set_ylabel("validation log2FC")
        ax.set_title("Cross-dataset DE concordance")
    else:
        ax.text(0.5, 0.5, "no validation overlap", ha="center")
    fig.tight_layout()
    fig.savefig(snakemake.output.fig, dpi=110)  # noqa: F821
except Exception as e:
    log.warning(f"concordance figure failed ({e}); writing placeholder")
    plt.figure()
    plt.text(0.5, 0.5, "figure unavailable", ha="center")
    plt.savefig(snakemake.output.fig, dpi=90)  # noqa: F821

log.info(f"repro: {len(out_df)} gene-compartment rows; "
         f"summary={pd.DataFrame(summary).to_dict('records')}")

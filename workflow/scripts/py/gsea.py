"""Stage 06a — GSEA on the donor-aware DE ranking (mechanism, Q7).

Ranks genes by the DESeq2 Wald stat and runs preranked GSEA against the curated fibrosis
gene sets (+ optional MSigDB Hallmark). Writes a per-compartment enrichment table.
Degrades gracefully when DE is empty or gene sets are tiny.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import ensure_parent, get_logger  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
seed = int(snakemake.params.seed)  # noqa: F821
genesets_dir = snakemake.params.genesets_dir  # noqa: F821
gsea_sets = list(snakemake.params.gsea_sets)  # noqa: F821
out = snakemake.output.gsea  # noqa: F821
ensure_parent(out)


def empty(note):
    pd.DataFrame(columns=["Term", "ES", "NES", "pval", "fdr", "lead_genes", "status"]).assign(
        status=note
    ).to_csv(out, sep="\t", index=False)
    log.info(f"empty GSEA: {note}")
    sys.exit(0)


def read_gmt(path):
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                sets[parts[0]] = parts[2:]
    return sets


de = pd.read_csv(snakemake.input.de, sep="\t")  # noqa: F821
if de.empty or "gene" not in de or de["gene"].isna().all():
    empty("no DE genes")

de = de.dropna(subset=["gene"])
rank_col = "stat" if de["stat"].notna().any() else "log2FoldChange"
rnk = (
    de[["gene", rank_col]].dropna().drop_duplicates("gene").set_index("gene")[rank_col].sort_values(ascending=False)
)
if rnk.shape[0] < 5:
    empty("too few ranked genes")

gene_sets = {}
if "fibrosis_core" in gsea_sets:
    gene_sets.update(read_gmt(os.path.join(genesets_dir, "fibrosis_core.gmt")))
if "hallmark" in gsea_sets:
    try:
        import gseapy as gp

        gene_sets.update(gp.get_library("MSigDB_Hallmark_2020"))
    except Exception as e:
        log.warning(f"Hallmark fetch failed (offline?): {e}")

try:
    import gseapy as gp

    pre = gp.prerank(
        rnk=rnk.reset_index(),
        gene_sets=gene_sets,
        min_size=2,
        max_size=1000,
        permutation_num=200,
        seed=seed,
        no_plot=True,
        threads=2,
        outdir=None,
    )
    res = pre.res2d.rename(
        columns={"NOM p-val": "pval", "FDR q-val": "fdr", "Lead_genes": "lead_genes"}
    )
    keep = [c for c in ["Term", "ES", "NES", "pval", "fdr", "lead_genes"] if c in res.columns]
    res = res[keep]
    res["status"] = "ok"
    res.to_csv(out, sep="\t", index=False)
    log.info(f"GSEA: {len(res)} sets tested ({rank_col} ranking, {rnk.shape[0]} genes)")
except Exception as e:
    log.warning(f"GSEA failed: {e}")
    empty(f"gsea_error:{e}")

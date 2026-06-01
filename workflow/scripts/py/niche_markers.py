"""Stage 05c — niche (subpopulation) biomarkers.

Compartment-level pseudobulk DE averages over the disease-emergent subpopulations where
fibrosis biology actually lives (scar-associated macrophages: SPP1/TREM2/GPNMB; scar-
associated endothelium: PLVAP/ACKR1; activated HSC / myofibroblasts), so genuine subset
markers come out non-significant at the compartment level. Here, within each compartment we:

  1. subcluster on the batch-corrected embedding (Leiden);
  2. flag subclusters ENRICHED for disease donors -- donor-aware, and valid despite FACS
     sorting because the comparison is WITHIN a compartment (one sort gate), not across
     whole-tissue proportions;
  3. take the up-regulated one-vs-rest markers of those disease-associated subclusters as
     niche biomarkers.

Emits a per-(compartment, subcluster) summary and a per-(gene, compartment) niche-marker
table (best disease-niche marker per gene) consumed by build_features for scoring.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import ensure_parent, get_logger, set_all_seeds  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
p = snakemake.params  # noqa: F821
set_all_seeds(int(p.seed))
contrast = dict(p.contrast)
grp, ref, tst = contrast["group"], contrast["ref"], contrast["test"]
res = float(p.resolution)
min_cells = int(p.min_cells_subcluster)
enrich_min = float(p.enrichment_min_log2fc)
marker_padj = float(p.marker_padj)
f2plus = int(p.f2plus_threshold)
compartments = list(p.compartments)
seed = int(p.seed)

adata = sc.read_h5ad(snakemake.input.h5ad)  # noqa: F821

# resolve the contrast group column (derive fibrosis_bin from the ordinal axis if needed)
if grp not in adata.obs.columns and grp == "fibrosis_bin" and "fibrosis_axis" in adata.obs:
    adata.obs["fibrosis_bin"] = np.where(adata.obs["fibrosis_axis"].astype(float) >= f2plus,
                                         "high", "low")

summary_rows, marker_rows = [], []
have_grp = grp in adata.obs.columns and adata.obs[grp].astype(str).nunique() >= 2
if not have_grp:
    log.warning(f"contrast group '{grp}' unusable; emitting empty niche tables")

adata.X = adata.layers.get("lognorm", adata.X)
emb = "X_pca_harmony" if "X_pca_harmony" in adata.obsm else ("X_emb" if "X_emb" in adata.obsm else None)

for comp in compartments:
    sub = adata[adata.obs["compartment"].astype(str) == comp].copy()
    if not have_grp or sub.n_obs < max(3 * min_cells, 60):
        continue
    # subcluster on the batch-corrected embedding (falls back to a fresh PCA)
    if emb is not None:
        sc.pp.neighbors(sub, use_rep=emb, n_neighbors=15, random_state=seed)
    else:
        sc.pp.pca(sub, n_comps=min(30, sub.n_obs - 1, sub.n_vars - 1), random_state=seed)
        sc.pp.neighbors(sub, random_state=seed)
    sc.tl.leiden(sub, resolution=res, key_added="subcluster", random_state=seed)
    sub.obs["subcluster"] = comp + ":" + sub.obs["subcluster"].astype(str)
    nsub = sub.obs["subcluster"].nunique()
    log.info(f"[{comp}] {sub.n_obs} cells -> {nsub} subclusters")

    # donor-aware disease enrichment: per-donor fraction of this compartment's cells in each
    # subcluster, averaged within each condition (robust to per-donor depth & to sorting).
    d = sub.obs[["donor_id", grp, "subcluster"]].astype(str)
    per_donor = (d.groupby(["donor_id", "subcluster"]).size()
                 / d.groupby("donor_id").size()).rename("frac").reset_index()
    donor_cond = d.drop_duplicates("donor_id").set_index("donor_id")[grp]
    per_donor["cond"] = per_donor["donor_id"].map(donor_cond)
    enrich = {}
    for sc_id, g in per_donor.groupby("subcluster"):
        f_tst = g.loc[g["cond"] == tst, "frac"].mean()
        f_ref = g.loc[g["cond"] == ref, "frac"].mean()
        f_tst = 0.0 if f_tst != f_tst else f_tst
        f_ref = 0.0 if f_ref != f_ref else f_ref
        log2fc = float(np.log2((f_tst + 1e-3) / (f_ref + 1e-3)))
        n_cells = int((sub.obs["subcluster"] == sc_id).sum())
        assoc = bool(log2fc >= enrich_min and n_cells >= min_cells)
        enrich[sc_id] = assoc
        summary_rows.append(dict(compartment=comp, subcluster=sc_id, n_cells=n_cells,
                                 frac_test=round(f_tst, 4), frac_ref=round(f_ref, 4),
                                 disease_log2fc=round(log2fc, 3), disease_associated=assoc))

    # one-vs-rest up-markers; keep only those marking a DISEASE-ASSOCIATED subcluster
    if nsub < 2 or not any(enrich.values()):
        continue
    sc.tl.rank_genes_groups(sub, "subcluster", method="wilcoxon", n_genes=300, use_raw=False)
    for sc_id, assoc in enrich.items():
        if not assoc:
            continue
        df = sc.get.rank_genes_groups_df(sub, group=sc_id)
        up = df[(df["logfoldchanges"] > 0) & (df["pvals_adj"] < marker_padj)]
        for _, r in up.iterrows():
            marker_rows.append(dict(gene=r["names"], compartment=comp, subcluster=sc_id,
                                    niche_lfc=float(r["logfoldchanges"]),
                                    niche_padj=float(r["pvals_adj"]),
                                    niche_score=float(r["scores"])))

summary = pd.DataFrame(summary_rows)
markers = pd.DataFrame(marker_rows)
ensure_parent(snakemake.output.summary)  # noqa: F821
summary.to_csv(snakemake.output.summary, sep="\t", index=False)

# collapse to the best disease-niche marker per (gene, compartment) for scoring
if not markers.empty:
    markers["abs_score"] = markers["niche_score"].abs()
    best = (markers.sort_values("abs_score", ascending=False)
            .drop_duplicates(["gene", "compartment"])
            .drop(columns="abs_score"))
    best["in_disease_niche"] = True
else:
    best = pd.DataFrame(columns=["gene", "compartment", "subcluster", "niche_lfc",
                                 "niche_padj", "niche_score", "in_disease_niche"])
ensure_parent(snakemake.output.markers)  # noqa: F821
best.to_csv(snakemake.output.markers, sep="\t", index=False)

n_assoc = int(summary["disease_associated"].sum()) if not summary.empty else 0
log.info(f"niche: {n_assoc} disease-associated subclusters; "
         f"{len(best)} (gene x compartment) disease-niche markers")

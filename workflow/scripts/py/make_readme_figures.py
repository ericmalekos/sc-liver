"""Render polished README figures from existing pipeline outputs into docs/figures/.

Standalone (no heavy recompute): re-plots the integration / annotation UMAPs, the scIB metrics,
a donor-aware DE volcano, and the ranked candidates with clean, non-overlapping legends suitable
for the repo README. Run after a pipeline run:

    python workflow/scripts/py/make_readme_figures.py
"""
import os

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

PRIMARY = "gse136103"
OUT = "docs/figures"
os.makedirs(OUT, exist_ok=True)

ann = ad.read_h5ad(f"results/04_annotate/{PRIMARY}/annotated.h5ad")
umap = ann.obsm["X_umap"]
TAB20 = list(plt.cm.tab20.colors)
TAB10 = list(plt.cm.tab10.colors)


def _scatter_cat(ax, values, title, palette, legend=True, legend_title=None, s=2):
    v = pd.Series(np.asarray(values).astype(str))
    cats = sorted(v.unique())
    for i, c in enumerate(cats):
        m = (v.values == c)
        ax.scatter(umap[m, 0], umap[m, 1], s=s, color=palette[i % len(palette)],
                   linewidths=0, rasterized=True)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("UMAP1", fontsize=8); ax.set_ylabel("UMAP2", fontsize=8)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    if legend:
        handles = [Line2D([0], [0], marker="o", linestyle="", markersize=5,
                          color=palette[i % len(palette)]) for i in range(len(cats))]
        ax.legend(handles, cats, title=legend_title, bbox_to_anchor=(1.01, 1),
                  loc="upper left", fontsize=7, title_fontsize=8, frameon=False,
                  handletextpad=0.3, borderaxespad=0)
    return cats


# --- Section 3: integration UMAP (sample = batches mix | condition = preserved) ---
fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
_scatter_cat(axes[0], ann.obs["sample_id"], "By sample — batches intermix", TAB20, legend=False)
_scatter_cat(axes[1], ann.obs["condition"], "By condition — preserved",
             ["#d73027", "#4575b4"], legend=True, legend_title="condition")  # cirrhotic, healthy
fig.subplots_adjust(wspace=0.25)
fig.savefig(f"{OUT}/03_integration_umap.png", dpi=140, bbox_inches="tight")
plt.close(fig)

# --- Section 3: scIB metrics (skinny bars) ---
m = pd.read_csv(f"results/03_integrate/{PRIMARY}/scib_metrics.tsv", sep="\t").iloc[0]
rows = [("batch_silhouette", "batch silhouette  (→0 mixed)", "#888888"),
        ("condition_silhouette", "condition silhouette", "#d73027"),
        ("condition_knn_accuracy", "condition kNN acc  (→1 kept)", "#1a9850")]
vals = [float(m[k]) for k, _, _ in rows]
fig, ax = plt.subplots(figsize=(5.2, 1.9))
ax.barh(range(len(vals)), vals, height=0.55, color=[c for _, _, c in rows])
ax.set_yticks(range(len(vals))); ax.set_yticklabels([l for _, l, _ in rows], fontsize=8)
ax.invert_yaxis(); ax.axvline(0, color="k", lw=0.6)
ax.set_xlim(min(-0.05, min(vals) - 0.1), max(vals) + 0.18)
for i, v in enumerate(vals):
    ax.text(v + (0.015 if v >= 0 else -0.015), i, f"{v:.2f}", va="center",
            ha="left" if v >= 0 else "right", fontsize=7.5)
ax.set_title(f"Integration metrics — {PRIMARY}", fontsize=10)
ax.tick_params(labelsize=7.5)
for sp in ["top", "right"]:
    ax.spines[sp].set_visible(False)
fig.savefig(f"{OUT}/03_integration_scib.png", dpi=140, bbox_inches="tight")
plt.close(fig)

# --- Section 4: annotation UMAP (compartments, legend outside) ---
fig, ax = plt.subplots(figsize=(7.4, 5))
_scatter_cat(ax, ann.obs["compartment"], "Cell compartments — primary cohort",
             TAB10, legend=True, legend_title="compartment")
fig.savefig(f"{OUT}/04_annotation_umap.png", dpi=140, bbox_inches="tight")
plt.close(fig)

# --- Section 5: donor-aware DE volcano (stellate/myofibroblast) ---
de = pd.read_csv(f"results/05_de/{PRIMARY}/stellate_myofibroblast/de_results.tsv", sep="\t")
de = de.dropna(subset=["padj", "log2FoldChange"]).copy()
de["nlp"] = -np.log10(de["padj"].clip(lower=1e-300))
sig = (de["padj"] < 0.05) & (de["log2FoldChange"].abs() >= 0.58)
fig, ax = plt.subplots(figsize=(6, 5))
ax.scatter(de.loc[~sig, "log2FoldChange"], de.loc[~sig, "nlp"], s=4, color="#d9d9d9",
           linewidths=0, rasterized=True)
up = sig & (de["log2FoldChange"] > 0); dn = sig & (de["log2FoldChange"] < 0)
ax.scatter(de.loc[up, "log2FoldChange"], de.loc[up, "nlp"], s=7, color="#d73027",
           linewidths=0, label=f"up ({int(up.sum())})", rasterized=True)
ax.scatter(de.loc[dn, "log2FoldChange"], de.loc[dn, "nlp"], s=7, color="#4575b4",
           linewidths=0, label=f"down ({int(dn.sum())})", rasterized=True)
for _, r in de[up].nlargest(8, "nlp").iterrows():
    ax.annotate(r["gene"], (r["log2FoldChange"], r["nlp"]), fontsize=7, ha="left", va="bottom")
ax.axvline(0.58, ls="--", lw=0.5, color="grey"); ax.axvline(-0.58, ls="--", lw=0.5, color="grey")
ax.axhline(-np.log10(0.05), ls="--", lw=0.5, color="grey")
ax.set_xlabel("log2 fold-change (cirrhotic vs healthy)", fontsize=9)
ax.set_ylabel("-log10 adjusted p", fontsize=9)
ax.set_title("Donor-aware pseudobulk DE — stellate / myofibroblast", fontsize=10)
ax.tick_params(labelsize=8); ax.legend(fontsize=8, frameon=False, loc="upper left")
for sp in ["top", "right"]:
    ax.spines[sp].set_visible(False)
fig.savefig(f"{OUT}/05_de_volcano.png", dpi=140, bbox_inches="tight")
plt.close(fig)

# --- Section 6: PROGENy pathway-activity shift in fibrosis (via decoupler) ---
act = pd.read_csv(f"results/06_pathway/{PRIMARY}/activity.tsv", sep="\t")
prog = act[act["level"] == "pathway_progeny"]
comps6 = sorted(prog["compartment"].unique())
paths = sorted(prog["source"].unique())
delta = np.full((len(comps6), len(paths)), np.nan)
for i, c in enumerate(comps6):
    for j, pth in enumerate(paths):
        sub = prog[(prog["compartment"] == c) & (prog["source"] == pth)]
        ci = sub.loc[sub["condition"] == "cirrhotic", "mean_activity"]
        he = sub.loc[sub["condition"] == "healthy", "mean_activity"]
        if len(ci) and len(he):
            delta[i, j] = ci.mean() - he.mean()
fig, ax = plt.subplots(figsize=(8.2, 3.0))
vmax = np.nanmax(np.abs(delta)) or 1.0
im = ax.imshow(delta, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks(range(len(paths))); ax.set_xticklabels(paths, rotation=45, ha="right", fontsize=8)
ax.set_yticks(range(len(comps6))); ax.set_yticklabels(comps6, fontsize=8)
cb = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
cb.set_label("activity shift\n(cirrhotic - healthy)", fontsize=7.5); cb.ax.tick_params(labelsize=7)
ax.set_title("PROGENy pathway-activity shift in fibrosis (decoupler), per compartment", fontsize=10)
fig.savefig(f"{OUT}/06_pathway_activity.png", dpi=140, bbox_inches="tight")
plt.close(fig)

# --- Section 7: cross-dataset concordance (per compartment, labeled contrasts) ---
from scipy.stats import spearmanr  # noqa: E402

rp = pd.read_csv("results/09_crossdataset/repro_scores.tsv", sep="\t")
rp = rp[rp["in_validation"].astype(str).str.lower() == "true"].dropna(
    subset=["primary_lfc", "valid_lfc"])
comps7 = sorted(rp["compartment"].unique())
pal7 = {c: TAB10[i % len(TAB10)] for i, c in enumerate(comps7)}
fig, ax = plt.subplots(figsize=(6.4, 5.8))
ax.axline((0, 0), (1, 1), color="#bbbbbb", lw=0.9, ls="--", zorder=0)  # y = x (perfect agreement)
ax.axhline(0, color="#e8e8e8", lw=0.6); ax.axvline(0, color="#e8e8e8", lw=0.6)
handles7 = []
for c in comps7:
    d = rp[rp["compartment"] == c]
    ax.scatter(d["primary_lfc"], d["valid_lfc"], s=6, color=pal7[c], alpha=0.5,
               linewidths=0, rasterized=True)
    rho = spearmanr(d["primary_lfc"], d["valid_lfc"], nan_policy="omit").correlation \
        if len(d) >= 5 else float("nan")
    handles7.append(Line2D([0], [0], marker="o", linestyle="", markersize=6, color=pal7[c],
                           label=f"{c}  (rho={rho:.2f}, n={len(d)})"))
ax.set_xlim(-5, 5); ax.set_ylim(-5, 5)
ax.set_xlabel("primary log2FC   (GSE136103: cirrhotic vs healthy, scRNA)", fontsize=8.5)
ax.set_ylabel("validation log2FC   (GSE244832: MASH F2+ vs low, snRNA)", fontsize=8.5)
ax.set_title("Cross-dataset DE concordance, matched per compartment", fontsize=10)
ax.legend(handles=handles7, title="compartment (Spearman rho)", fontsize=7,
          title_fontsize=7.5, frameon=False, loc="upper left")
ax.tick_params(labelsize=8)
for sp in ["top", "right"]:
    ax.spines[sp].set_visible(False)
fig.savefig(f"{OUT}/07_crossdataset_concordance.png", dpi=140, bbox_inches="tight")
plt.close(fig)

# --- Section 8: ranked candidates with compartment color legend ---
cs = pd.read_csv("results/08_score/candidate_scores.tsv", sep="\t")
sel = cs[cs["selected"].astype(str).str.lower() == "true"].sort_values("composite")
comps = sorted(sel["compartment"].unique())
cmap = {c: TAB10[i % len(TAB10)] for i, c in enumerate(comps)}
fig, ax = plt.subplots(figsize=(8, 6.4))
ax.barh(range(len(sel)), sel["composite"], color=[cmap[c] for c in sel["compartment"]], height=0.72)
ax.set_yticks(range(len(sel))); ax.set_yticklabels(sel["gene"], fontsize=8)
ax.set_xlabel("composite score", fontsize=9)
ax.set_title("Top prioritized biomarker / target candidates", fontsize=11)
ax.tick_params(labelsize=8)
handles = [Line2D([0], [0], marker="s", linestyle="", markersize=8, color=cmap[c]) for c in comps]
ax.legend(handles, comps, title="compartment", fontsize=8, title_fontsize=8,
          loc="lower right", frameon=True)
for sp in ["top", "right"]:
    ax.spines[sp].set_visible(False)
fig.savefig(f"{OUT}/08_top_candidates.png", dpi=140, bbox_inches="tight")
plt.close(fig)

print("README figures written to", OUT)

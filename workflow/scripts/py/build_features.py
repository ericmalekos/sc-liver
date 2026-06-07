"""Stage 08a — assemble the (gene x compartment) feature matrix for scoring (Q6).

Features: donor-aware DE effect/significance, cell-type specificity (tau), cross-dataset
reproducibility, druggability (cached Open Targets/DGIdb snapshot), biomarker accessibility
(secretome/surfaceome), and fibrosis-geneset membership. The ML/SHAP signal is merged later.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import ensure_parent, get_logger  # noqa: E402

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
p = snakemake.params  # noqa: F821
TOP_PER_COMP = 1000  # cap candidates/compartment by |stat| to bound feature-matrix size


def tau_specificity(adata):
    """tau index per gene across compartment mean expression (1 = compartment-specific)."""
    adata.X = adata.layers.get("lognorm", adata.X)
    comps = adata.obs["compartment"].astype(str)
    means = {}
    for c in comps.unique():
        m = np.asarray(adata[comps == c].X.mean(axis=0)).ravel()
        means[c] = m
    M = pd.DataFrame(means, index=adata.var_names)  # genes x compartments
    xmax = M.max(axis=1).replace(0, np.nan)
    tau = ((1 - M.div(xmax, axis=0)).sum(axis=1)) / (M.shape[1] - 1)
    home = M.idxmax(axis=1)
    ratio = M.div(xmax, axis=0)  # genes x compartments; 1.0 at the home compartment
    return tau.fillna(0.0), home, ratio


# --- DE per compartment ---
de_rows = []
for f in list(snakemake.input.de):  # noqa: F821
    comp = os.path.basename(os.path.dirname(f))
    d = pd.read_csv(f, sep="\t")
    if d.empty or "gene" not in d:
        continue
    d = d.dropna(subset=["gene"]).copy()
    d["compartment"] = comp
    d["abs_stat"] = d["stat"].abs()
    d = d.sort_values("abs_stat", ascending=False).head(TOP_PER_COMP)
    de_rows.append(d)
de = (
    pd.concat(de_rows, ignore_index=True)
    if de_rows
    else pd.DataFrame(
        columns=["gene", "compartment", "log2FoldChange", "stat", "padj", "baseMean"]
    )
)
de["direction"] = np.where(de.get("log2FoldChange", 0) >= 0, "up", "down")

# --- specificity ---
adata = ad.read_h5ad(snakemake.input.annotated)  # noqa: F821
tau, home, ratio = tau_specificity(adata)
de["tau"] = de["gene"].map(tau).fillna(0.0)
de["home_compartment"] = de["gene"].map(home)
de["is_home"] = de["home_compartment"] == de["compartment"]


# per-(gene, compartment) relative expression (mean here / gene's peak compartment mean):
# 1.0 at the home compartment, lower elsewhere. Distinguishes a genuinely multi-compartment
# gene (MMP2 in stellate ~0.6) from ambient leakage into a compartment (SPP1 in lymphoid ~0.15).
def _ratio(g, c):
    try:
        v = ratio.at[g, c]
        return float(v) if v == v else 0.0
    except (KeyError, TypeError, ValueError):
        return 0.0


de["specificity_ratio"] = [_ratio(g, c) for g, c in zip(de["gene"], de["compartment"])]

# --- reproducibility ---
if "repro" in snakemake.input.keys():  # noqa: F821
    repro = pd.read_csv(snakemake.input.repro, sep="\t")  # noqa: F821
    de = de.merge(
        repro[
            ["gene", "compartment", "repro_score", "sign_concordant", "in_validation"]
        ],
        on=["gene", "compartment"],
        how="left",
    )
de["repro_score"] = de.get("repro_score", pd.Series(index=de.index)).fillna(0.3)

# --- niche (subpopulation) markers: disease-associated-subcluster signal per gene ---
if "niche" in snakemake.input.keys():  # noqa: F821
    nm = pd.read_csv(snakemake.input.niche, sep="\t")  # noqa: F821
    if not nm.empty:
        de = de.merge(
            nm[["gene", "compartment", "niche_lfc", "niche_padj", "in_disease_niche"]],
            on=["gene", "compartment"],
            how="left",
        )
de["in_disease_niche"] = (
    de.get("in_disease_niche", pd.Series(index=de.index)).fillna(False).astype(bool)
)
de["niche_lfc"] = de.get("niche_lfc", pd.Series(index=de.index)).fillna(0.0)
de["niche_padj"] = de.get("niche_padj", pd.Series(index=de.index)).fillna(1.0)

# --- druggability (cached snapshot; fallback default low) ---
drug = {}
cache = p.druggability_cache
if os.path.exists(cache):
    dd = pd.read_csv(cache, sep="\t")
    drug = dict(zip(dd["gene"], dd["druggability"]))
de["druggability"] = de["gene"].map(drug).fillna(0.10)


# --- accessibility (secretome > surfaceome > intracellular) ---
def load_list(path):
    return set(pd.read_csv(path)["gene"]) if path and os.path.exists(path) else set()


secret = load_list(p.accessibility.get("secretome"))
surf = load_list(p.accessibility.get("surfaceome"))


def access(g):
    if g in secret:
        return ("secreted", 1.0)
    if g in surf:
        return ("surface", 0.8)
    return ("intracellular", 0.3)


acc = de["gene"].map(access)
de["accessibility_class"] = [a[0] for a in acc]
de["accessibility_score"] = [a[1] for a in acc]

# --- fibrosis geneset membership (metadata) ---
fib_genes = set()
gmt = os.path.join("config/genesets/fibrosis_core.gmt")
if os.path.exists(gmt):
    for line in open(gmt):
        parts = line.rstrip("\n").split("\t")
        fib_genes |= set(parts[2:])
de["fibrosis_geneset_member"] = de["gene"].isin(fib_genes).astype(int)

ensure_parent(snakemake.output.features)  # noqa: F821
keep = [
    "gene",
    "compartment",
    "log2FoldChange",
    "stat",
    "padj",
    "baseMean",
    "direction",
    "tau",
    "home_compartment",
    "is_home",
    "specificity_ratio",
    "repro_score",
    "druggability",
    "accessibility_class",
    "accessibility_score",
    "fibrosis_geneset_member",
    "niche_lfc",
    "niche_padj",
    "in_disease_niche",
]
de[[c for c in keep if c in de.columns]].to_csv(
    snakemake.output.features, sep="\t", index=False
)  # noqa: F821
log.info(
    f"feature matrix: {len(de)} (gene x compartment) rows across "
    f"{de['compartment'].nunique()} compartments"
)

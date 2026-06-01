"""Stage 08c — composite biomarker prioritization score + ranked candidate list (Q6).

Each component is min-max normalized WITHIN compartment, then combined with the configured
weights. The headline list takes the top-N per required compartment (so each disease
compartment is represented) and fills to top-N overall. Output includes every component and a
one-line rationale per candidate for transparency / explainability.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import ensure_parent, get_logger, use_agg  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

use_agg()
import matplotlib.pyplot as plt  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
p = snakemake.params  # noqa: F821
W = dict(p.weights)
required = list(p.required_compartments)
padj_thr = float(p.padj_threshold)
lfc_thr = float(p.lfc_threshold)
require_sig = bool(p.require_de_significance)
require_home = bool(p.require_home_compartment)
spec_floor = float(p.specificity_ratio_floor)

feat = pd.read_csv(snakemake.input.features, sep="\t")  # noqa: F821
shap = pd.read_csv(snakemake.input.shap, sep="\t")  # noqa: F821
if not shap.empty:
    feat = feat.merge(shap[["gene", "compartment", "shap_importance"]],
                      on=["gene", "compartment"], how="left")
feat["shap_importance"] = feat.get("shap_importance", pd.Series(index=feat.index)).fillna(0.0)

# raw component signals
feat["padj"] = feat["padj"].fillna(1.0)
feat["de_signal"] = feat["log2FoldChange"].abs() * -np.log10(feat["padj"].clip(lower=1e-300))
comp_raw = {
    "de": "de_signal",
    "specificity": "tau",
    "reproducibility": "repro_score",
    "druggability": "druggability",
    "accessibility": "accessibility_score",
    "ml_shap": "shap_importance",
}


def minmax_within(df, col):
    def mm(x):
        rng = x.max() - x.min()
        return (x - x.min()) / rng if rng > 0 else pd.Series(0.0, index=x.index)
    return df.groupby("compartment")[col].transform(mm)


for k, raw in comp_raw.items():
    feat[f"S_{k}"] = minmax_within(feat, raw)

feat["composite"] = sum(W[k] * feat[f"S_{k}"] for k in W)

# rationale string
def rationale(r):
    bits = [f"{r['direction']} in fibrosis (LFC={r['log2FoldChange']:.2f}, padj={r['padj']:.2g})"]
    if r["tau"] >= 0.7:
        bits.append(f"{r['home_compartment']}-specific (tau={r['tau']:.2f})")
    if r.get("repro_score", 0) >= 0.6:
        bits.append("reproduced in validation")
    if r["druggability"] >= 0.6:
        bits.append("druggable")
    bits.append(r["accessibility_class"])
    if r.get("shap_importance", 0) >= 0.5:
        bits.append("high ML importance")
    return "; ".join(bits)


feat["rationale"] = feat.apply(rationale, axis=1)

# candidates = up-regulated, DE-tested
cand = feat[(feat["direction"] == "up") & feat["padj"].notna()].copy()

# DE significance gate: only genes that genuinely change in disease are ELIGIBLE for the
# headline list. Specificity/druggability/accessibility must not float a non-DE lineage
# marker (e.g. padj~1) into the ranked candidates. Every up-regulated gene is still scored
# and written out (flagged by `de_significant`) for transparency.
cand["de_significant"] = (cand["padj"] < padj_thr) & (cand["log2FoldChange"].abs() >= lfc_thr)
pool = cand[cand["de_significant"]] if require_sig else cand
# specificity gate: a gene may be ranked in a compartment only where it is genuinely
# expressed -- relative expression (specificity_ratio) >= floor. The home compartment is
# always 1.0; a multi-compartment fibrosis gene (MMP2 in stellate ~0.6) clears the floor,
# while ambient leakage (a macrophage gene at ~0.15 in lymphoid) does not.
if require_home and "specificity_ratio" in pool.columns:
    pool = pool[pool["specificity_ratio"].astype(float) >= spec_floor]
elif require_home and "is_home" in pool.columns:  # fallback if ratio column absent
    pool = pool[pool["is_home"].astype(str).str.lower().isin(["true", "1", "1.0"])]
pool = pool.sort_values(["composite"] + [f"S_{t}" for t in p.tie_breaker], ascending=False)

# selection: top-N per required compartment, then fill to top-N overall (from eligible pool)
selected = []
for comp in required:
    sub = pool[pool["compartment"] == comp].head(int(p.top_n_per_compartment))
    selected.append(sub)
sel = pd.concat(selected) if selected else pool.head(0)
remaining = pool[~pool.index.isin(sel.index)]
fill = remaining.head(max(0, int(p.top_n_overall) - len(sel)))
sel = pd.concat([sel, fill]).sort_values("composite", ascending=False)
sel = sel.head(int(p.top_n_overall))
cand["selected"] = cand.index.isin(sel.index)
cand = cand.sort_values("composite", ascending=False).reset_index(drop=True)
cand.loc[cand["selected"], "rank"] = range(1, int(cand["selected"].sum()) + 1)

out_cols = ["rank", "gene", "compartment", "composite", "direction", "log2FoldChange", "padj",
            "S_de", "S_specificity", "S_reproducibility", "S_druggability", "S_accessibility",
            "S_ml_shap", "tau", "home_compartment", "is_home", "specificity_ratio",
            "druggability", "accessibility_class", "repro_score", "de_significant",
            "selected", "rationale"]
ensure_parent(snakemake.output.scores)  # noqa: F821
cand[[c for c in out_cols if c in cand.columns]].to_csv(snakemake.output.scores, sep="\t", index=False)  # noqa: F821

# figure: top selected candidates
ensure_parent(snakemake.output.fig)  # noqa: F821
top = cand[cand["selected"]].sort_values("composite")
fig, ax = plt.subplots(figsize=(7, max(3, 0.32 * len(top))))
if not top.empty:
    comps = sorted(top["compartment"].unique())
    cmap = {c: plt.cm.tab10(i) for i, c in enumerate(comps)}
    ax.barh(top["gene"] + " (" + top["compartment"].str[:4] + ")", top["composite"],
            color=[cmap[c] for c in top["compartment"]])
    ax.set_xlabel("composite score")
    ax.set_title("Top prioritized biomarker / target candidates")
else:
    ax.text(0.5, 0.5, "no candidates", ha="center")
fig.tight_layout(); fig.savefig(snakemake.output.fig, dpi=120)  # noqa: F821

n_sel = int(cand["selected"].sum())
log.info(f"scored {len(cand)} up-regulated candidates; selected top {n_sel} "
         f"(per-compartment {dict(top['compartment'].value_counts())})")

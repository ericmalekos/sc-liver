"""Stage 08b — explainable ML signal for prioritization (Q6).

Per required compartment, trains a fibrotic-vs-healthy classifier on DONOR-level pseudobulk
(no cell leakage) with DONOR-GROUPED cross-validation, then extracts SHAP feature importances
to rank genes by their contribution. Robust to tiny/degenerate compartments (importance 0).
"""

import json
import os
import pickle
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import ensure_parent, get_logger, set_all_seeds  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.model_selection import GroupKFold, cross_val_score  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
p = snakemake.params  # noqa: F821
set_all_seeds(int(p.seed))
contrast = dict(p.contrast)
group_col, ref, tst = contrast["group"], contrast["ref"], contrast["test"]

features = pd.read_csv(snakemake.input.features, sep="\t")  # noqa: F821
counts_files = list(snakemake.input.counts)  # noqa: F821
coldata_files = list(snakemake.input.coldata)  # noqa: F821
by_comp = {os.path.basename(os.path.dirname(f)): f for f in counts_files}
col_by_comp = {os.path.basename(os.path.dirname(f)): f for f in coldata_files}


def make_model():
    if p.model == "random_forest":
        from sklearn.ensemble import RandomForestClassifier

        return RandomForestClassifier(n_estimators=300, random_state=int(p.seed))
    from xgboost import XGBClassifier

    return XGBClassifier(
        n_estimators=200,
        max_depth=3,
        learning_rate=0.1,
        subsample=0.8,
        eval_metric="logloss",
        random_state=int(p.seed),
    )


rows, metrics, models = [], {}, {}
for comp, cfile in by_comp.items():
    counts = pd.read_csv(cfile, sep="\t", index_col=0)
    coldata = pd.read_csv(col_by_comp[comp], sep="\t", index_col=0)
    if counts.shape[1] < 4 or coldata.empty or group_col not in coldata:
        metrics[comp] = {"status": "skipped_small", "n_donors": int(counts.shape[1])}
        continue
    y = coldata[group_col].map({ref: 0, tst: 1})
    mask = y.notna()
    counts, y, coldata = counts.loc[:, mask.values], y[mask], coldata[mask.values]
    if y.nunique() < 2 or y.value_counts().min() < 2:
        metrics[comp] = {
            "status": "skipped_one_class",
            "n_donors": int(counts.shape[1]),
        }
        continue

    # candidate genes for this compartment (fallback to all)
    cand = features.loc[features["compartment"] == comp, "gene"].tolist()
    genes = [g for g in cand if g in counts.index] or list(counts.index)
    cpm = counts.loc[genes].T  # donors x genes
    cpm = np.log1p(cpm.div(cpm.sum(axis=1).replace(0, 1), axis=0) * 1e4)
    X, yv = cpm.values, y.values
    groups = coldata["donor_id"].values if "donor_id" in coldata else np.arange(len(yv))

    n_splits = int(
        min(int(p.cv_folds), len(np.unique(groups)), pd.Series(yv).value_counts().min())
    )
    cv_score = np.nan
    try:
        if n_splits >= 2:
            cv = GroupKFold(n_splits=n_splits)
            cv_score = float(
                np.mean(cross_val_score(make_model(), X, yv, groups=groups, cv=cv))
            )
    except Exception as e:
        log.warning(f"[{comp}] CV failed: {e}")

    model = make_model()
    model.fit(X, yv)
    models[comp] = model
    try:
        import shap

        expl = shap.TreeExplainer(model)
        sv = expl.shap_values(X)
        sv = sv[1] if isinstance(sv, list) else sv
        imp = np.abs(sv).mean(axis=0)
    except Exception as e:
        log.warning(f"[{comp}] SHAP fallback to feature_importances_: {e}")
        imp = getattr(model, "feature_importances_", np.zeros(len(genes)))
    imp = np.asarray(imp, dtype=float).ravel()
    rng = imp.max() - imp.min()
    norm = (imp - imp.min()) / rng if rng > 0 else np.zeros_like(imp)
    for g, s in zip(genes, norm):
        rows.append(
            {
                "gene": g,
                "compartment": comp,
                "shap_importance": float(s),
                "cv_score": cv_score,
            }
        )
    metrics[comp] = {
        "status": "ok",
        "n_donors": int(len(yv)),
        "n_features": len(genes),
        "cv_accuracy": cv_score,
    }
    log.info(f"[{comp}] trained on {len(yv)} donors, {len(genes)} genes, cv={cv_score}")

shap_df = pd.DataFrame(
    rows, columns=["gene", "compartment", "shap_importance", "cv_score"]
)
ensure_parent(snakemake.output.shap)  # noqa: F821
shap_df.to_csv(snakemake.output.shap, sep="\t", index=False)  # noqa: F821
with open(snakemake.output.model, "wb") as fh:  # noqa: F821
    pickle.dump(models, fh)
with open(snakemake.output.metrics, "w") as fh:  # noqa: F821
    json.dump(metrics, fh, indent=2)
log.info(f"ML done: {len(shap_df)} gene importances across {len(models)} compartments")

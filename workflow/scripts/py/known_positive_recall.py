"""Stage 08e — known-positive recall benchmark (Q6 sanity check).

Asks the blunt question a skeptic asks first: "does this pipeline recover the biomarkers we
already KNOW are real?" For each curated literature positive (resources/known_positive_markers.tsv)
it classifies the outcome AND the reason, so a miss is actionable rather than mysterious:

  selected                  -> in the final ranked candidate list (full success)
  scored_not_selected       -> scored as a candidate but below the top-N cutoff
  DE_significant_off_list   -> significant in primary DE but absent from the candidate list
  DE_significant_wrong_dir  -> significant in primary DE but opposite the expected direction
  DE_tested_nonsignificant  -> tested in primary DE, padj >= threshold (underpowered/diluted)
  DE_padj_NA                -> tested but padj=NA (DESeq2 independent filtering / low counts)
  detected_filtered_pre_DE  -> in the pseudobulk matrix but dropped by the expression filter
  absent_primary            -> not in the primary pseudobulk at all (detection floor; e.g. SMOC2)

Each row is independently cross-referenced against the validation dataset(s), so a known
positive that the primary misses but the validation confirms (the SMOC2 pattern) is flagged
explicitly as `rescued_by_validation` — the clearest signal of a primary-anchoring blind spot.
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
P = snakemake.params  # noqa: F821
padj_thr = float(P.padj_threshold)


def _read(f):
    try:
        df = pd.read_csv(f, sep="\t")
    except Exception:
        return pd.DataFrame()
    return df


def _comp_of(path):
    return os.path.basename(os.path.dirname(path))


def _by_comp(files):
    """{compartment: dataframe} keyed on the path's parent dir, skipping empties."""
    out = {}
    for f in files:
        df = _read(f)
        if not df.empty and "gene" in df.columns:
            out[_comp_of(f)] = df.set_index("gene")
    return out


# ---- load the pipeline's own outputs -----------------------------------------
panel = pd.read_csv(P.panel, sep="\t", comment="#")
cand = _read(snakemake.input.candidates)  # noqa: F821
if not cand.empty:
    cand["selected"] = cand["selected"].astype(str).str.lower().eq("true")
primary_de = _by_comp(list(snakemake.input.primary_de))  # noqa: F821
primary_counts = _by_comp(list(snakemake.input.primary_counts))  # noqa: F821
valid_de = _by_comp(list(snakemake.input.valid_de)) if "valid_de" in snakemake.input.keys() else {}  # noqa: F821


def _padj(row):
    v = row.get("padj", np.nan)
    return v if (v == v) else np.nan  # NaN-safe


def classify(gene, comp, direction):
    """Return (status, detail dict) for one known positive in its expected compartment."""
    want_up = str(direction).lower() == "up"
    d = dict(rank="", composite=np.nan, primary_lfc=np.nan, primary_padj=np.nan)

    # 1) candidate list (the headline output)
    if not cand.empty:
        hit = cand[(cand["gene"] == gene) & (cand["compartment"] == comp)]
        if not hit.empty:
            r = hit.iloc[0]
            d.update(composite=r.get("composite", np.nan),
                     primary_lfc=r.get("log2FoldChange", np.nan),
                     primary_padj=r.get("padj", np.nan))
            if bool(r["selected"]):
                d["rank"] = r.get("rank", "")
                return "selected", d
            return "scored_not_selected", d

    # 2) primary DE table for the compartment
    de = primary_de.get(comp)
    if de is not None and gene in de.index:
        row = de.loc[gene]
        lfc = row.get("log2FoldChange", np.nan)
        padj = _padj(row)
        d.update(primary_lfc=lfc, primary_padj=padj)
        if padj != padj:
            return "DE_padj_NA", d
        if padj < padj_thr:
            if (lfc > 0) == want_up:
                return "DE_significant_off_list", d
            return "DE_significant_wrong_dir", d
        return "DE_tested_nonsignificant", d

    # 3) was it even in the pseudobulk matrix? (filtered vs detection floor)
    counts = primary_counts.get(comp)
    if counts is not None and gene in counts.index:
        return "detected_filtered_pre_DE", d
    return "absent_primary", d


def validation_status(gene, comp, direction):
    """Best validation outcome for this gene across all validation datasets."""
    want_up = str(direction).lower() == "up"
    de = valid_de.get(comp)
    if de is None or gene not in de.index:
        return "val_absent", np.nan, np.nan
    row = de.loc[gene]
    lfc = row.get("log2FoldChange", np.nan)
    padj = _padj(row)
    if padj == padj and padj < padj_thr:
        return ("validated" if (lfc > 0) == want_up else "validated_wrong_dir"), lfc, padj
    return "val_tested_nonsig", lfc, padj


MISS_IN_PRIMARY = {"DE_tested_nonsignificant", "DE_padj_NA",
                   "detected_filtered_pre_DE", "absent_primary"}

rows = []
for _, m in panel.iterrows():
    g, comp, direction = m["gene"], m["compartment"], m["direction"]
    status, d = classify(g, comp, direction)
    vstat, vlfc, vpadj = validation_status(g, comp, direction)
    rescued = status in MISS_IN_PRIMARY and vstat == "validated"
    rows.append(dict(
        gene=g, compartment=comp, expected_direction=direction, source=m.get("source", ""),
        found_status=status, rank=d["rank"], composite=d["composite"],
        primary_lfc=d["primary_lfc"], primary_padj=d["primary_padj"],
        validation_status=vstat, valid_lfc=vlfc, valid_padj=vpadj,
        rescued_by_validation=rescued,
    ))

res = pd.DataFrame(rows)
ensure_parent(snakemake.output.recall)  # noqa: F821
res.to_csv(snakemake.output.recall, sep="\t", index=False)  # noqa: F821

# ---- summary -----------------------------------------------------------------
n = len(res)
n_selected = int((res["found_status"] == "selected").sum())
n_surfaced = int(res["found_status"].isin(["selected", "scored_not_selected"]).sum())
n_sig_primary = int(res["found_status"].str.startswith("DE_significant").sum()) + n_surfaced
n_rescued = int(res["rescued_by_validation"].sum())
cat_counts = res["found_status"].value_counts()

summary = pd.DataFrame({
    "metric": ["n_known_positives", "recall_selected", "recall_selected_frac",
               "surfaced_as_candidate", "detected_significant_in_primary",
               "rescued_by_validation_only"],
    "value": [n, n_selected, round(n_selected / n, 3) if n else 0,
              n_surfaced, n_sig_primary, n_rescued],
})
ensure_parent(snakemake.output.summary)  # noqa: F821
summary.to_csv(snakemake.output.summary, sep="\t", index=False)  # noqa: F821

# ---- figure: outcome breakdown -----------------------------------------------
ensure_parent(snakemake.output.fig)  # noqa: F821
order = ["selected", "scored_not_selected", "DE_significant_off_list",
         "DE_significant_wrong_dir", "DE_tested_nonsignificant", "DE_padj_NA",
         "detected_filtered_pre_DE", "absent_primary"]
counts = [int(cat_counts.get(k, 0)) for k in order]
colors = ["#1a9850", "#91cf60", "#d9ef8b", "#fee08b", "#fdae61", "#f46d43",
          "#d73027", "#a50026"]
fig, ax = plt.subplots(figsize=(7.5, 4.2))
ax.barh(order[::-1], counts[::-1], color=colors[::-1])
ax.set_xlabel("number of known positives")
ax.set_title(f"Known-positive recall: {n_selected}/{n} selected "
             f"({100 * n_selected / n:.0f}%); {n_rescued} rescued by validation only")
for i, v in enumerate(counts[::-1]):
    if v:
        ax.text(v + 0.05, i, str(v), va="center", fontsize=8)
fig.tight_layout()
fig.savefig(snakemake.output.fig, dpi=120)  # noqa: F821

rescued_genes = res.loc[res["rescued_by_validation"], "gene"].tolist()
log.info(f"known-positive recall: {n_selected}/{n} selected ({100*n_selected/n:.0f}%); "
         f"{n_surfaced} surfaced; {n_sig_primary} significant in primary; "
         f"{n_rescued} rescued-by-validation-only ({rescued_genes}); "
         f"breakdown={dict(cat_counts)}")

"""Stage 07b — condition-differential communication + downstream-target corroboration (Q7).

Overinterpretation guards implemented here:
  1. specificity filter on the LIANA consensus (keep specific interactions only);
  2. condition-differential ranking — keep interactions enriched in the fibrotic condition,
     not those present everywhere;
  3. downstream-target corroboration (NicheNet-style): a predicted ligand->receiver link is
     only "corroborated" if the ligand's canonical target genes are concordantly UP in the
     receiver compartment's donor-aware DE. This requires a transcriptional response, guarding
     against expression-only false positives.
(Full NicheNet is available via workflow/scripts/R/nichenet_corroborate.R for users who want
the licensed prior model; this offline corroboration keeps the default pipeline self-contained.)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import ensure_parent, get_logger  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
p = snakemake.params  # noqa: F821
spec_cut = float(p.specificity_cutoff)  # noqa: F821
min_targets = int(p.min_targets)  # noqa: F821

# canonical downstream targets of key fibrotic ligands (receiver transcriptional response)
LIGAND_TARGETS = {
    "TGFB1": ["SERPINE1", "COL1A1", "COL3A1", "CTGF", "CCN2", "TIMP1", "ACTA2", "TAGLN"],
    "PDGFB": ["COL1A1", "ACTA2", "TIMP1", "MKI67", "PDGFRB"],
    "PDGFD": ["COL1A1", "ACTA2", "PDGFRB"],
    "SPP1": ["COL1A1", "TIMP1", "ACTA2", "TGFB1"],
    "TNFSF12": ["CCL2", "MMP9", "COL1A1", "TIMP1"],
    "JAG1": ["HEY1", "HES1", "NOTCH3"],
    "DLL4": ["HEY1", "HES1", "NOTCH3"],
    "CCL2": ["CCR2", "IL1B", "CXCL8"],
    "IL1B": ["CCL2", "CXCL8", "MMP9", "NFKB1"],
}

liana = pd.read_csv(snakemake.input.liana, sep="\t")  # noqa: F821
de_files = list(snakemake.input.de)  # noqa: F821

ensure_parent(snakemake.output.diff)  # noqa: F821
ensure_parent(snakemake.output.nichenet)  # noqa: F821

diff_cols = ["source", "target", "ligand_complex", "receptor_complex",
             "fibrotic_condition", "delta_specificity", "enriched_in_fibrosis"]
nn_cols = ["ligand", "receiver_compartment", "n_targets_up", "frac_targets_up",
           "corroborated", "targets_up"]

if liana.empty or "specificity_rank" not in liana:
    pd.DataFrame(columns=diff_cols).to_csv(snakemake.output.diff, sep="\t", index=False)  # noqa: F821
    pd.DataFrame(columns=nn_cols).to_csv(snakemake.output.nichenet, sep="\t", index=False)  # noqa: F821
    log.warning("No LIANA input; wrote empty differential + corroboration tables")
    sys.exit(0)

# guard 1: keep specific interactions
liana = liana[liana["specificity_rank"] <= max(spec_cut, liana["specificity_rank"].min())]

# guard 2: condition-differential. specificity_rank: lower = more specific -> use -rank as score
liana["score"] = -np.log10(liana["specificity_rank"].clip(lower=1e-6))
conds = [c for c in liana["condition"].unique() if c != "all"]
fibrotic = next((c for c in conds if str(c).lower() in
                 ("cirrhotic", "high", "mash", "fibrotic")), conds[0] if conds else None)

key = ["source", "target", "ligand_complex", "receptor_complex"]
if fibrotic is not None and liana["condition"].nunique() >= 2:
    piv = liana.pivot_table(index=key, columns="condition", values="score", aggfunc="mean").fillna(0)
    other = [c for c in piv.columns if c != fibrotic]
    piv["delta_specificity"] = piv[fibrotic] - piv[other].mean(axis=1)
    diff = piv.reset_index()
    diff["fibrotic_condition"] = fibrotic
    diff["enriched_in_fibrosis"] = diff["delta_specificity"] > 0
else:
    diff = liana.groupby(key, as_index=False)["score"].mean()
    diff["fibrotic_condition"] = fibrotic if fibrotic else "NA"
    diff["delta_specificity"] = diff["score"]
    diff["enriched_in_fibrosis"] = True

diff = diff.sort_values("delta_specificity", ascending=False)
diff[[c for c in diff_cols if c in diff.columns]].to_csv(snakemake.output.diff, sep="\t", index=False)  # noqa: F821

# guard 3: downstream-target corroboration against receiver DE
# load per-compartment DE (up genes in fibrosis), keyed by compartment from the file path
de_up = {}
for f in de_files:
    comp = os.path.basename(os.path.dirname(f))
    d = pd.read_csv(f, sep="\t")
    if d.empty or "log2FoldChange" not in d:
        de_up[comp] = set()
        continue
    sig = d[(d["log2FoldChange"] > 0.25) & (d.get("padj", 1).fillna(1) < 0.25)]
    de_up[comp] = set(sig["gene"].astype(str))

rows = []
fib_diff = diff[diff.get("enriched_in_fibrosis", True)] if "enriched_in_fibrosis" in diff else diff
for _, r in fib_diff.iterrows():
    ligand = str(r["ligand_complex"]).split("_")[0]
    receiver = str(r["target"])
    targets = LIGAND_TARGETS.get(ligand)
    if not targets:
        continue
    up = de_up.get(receiver, set())
    hit = sorted(set(targets) & up)
    rows.append(dict(
        ligand=ligand, receiver_compartment=receiver, n_targets_up=len(hit),
        frac_targets_up=round(len(hit) / len(targets), 3),
        corroborated=len(hit) >= min_targets, targets_up=",".join(hit),
    ))
nn = pd.DataFrame(rows, columns=nn_cols).drop_duplicates().sort_values(
    "n_targets_up", ascending=False) if rows else pd.DataFrame(columns=nn_cols)
nn.to_csv(snakemake.output.nichenet, sep="\t", index=False)  # noqa: F821
log.info(f"differential CCC: {len(diff)} interactions; corroborated ligand-receiver links: "
         f"{int(nn['corroborated'].sum()) if not nn.empty else 0}")

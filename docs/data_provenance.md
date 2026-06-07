# Data provenance & dataset choice

## Primary dataset — GSE136103
- **Study:** Ramachandran et al. 2019, *Nature*, "Resolving the fibrotic niche of human liver
  cirrhosis at single-cell level." doi:10.1038/s41586-019-1631-3.
- **Content:** human liver scRNA-seq (10x 3′ v2), healthy vs cirrhotic, CD45+/CD45− FACS-sorted;
  defines the disease states this project validates — scar-associated macrophages (TREM2/CD9/
  SPP1), scar-associated endothelium (ACKR1/PLVAP), fibrogenic mesenchyme/myofibroblasts
  (PDGFRA/PDGFRB, collagen).
- **Access:** GEO `GSE136103` (per-sample 10x MEX in `GSE136103_RAW.tar`); raw reads at SRA
  `SRP218975`; a processed Seurat object at Edinburgh DataShare (handle `10283/3433`, CC-BY 4.0).
- **Key caveat:** GEO's public sample metadata is essentially **binary (healthy / cirrhotic)** —
  there is **no fibrosis staging** — so staged contrasts (F2+) must come from the validation arm.

## Validation dataset — GSE244832 (chosen)
- **Study:** human hepatic stellate cell / MASLD–MASH target-discovery dataset (snRNA-seq +
  snATAC), 18 livers spanning Normal / MASL / MASH **F2–F4** (Kim/Brenner/Kisseleva 2024,
  *J. Hepatology*). Open on GEO.

### How it was obtained (a real Q1 metadata-curation case)
- GEO lists 36 samples = **18 snATAC + 18 snRNA** (the snRNA GSMs are `GSM7830559–7830576`).
  Only the 18 snRNA are used for expression validation.
- **The GEO per-sample metadata contains no disease/stage labels** — each sample's
  characteristics are literally `tissue: liver` with an internal code (`JB288`). The graded
  fibrosis grouping is **not recoverable from GEO's machine-readable metadata** — a textbook
  example of the "incomplete clinical metadata" challenge (screening Q1).
- The labels live in the processed bundle `GSE244832_hLIVER_processed_files.tar.gz` (one
  combined matrix `128,223 cells × 27,200 genes` + per-cell metadata with a `condition` column),
  which recovers the paper's cohort exactly: **5 NORMAL + 4 NAFL/MASL + 9 NASH/MASH**.
- `workflow/scripts/py/prep_gse244832.py` splits that bundle into 18 per-sample 10x triplet
  dirs and writes the harmonized validation samplesheet, mapping `NORMAL→0, NAFL→1, NASH→3`
  onto the ordinal axis → a balanced **9-high / 9-low F2+ contrast**.

### Demonstration-run QC note
The committed `config.yaml` documents decontX (ambient) and Scrublet/scDblFinder (doublets) as
the principled defaults, but the **full real run uses `none` for both**: GSE136103 is
CellRanger-cell-filtered and GSE244832 is author-processed (so ambient/doublet removal is
low-value), and the decontX R path (zellkonverter→basilisk) is prohibitively heavy. Both methods
remain wired in and config-selectable.
- **Why this one (the decision):**
  - **Independence.** It shares no samples with the primary, so it is a genuine out-of-sample
    test (see the SCP2154 note below for why this matters).
  - **Adds a fibrosis-severity gradient** (F0→F4) that the binary primary dataset lacks — letting
    us test whether candidate biomarkers *track fibrosis progression*, the core translational claim.
  - **Orthogonal etiology** (metabolic MASH vs the primary's mixed/cirrhotic) and **orthogonal
    modality** (snRNA vs scRNA) — a biomarker that replicates across both is far more credible.
  - **Deep coverage of the central fibrogenic compartment** (hepatic stellate / myofibroblast),
    the richest source of cell-type-specific fibrosis targets; the multi-omic (ATAC) layer also
    supports TF/mechanism corroboration.
  - **Open and programmatically accessible** (GEO), so the pipeline can fetch it reproducibly.
  - *Accepted trade-off:* shallower macrophage/endothelial coverage — acceptable because the
    primary already profiles those compartments deeply via CD45+ sorting.

## Candidates considered and not chosen
- **GSE207310 (SMOC2 / NAFLD severity).** Uniquely offers **plasma-protein** biomarker validation
  (SMOC2, AUROC ≈0.88) and a clean NAFLD severity gradient — but its **primary modality is bulk
  RNA-seq**, so it is weak for *cell-type-specific* validation. Retained as an optional **targeted
  secondary** check if a top candidate is secreted (e.g., SMOC2) and we want plasma-level
  confirmation.
- **Broad SCP2154 (macrophage atlas).** Broadest compartment coverage, **but it is a meta-atlas
  that integrates GSE136103 itself** → using it to "validate" GSE136103 findings would be
  **circular (data leakage)**, and download is login-gated. Used only as a published-signature
  cross-reference for the scar-associated-macrophage program, never as the validation dataset.

## How reproducibility handles the sc/sn + etiology gap
Cross-dataset comparison (`workflow/scripts/py/crossdataset_repro.py`) is done at the level of DE
**direction/rank** (sign-concordance, Spearman of the Wald statistic), which is robust to absolute
expression differences between scRNA and snRNA and between etiologies; reproducibility is one
**weighted component** of the biomarker score, not a hard filter.

## Annotation resource snapshots (offline, deterministic)
The biomarker score (`workflow/scripts/py/build_features.py`) annotates each gene against dated
reference snapshots committed under `resources/`, so a clean checkout scores identically with no
network access. Each was downloaded once (2026-06-07) and frozen:

| Snapshot | Source | Genes | Used for |
|---|---|---|---|
| `secretome_hpa_2019.csv` | Human Protein Atlas, **Predicted secreted proteins** protein class (Uhlen et al. 2019, *Sci. Signal.*, "The human secretome") | 1,902 | accessibility = `secreted` |
| `surfaceome_surfy_2018.csv` | **SURFY** in-silico human surfaceome (Bausch-Fluck et al. 2018, *PNAS*), sheet "in silico surfaceome only" | 2,799 | accessibility = `surface` |
| `druggability_snapshot.tsv` | Open Targets tractability + DGIdb (cached) | (panel) | druggability component |

**Accessibility precedence is surface > secreted > unknown.** A gene in the surfaceome is called
`surface` even if it is also in the secretome, because a membrane protein that is also predicted-
secreted is a **shed receptor** (e.g. soluble TREM2) whose structural localization is the cell
surface. A gene in neither reference is labeled **`unknown`**, not `intracellular`: absence from
these (non-exhaustive) sets is not positive evidence of intracellular localization, so the pipeline
does not assert a localization it has not verified. This replaced an earlier pair of ~30-gene
hand-curated lists whose unknown-defaults mislabeled bona fide secreted proteins (e.g. MMP7, SPP2,
A2M) as intracellular.

# Methods & rationale

Per-stage method choices and why, with citations. This complements `written_answers.md`
(which frames the same choices against the screening questions) and is the narrative behind
the `config/` defaults.

## Ingestion & harmonization (Q1)
Per-sample 10x MEX → AnnData (`scanpy.read_10x_mtx`); raw counts stashed in `layers["counts"]`.
Heterogeneous fibrosis labels are mapped to a common ordinal axis (0–4) with the original label
retained; disease terms follow the 2023 MASLD/MASH nomenclature (Rinella et al. 2023,
*Hepatology*, PMID 37363821). Metadata follow CZ CELLxGENE / HCA tier-1 conventions.

## QC (Q2)
Heumos et al. 2023 (*Nat Rev Genet*; sc-best-practices) consensus: **adaptive median±n·MAD**
thresholds on log counts/genes; **modality- and liver-aware** mitochondrial *upper* bound only
(hepatocytes are mito-high; snRNA is mito-low); ambient RNA via **decontX** for filtered matrices
(CellBender/SoupX need raw droplets); doublets via **scDblFinder** (Germain 2022, best AUPRC) or
Scrublet. Stressed/diseased cells are **flagged, not filtered** (`keep_stressed`).

## Integration (Q3)
Batch-correct **over donor/sample, not condition** — Harmony (Korsunsky 2019) by default, scVI/
scANVI (scvi-tools) when GPU is enabled. Quality is judged with **scIB** metrics (Luecken 2022,
*Nat Methods*, PMID 34949812), separating batch-removal from bio-conservation; a **fibrosis-signal
guard** fails the run if condition becomes inseparable post-integration.

## Annotation (Q4)
Leiden clustering + curated marker scoring + optional CellTypist cross-check, rolled up to
compartments. The mesenchymal cluster is disambiguated (activated HSC vs portal fibroblast vs
myofibroblast vs VSMC) with lineage-discriminating panels (Ramachandran 2019; Dobie 2019) reported
in `compartment_validation.tsv`.

## Differential expression (Q5)
**Donor-aware pseudobulk** is the crux: aggregate counts to donor level, then DESeq2/edgeR/muscat
with the donor as replication unit — avoiding pseudoreplication-inflated false discoveries (Squair
et al. 2021, PMID 34594038; Zimmerman 2021; Murphy & Skene 2022; Crowell muscat 2020). Min cells/
donor and min donors/group gates; covariate-adjusted design; ashr-shrunken fold-changes for
ranking with the Wald statistic retained.

## Pathway & mechanism (Q7)
decoupler PROGENy (pathway activity) + CollecTRI (TF activity); preranked GSEA (gseapy) on the DE
statistic against curated fibrosis sets (TGF-β, PDGF, Notch, Hedgehog, YAP/TAZ, Wnt, TWEAK/Fn14)
and MSigDB Hallmark; AUCell per-cell signature scores.

## Cell-cell communication (Q7)
**LIANA** consensus ligand-receptor (Dimitrov 2022/2024), run per condition. Overinterpretation
guards: expression-proportion filter, specificity ranking, **condition-differential** prioritization,
and **downstream-target corroboration** (a link is kept only if the ligand's known targets are
concordantly induced in the receiver's DE — NicheNet logic, Browaeys 2019). Full NicheNet is
available in `scripts/R/nichenet_corroborate.R`.

## Biomarker prioritization (Q6)
Transparent composite over six normalized axes — DE evidence, cell-type specificity (**tau** index;
Kryuchkova-Mostacci & Robinson-Rechavi 2017), cross-dataset reproducibility, druggability (Open
Targets/DGIdb), biomarker accessibility (secretome/surfaceome), and an explainable **XGBoost+SHAP**
signal trained on donor-level pseudobulk with **donor-grouped CV**. Top-N per required compartment
then fill to top-N overall; weights are config-driven and a sensitivity sweep is reported.

## Reproducibility & engineering (Q8)
Snakemake + Scanpy, config-driven (JSON-schema-validated), per-rule pinned conda envs (built fast
via micromamba), Docker→Singularity containers, local + SLURM profiles, central seeds, provenance
dump, Quarto reporting, and CI (lint → dry-run → tiny end-to-end on the committed fixture).

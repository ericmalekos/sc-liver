# Written Screening Questions — Responses

Each answer reflects current best practice (citations
inline) and points to where the practical pipeline in this repository implements it. Gene
symbols are human (HGNC). Abbreviations: HSC = hepatic stellate cell; SAMac = scar-associated
macrophage; DE = differential expression; sc/sn = single-cell / single-nucleus.

---

## Q1 — Dataset curation and fibrosis-stage harmonization

**The core problem** is that "fibrosis stage" is recorded on incompatible scales across studies
(METAVIR F0–F4, Ishak 0–6, NASH-CRN/Kleiner 0–4, or a coarse cirrhosis/non-cirrhosis binary),
disease nomenclature shifted in 2023 (NAFLD→**MASLD**, NASH→**MASH**; Rinella et al. 2023,
*Hepatology*, PMID 37363821), and clinical metadata are frequently incomplete. Naively pooling
these labels produces a confounded outcome variable.

**Approach.** (1) **Map every study onto a common ordinal axis (0–4)** while *retaining the
original label* for provenance — periportal/F2 ≈ axis 2, bridging/F3 ≈ axis 3, cirrhosis/F4 ≈
axis 4 (crosswalk in `docs/reference_liver_cell_types_fibrosis.md`). (2) **Standardize disease
terms to MASLD/MASH** and record etiology (metabolic vs ALD vs viral vs cholestatic) as a
separate covariate, because etiology, not just stage, drives composition. (3) **Adopt a
structured metadata schema** — the CZ CELLxGENE schema and HCA tier-1 fields — so donor, assay,
tissue, and disease are captured consistently and the object is interoperable. (4) **Handle
missing metadata explicitly**: treat stage as *ordinal* (and analyze with a binarized "F2+/
significant fibrosis" contrast when fine stages are unavailable), flag imputed values, and run
**sensitivity analyses** dropping incompletely-annotated donors. (5) **Validate harmonization
biologically, not just syntactically**: confirm that canonical signatures recover the expected
gradient (e.g., collagen/ACTA2 and the SAMac program increase monotonically with axis), that
label transfer between datasets is concordant, and that technical covariates (assay, chemistry,
sc vs sn) are not confounded with stage.

**In this pipeline.** `workflow/scripts/py/make_samplesheet.py` parses GEO characteristics and
maps labels to the axis via `config.metadata.fibrosis_mapping`; `harmonize_metadata.py` emits a
per-dataset harmonized samplesheet with the axis + a binary F2+ split + cell counts; and the
**cross-dataset reproducibility** step (`crossdataset_repro.py`) is the biological validation —
it checks that DE *direction* agrees between the primary (GSE136103) and validation (GSE244832)
cohorts despite their different scales, etiologies, and modalities. A key honest caveat we
surface: **GSE136103's public metadata is only healthy-vs-cirrhotic (binary)**, so staged
contrasts come from the graded validation cohort. We hit this problem for real with the
validation set: **GSE244832's GEO sample metadata carries no disease/stage labels at all**
(just `tissue: liver` + an internal code), so the cohort grouping had to be recovered from the
authors' processed bundle (`condition` column → 5 NORMAL / 4 MASL / 9 MASH) and mapped onto the
axis — exactly the kind of "incomplete metadata" curation this question describes (see
`docs/data_provenance.md` and `workflow/scripts/py/prep_gse244832.py`).

---

## Q2 — QC and preprocessing for liver sc/snRNA-seq

**The tension** is removing genuine artifacts (empty/dying droplets, ambient RNA, doublets,
batch effects) **without** deleting the biologically interesting stressed/activated/diseased
cells, which can superficially resemble low-quality cells. The current consensus reference is
the single-cell best-practices book (Heumos et al. 2023, *Nat Rev Genet*).

**Steps and rationale.**
- **Ambient RNA.** CellBender (Fleming et al. 2023) and SoupX need the *raw, unfiltered* droplet
  matrix and are ideal when available; for the **filtered** matrices GEO typically ships,
  **decontX** (Yang et al., celda) estimates contamination from the expression itself and is the
  right default — implemented in `workflow/scripts/R/qc_ambient.R`.
- **Doublets.** **scDblFinder** has the best benchmarked AUPRC (Germain et al. 2022); Scrublet is
  a solid Python option. We support both (`scdblfinder.R` / `qc_doublets.py`).
- **Adaptive thresholds, not fixed cutoffs.** Filter on **median ± n·MAD** (n=3–5) of
  log-transformed counts and detected genes rather than hard universal cutoffs (Heumos 2023);
  fixed thresholds discard valid cells in low-RNA compartments.
- **Liver- and modality-aware mitochondrial QC.** Hepatocytes are legitimately **mito-high**, so
  an aggressive global mito cutoff erases them; **snRNA-seq is mito-low**, so the ceiling must be
  modality-specific. We use only an **upper** mito bound (`config.qc.mito_pct_max` = {scRNA:8,
  snRNA:5}); a probabilistic alternative is **miQC** (Hippen et al. 2021).
- **Protecting stressed/diseased cells.** We **do not filter on stress/IEG or dissociation
  signatures** (`keep_stressed: true`), because activated HSCs, SAMacs, and dying parenchyma
  carry real stress programs; instead we *flag* them. snRNA-seq is also preferable for fibrotic/
  cirrhotic tissue because it avoids dissociation-induced artifacts in fragile parenchyma.
- **Batch artifacts** are handled at integration (Q3), not by deletion.

**In this pipeline.** Stage 02 (`qc_ambient → qc_doublets → qc_filter`) applies adaptive MAD
bounds + a modality-aware mito ceiling, drops flagged doublets, keeps stressed cells, and writes
a per-sample metrics JSON + before/after figures aggregated by `qc_summary.py`.

---

## Q3 — Integration without erasing fibrosis biology

**The risk** is that batch correction, applied too aggressively or over the wrong variable, also
removes the disease signal we are trying to study. The governing principle (Heumos 2023; Luecken
et al. 2022, *Nat Methods*, scIB benchmark, PMID 34949812) is: **correct over the technical/
nuisance variable (donor, sample, assay) — never over the biological condition.**

**Approach.** (1) **Integrate over donor/sample**, with disease condition left as a free
biological axis (`config.integration.batch_key = sample_id`). (2) **Choose the method to the
problem**: Harmony (Korsunsky 2019) is fast and a strong CPU default; **scVI/scANVI**
(scvi-tools) handle complex multi-batch settings and, in the *semi-supervised* scANVI form, can
preserve labelled biology. (3) **Benchmark with scIB metrics** that explicitly separate
*batch-removal* from *bio-conservation*, and pick the operating point that mixes batches **while
keeping bio-conservation high** rather than maximizing mixing. (4) **Verify the disease signal
survived**: the condition should remain separable in the integrated embedding, fibrosis
signatures should still track the stage axis, and downstream pseudobulk DE results should be
stable pre/post integration.

**In this pipeline.** `integrate.py` batch-corrects over `sample_id` (Harmony default, scVI/
scANVI when `gpu.enabled`) and **retains the full gene set** (PCA on HVGs only) so DE later uses
all genes. `scib_benchmark.py` reports batch-mixing vs. condition-separation metrics and enforces
a **fibrosis-signal guard** that fails the run if condition is no longer separable post-
integration — a concrete safeguard against over-correction.

---

## Q4 — Cell-type annotation and validation in fibrotic liver

A cluster expressing COL1A1, COL3A1, ACTA2, TAGLN, PDGFRB, LUM, DCN and enriched in F3/F4 is
**collagen-producing mesenchyme** — but "fibroblast" is too coarse. These markers are shared
across activated HSCs, portal fibroblasts, and terminal myofibroblasts, which is exactly why
automated labels must be validated, not trusted.

**How to disambiguate.**
- **Subcluster the mesenchyme** and score lineage-discriminating panels (not the shared pan-
  fibroblast genes): **HSC origin** → RGS5, DES, LRAT, RBP1, NGFR + retinoid/quiescence remnants;
  **portal fibroblast** → THY1(CD90), CD34, ELN, FBLN1, MSLN, SLIT2, GPC3, PDGFRA; **terminal
  myofibroblast** → ACTA2-high, TAGLN, TIMP1, TNC, POSTN with loss of quiescence markers;
  **vascular smooth muscle** → MYH11, NOTCH3, PLN (perivascular, not scar). (Markers curated in
  `config/markers_liver.yaml`; biology in the reference doc, grounded in Ramachandran 2019,
  Dobie et al. 2019.)
- **Use trajectory/ordering** (diffusion/PAGA) to test whether the cluster lies on a quiescent-
  HSC → activated-myofibroblast continuum, which favors an HSC origin over portal fibroblast.
- **Use spatial/zonation priors**: HSCs are sinusoidal/pericentral; portal fibroblasts are
  periportal. Spatial transcriptomics or IHC localizes the population in the scar.
- **Cross-reference reference atlases / label transfer** (CellTypist, Azimuth liver, SingleR) and
  require agreement across methods.
- **Acknowledge mixtures**: in advanced fibrosis the population is often a *continuum/mixed
  state*, not a discrete type — report it as such.

**In this pipeline.** `cluster_annotate.py` assigns clusters by marker scoring; `celltypist_-
annotate.py` finalizes types, rolls them to compartments, optionally cross-checks with CellTypist,
and writes `compartment_validation.tsv` containing, per cluster, the cell-type scores, the
**disease-state** scores (activated-HSC, SAMac, scar-endothelium), the **stromal discriminators**
(HSC vs portal-fibroblast vs myofibroblast vs VSMC), and the fraction of cells from F3+ samples —
i.e., the explicit evidence to make this call.

---

## Q5 — Donor-aware differential expression and biomarker discovery

**Why cell-level DE is dangerous.** Treating individual cells as independent replicates is
**pseudoreplication**: cells from one donor are correlated, so the effective sample size is the
number of *donors*, not cells. Cell-level tests (e.g., Wilcoxon on tens of thousands of cells)
therefore produce drastically inflated false positives and frequently report *donor* identity or
a single outlier donor as if it were disease biology. This is well documented: Squair et al. 2021
(*Nat Commun*, PMID 34594038), Zimmerman et al. 2021, and Murphy & Skene 2022 all show
**pseudobulk methods control type-I error far better** than naive single-cell tests.

**The correct strategy — donor-aware pseudobulk.** Aggregate counts **to the donor × cell-type
level** (sum), then apply a bulk DE engine (**DESeq2 / edgeR / limma-voom**, or **muscat**;
Crowell et al. 2020) with the **donor as the replication unit**. Specifics that matter for an
F2+ contrast in macrophages and endothelium: require a **minimum number of cells per donor** (so
a pseudobulk profile is stable) and a **minimum number of donors per group**; **adjust for
covariates** (sex, sort gate, etiology, batch); shrink fold-changes (apeglm/ashr) for ranking;
and, where a mixed-model is preferred over aggregation, use **NEBULA** (donor as random effect).
Report effect sizes and CIs, not just p-values.

**In this pipeline.** `pseudobulk_aggregate.py` sums counts to donor level per compartment with a
`min_cells_per_donor` gate; `deseq2_pseudobulk.R` runs DESeq2 (default) / edgeR / muscat with a
covariate-aware design, **writing an explicit "insufficient donors" status rather than a spurious
result** when a group is underpowered. This donor→DE boundary is the scientific crux of the whole
project.

---

## Q6 — AI/ML-based biomarker prioritization

With ~300 candidates across HSC, macrophage, endothelial, and cholangiocyte compartments, the
goal is to convert a long DE list into a short, *defensible, translationally-ranked* list — and
to do so **explainably**, since a black-box ranking is not actionable for target/biomarker
decisions.

**Approach — a transparent composite score plus an explainable ML signal.** Score each
**(gene, compartment)** candidate on orthogonal axes, each mapped to a translational criterion:
1. **DE evidence** (effect size × significance from the donor-aware test);
2. **Cell-type specificity** — the **tau index** (most robust specificity metric; Kryuchkova-
   Mostacci & Robinson-Rechavi 2017), so we favor markers that are *specific* to a compartment;
3. **Cross-dataset reproducibility** — sign/rank concordance between primary and validation,
   robust to the sc/sn modality gap;
4. **Druggability/tractability** — Open Targets tractability + DGIdb interactions;
5. **Biomarker accessibility** — secreted (best for circulating diagnostics) > cell-surface
   (targetable) > intracellular, via surfaceome/secretome membership;
6. **An explainable ML signal** — an **XGBoost** (or random-forest) classifier of fibrotic-vs-
   healthy trained on **donor-level pseudobulk with donor-grouped cross-validation** (no cell or
   donor leakage), with **SHAP** values giving per-gene importance.

The ML layer is **additive and audited**, not a black box driving the ranking: a gene with high
SHAP but no DE support is surfaced separately as an "ML-only hypothesis." We report a **weight-
sensitivity analysis** so the ranking's robustness is visible, and select **top-N per required
compartment first** (so each disease compartment is represented) before filling to the overall
top-20.

**In this pipeline.** `build_features.py` (features), `ml_classifier_shap.py` (donor-grouped
XGBoost + SHAP), and `compute_score.py` (within-compartment normalized composite, configurable
weights, per-candidate rationale) produce `results/08_score/candidate_scores.tsv`.

---

## Q7 — Cell–cell interaction and pathway mechanism discovery

**Goal.** Recover the fibrotic crosstalk — SAMac ↔ activated HSC ↔ remodeling endothelium — while
not over-interpreting ligand-receptor (LR) predictions, which are notorious for plausible-looking
false positives.

**Method.** Use a **consensus** LR method — **LIANA/LIANA+** (Dimitrov et al. 2022, 2024) — which
aggregates CellPhoneDB, CellChat, NATMI, Connectome, etc., rather than trusting one tool. Layer
**pathway and TF activity** with **decoupler** (PROGENy for pathways, CollecTRI for TFs) and
**GSEA** on the DE ranking, plus per-cell signature scores (AUCell/UCell), focusing on the
fibrosis pathways TGF-β, PDGF, Notch, Hedgehog, YAP/TAZ, Wnt, TWEAK/Fn14.

**Guarding against over-interpretation** (the crux of the question):
1. **Expression filters** — require the ligand and receptor to be expressed in a minimum fraction
   of cells in each population (drops noise-driven hits).
2. **Specificity / permutation nulls** — keep interactions that are *specific* to a population
   pair, not ubiquitously "significant".
3. **Condition-differential analysis** — prioritize interactions **enriched in fibrosis**, not
   those present everywhere.
4. **Require a downstream transcriptional response** — a predicted ligand→receiver link is only
   credible if the ligand's known target genes are concordantly induced in the receiver
   (NicheNet's logic; Browaeys et al. 2019). LR co-expression alone is *not* evidence of signaling.
5. **Orthogonal validation** — spatial transcriptomics/IHC for co-localization, or perturbation/
   recombinant-ligand assays, before any causal claim.

Expected, literature-supported hits to look for: SPP1/PDGFB/TNFSF12 (SAMac) → CD44/PDGFRB/
TNFRSF12A (HSC); JAG1/DLL4 (endothelium) → NOTCH3 (HSC); TGFB1 → TGFBR (HSC) (Ramachandran 2019).

**In this pipeline.** `pathway_decoupler.py` + `gsea.py` (mechanism) and `ccc_liana.py` (consensus
LR, per condition) → `ccc_differential.py`, which applies the specificity filter, computes
condition-differential communication, and **corroborates each link against the receiver's
donor-aware DE** (offline NicheNet-style); full NicheNet is available in `scripts/R/nichenet_-
corroborate.R`.

---

## Q8 — Reproducible pipeline, GitHub workflow, and delivery plan (12–16 weeks)

**This repository is the answer.** It is a modular, config-driven **Snakemake + Scanpy** pipeline
built with nf-core conventions: one rule file per stage (`workflow/rules/00…10`), one script per
task (`workflow/scripts/{py,R}`), all parameters in `config/` (validated against a JSON schema),
**per-rule pinned conda environments** (`envs/`), **Docker→Singularity containers** for the
cluster, **local and SLURM profiles** (`profiles/`), a committed **downsampled test fixture** +
**CI** (lint → dry-run → tiny end-to-end), central **seeds** and a **provenance** dump (resolved
config, git SHA, tool versions) for determinism.

**Repository structure / tools.** AnnData (`.h5ad`) is the currency between Python stages; the R
boundary (DESeq2, scDblFinder, decontX) is plain TSV/RDS. Reporting is **Quarto** (HTML + PDF).
Reproducibility is enforced by containers + pinned envs + seeds + schema-validated config; CI
proves the *whole* pipeline runs on the fixture, not just that it parses.

**Milestones, quality gates, deliverables** (one biweekly block per screening-question theme):

| Weeks | Milestone | Quality gate | Deliverable |
|---|---|---|---|
| 1–2 | Scaffold + reproducibility spine + reference doc | CI green (lint + dry-run); fixture builds | repo skeleton, README, CI |
| 3–4 | Ingestion + metadata harmonization (Q1) | md5-verified data; harmonized axis reviewed | dataset/metadata summary |
| 5–6 | QC/preprocessing (Q2) | stressed/diseased cells retained; QC figures sane | **QC summary** |
| 7–8 | Integration + annotation (Q3, Q4) | scIB + fibrosis-signal guard pass; 3 compartments validated | **annotation figures** |
| 9–10 | Donor-aware DE + pathway (Q5, Q7) | ≥ min donors/group; fibrosis pathways recovered | **DE + pathway results** |
| 11–12 | Cell-cell communication (Q7) | LR pairs pass filters + downstream corroboration | CCC tables/figures |
| 13–14 | ML prioritization + cross-dataset (Q6, Q1) | donor-grouped CV; weights stable; arms concordant | **ranked 10–20 table** |
| 15–16 | Full SLURM run, report, hardening (Q8) | reproducible from clean checkout; report renders | **report + executive summary** |

**Final deliverables**: GitHub repo + README; QC summary; cell-annotation figures; DE + pathway
results; ranked biomarker/target table; 1–2 page executive summary; and these written answers.

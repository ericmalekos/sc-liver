# Risk register

| # | Risk | Impact | Mitigation (where) |
|---|---|---|---|
| 1 | **GSE136103 has no fibrosis staging** in public metadata (binary healthy/cirrhotic) | Can't make staged (F2+) contrasts on the primary | Primary contrast = cirrhotic-vs-healthy; staged F2+ logic lives in the validation arm; stated as the headline biological caveat (`data_provenance.md`, report) |
| 2 | **Cell-level DE pseudoreplication** | Inflated false positives, donor effects mistaken for disease | Donor-aware pseudobulk (`pseudobulk_aggregate.py` → `deseq2_pseudobulk.R`); min cells/donor + min donors/group gates; covariate-adjusted design (Squair 2021) |
| 3 | **Integration erases fibrosis biology** | Disease signal lost to over-correction | Batch_key = donor/sample, never condition; `scib_benchmark.py` **fibrosis-signal guard** fails the run if condition becomes inseparable |
| 4 | **sc vs sn + etiology gap** between primary and validation | Naive expression comparison misleading | Compare DE **direction/rank** (sign-concordance), snRNA-aware mito QC; reproducibility is a weighted component, not a hard filter |
| 5 | **Removing real stressed/diseased cells in QC** | Loss of activated HSCs / SAMacs / dying parenchyma | Adaptive MAD (not fixed cutoffs); modality- & liver-aware mito ceiling (upper-only); `keep_stressed: true` |
| 6 | **CellBender/SoupX need raw droplets**; GEO ships filtered matrices | Ambient method may not apply | Default ambient = **decontX** (works on filtered data); CellBender path documented as GPU/raw-droplet only |
| 7 | **Ligand-receptor over-interpretation** | Plausible-but-false signaling claims | Consensus LIANA + expression/specificity filters + condition-differential + **downstream-target corroboration** (`ccc_differential.py`) |
| 8 | **ML on few donors overfits** | Unstable SHAP / inflated CV | Donor-grouped CV (no leakage), modest model depth, SHAP as an *audited additive* component (gate, not black-box driver); weight-sensitivity reported |
| 9 | **Live Open Targets/DGIdb API flakiness/non-determinism** | Irreproducible druggability scores | Ship a **dated cached snapshot** (`resources/druggability_snapshot.tsv`); regeneration documented |
| 10 | **OmniPath (PROGENy/CollecTRI) unreachable offline** | Pathway/TF activity missing | `pathway_decoupler.py` wraps fetches and falls back to local-GMT AUCell scoring |
| 11 | **Quarto/TinyTeX not on PATH; PDF toolchain fragile** | Report fails to render | Quarto+TinyTeX inside `envs/report.yaml`; `render_report.py` falls back to self-contained HTML so the run never hard-fails |
| 12 | **No `mamba` on the cluster** | Default conda-frontend errors | `--conda-frontend conda` set in both profiles |
| 13 | **GPU not guaranteed on SLURM** | scVI/CellBender unavailable | GPU-**optional**: CPU defaults (Harmony, decontX, Scrublet, CPU XGBoost) give a complete result; GPU methods additive via `config.gpu.enabled` |
| 14 | **GEO RAW tar layout varies** (flat GSM-prefixed files vs folders) | Ingestion may mis-organize | `download_geo.py` regex-groups GSM files into triplet dirs; verify layout per dataset before the full run |
| 15 | **Large memory at integration/DE on ~100k cells** | OOM on full data | SLURM profile sets generous `mem_mb` + attempt-scaled retries on the ~2 TB nodes |

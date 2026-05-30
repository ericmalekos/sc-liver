# ADR-0001 — Core architecture decisions

Status: accepted · Date: 2026-05-29

## Context
Greenfield mini-pipeline to discover/prioritize cell-type-specific biomarkers in human liver
fibrosis (primary GSE136103), to be reproducible and modular and to answer 8 screening questions.

## Decisions
1. **Orchestrator = Snakemake; analysis = Python/Scanpy** (R only for pseudobulk DE, decontX,
   scDblFinder, NicheNet). Snakemake is installed on the target cluster; Python/Scanpy + scvi-tools
   is the strongest ecosystem for the AI/ML prioritization. nf-core conventions are adapted to
   Snakemake (per-rule envs, config schema, profiles, test fixture, CI).
2. **Validation dataset = GSE244832**, not SCP2154 (which *contains* GSE136103 → circular) nor
   GSE207310 (bulk, not single-cell). GSE244832 is independent, adds a fibrosis-severity gradient
   in an orthogonal MASH etiology and snRNA modality. (See `data_provenance.md`.)
3. **Donor-aware pseudobulk DE** is mandatory (donor = replication unit) to avoid pseudoreplication.
4. **Integrate over donor, never condition**, with a scIB fibrosis-signal guard.
5. **GPU-optional**: CPU defaults (Harmony, decontX, Scrublet, CPU XGBoost) give a complete result;
   scVI/scANVI/CellBender are additive via `config.gpu.enabled`.
6. **AnnData is the inter-stage currency; the R boundary is TSV/RDS**, never raw .h5ad into R.
7. **micromamba** (pinned 1.5.8, via a `mamba` shim) is the recommended env-creation frontend —
   far faster/more reliable than classic conda for the ~10 per-rule envs. conda remains the
   fallback default in profiles.
8. **Reporting = Quarto** with a self-contained HTML fallback so the pipeline never hard-fails on
   the report step.

## Consequences
- The pipeline runs end-to-end on CPU with no GPU and no external services (cached druggability
  snapshot, local gene sets) — important for reproducibility and offline/CI runs.
- The biggest correctness-sensitive boundary is `pseudobulk_aggregate.py` → `deseq2_pseudobulk.R`;
  it is exercised by the committed fixture and CI.
- Switching the validation dataset or DE engine is a config change, not a code change.

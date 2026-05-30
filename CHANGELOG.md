# Changelog

All notable changes to this project are documented here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/); versions follow SemVer.

## [Unreleased]
### Added
- Initial nf-core-style Snakemake + Scanpy pipeline scaffold.
- Reproducibility spine: per-rule conda envs, local + SLURM profiles, containers, CI.
- Curated liver cell-type / fibrosis reference (`docs/reference_liver_cell_types_fibrosis.md`).
- End-to-end modules: ingest → QC → integrate → annotate → donor-aware DE →
  pathway → cell-cell communication → ML biomarker prioritization → report.
- Written answers to the 8 screening questions (`docs/written_answers.md`).

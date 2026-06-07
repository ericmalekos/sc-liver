# =============================================================================
# Shared setup: config validation, samplesheet loading, dataset/sample/compartment
# helpers, wildcard constraints, and target builders used by `rule all`.
# =============================================================================
import os
import pandas as pd
from snakemake.utils import validate, min_version

min_version("7.0")

# ---- config validation -------------------------------------------------------
validate(config, "../../config/schemas/config.schema.yaml")

# ---- samplesheet -------------------------------------------------------------
SAMPLES = pd.read_csv(config["samples_tsv"], sep="\t", comment="#", dtype=str).dropna(
    how="all"
)
# tolerate an (almost) empty real samplesheet before make_samplesheet.py is run
if not SAMPLES.empty:
    SAMPLES["fibrosis_axis"] = SAMPLES["fibrosis_axis"].astype(int)
    for _, row in SAMPLES.iterrows():
        validate(row.to_dict(), "../../config/schemas/samples.schema.yaml")
    SAMPLES = SAMPLES.set_index("sample_id", drop=False)

# ---- dataset / sample / compartment views ------------------------------------
DATASETS = list(config["datasets"].keys())
PRIMARY_DS = next(d for d, v in config["datasets"].items() if v["role"] == "primary")
VALID_DS = [d for d, v in config["datasets"].items() if v["role"] == "validation"]


def samples_of(ds):
    """sample_ids belonging to dataset `ds` (empty list if none yet)."""
    if SAMPLES.empty:
        return []
    return SAMPLES.index[SAMPLES["dataset"] == ds].tolist()


def datasets_with_samples():
    return [d for d in DATASETS if samples_of(d)]


# analysis compartments come from the marker map (config-driven, excludes catch-all)
import yaml as _yaml

with open(config["annotation"]["markers_file"]) as _fh:
    _MARKERS = _yaml.safe_load(_fh)
COMPARTMENTS = [c for c in _MARKERS["compartments"].keys() if c != "myeloid_other"]


def contrast_name(ds):
    return config["de"]["contrasts"][ds]["name"]


# ---- ingest input: gate on download sentinel for real runs, source dir for fixture
def ingest_input(wildcards):
    dp = str(SAMPLES.loc[wildcards.sample, "data_path"])
    # already staged (committed fixture, or pre-split datasets like GSE244832) -> use directly
    if os.path.exists(os.path.join(dp, "matrix.mtx.gz")):
        return dp
    # otherwise gate ingestion on the dataset download sentinel
    if dp.startswith("results/00_download/"):
        return f"results/00_download/{wildcards.ds}/DONE"
    return dp


def sample_data_path(wildcards):
    return str(SAMPLES.loc[wildcards.sample, "data_path"])


# ---- wildcard hygiene --------------------------------------------------------
wildcard_constraints:
    ds="|".join(DATASETS),
    sample="|".join(SAMPLES.index.tolist()) if not SAMPLES.empty else "NOSAMPLES",
    compartment="|".join(COMPARTMENTS),


# ---- path helpers ------------------------------------------------------------
# `script:` and `conda:` paths resolve relative to the rule file (workflow/rules/),
# so scripts (workflow/scripts/) and envs (envs/) need the prefixes below.
def script_py(name):
    return f"../scripts/py/{name}"


def script_r(name):
    return f"../scripts/R/{name}"


def env(name):
    return f"../../envs/{name}"


# ---- QC method dispatch (pick script + conda env by config, at parse time) ----
# decontX/SoupX are R and work on FILTERED matrices (what GEO ships); CellBender is
# Python and needs raw droplets; `none` is a Python passthrough (fixture/CI).
def ambient_script():
    m = config["qc"]["ambient"]["method"]
    if m == "decontx":
        return script_py("qc_ambient_decontx.py")  # light Python<->R hand-off (no basilisk)
    if m == "soupx":
        return script_r("qc_ambient.R")  # needs raw droplets; not for filtered GEO data
    return script_py("qc_ambient.py")


def ambient_env():
    m = config["qc"]["ambient"]["method"]
    if m == "decontx":
        return env("qc_decontx.yaml")
    if m == "soupx":
        return env("qc_r.yaml")
    if m == "cellbender":
        return env("cellbender.yaml")
    return env("qc_py.yaml")


def doublet_script():
    return (
        script_r("scdblfinder.R")
        if config["qc"]["doublets"]["method"] == "scdblfinder"
        else script_py("qc_doublets.py")
    )


def doublet_env():
    return (
        env("qc_r.yaml")
        if config["qc"]["doublets"]["method"] == "scdblfinder"
        else env("qc_py.yaml")
    )


# ---- final targets for `rule all` --------------------------------------------
def final_targets():
    t = []
    ds_ok = datasets_with_samples()
    # QC + integration + annotation per dataset
    for ds in ds_ok:
        t.append(f"results/02_qc/{ds}/qc_summary.tsv")
        t.append(f"results/03_integrate/{ds}/scib_metrics.tsv")
        t.append(f"results/04_annotate/{ds}/compartment_validation.tsv")
        if config["niche"]["enabled"]:
            t.append(f"results/05_de/{ds}/niche_summary.tsv")
        for comp in COMPARTMENTS:
            t.append(f"results/05_de/{ds}/{comp}/de_results.tsv")
            t.append(f"results/06_pathway/{ds}/{comp}/gsea.tsv")
        t.append(f"results/06_pathway/{ds}/activity.tsv")
        t.append(f"results/07_ccc/{ds}/liana_consensus.tsv")
        t.append(f"results/07_ccc/{ds}/differential_ccc.tsv")
    # convergence (only if the primary dataset has samples)
    if PRIMARY_DS in ds_ok:
        t.append("results/08_score/candidate_scores.tsv")
        t.append("results/08_score/known_positive_recall.tsv")
        if config["crossdataset"]["enabled"]:
            t.append("results/09_crossdataset/repro_scores.tsv")
        t.append("results/10_report/report.html")
        if "pdf" in config["report"]["formats"]:
            t.append("results/10_report/report.pdf")
            t.append("results/10_report/exec_summary.pdf")
    t.append("results/provenance/run_metadata.json")
    return t

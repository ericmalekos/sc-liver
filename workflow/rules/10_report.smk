# Stage 10 — Quarto report + executive summary, and run provenance.
_report_out = [
    "results/10_report/report.html",
    "results/10_report/exec_summary.html",
]
if "pdf" in config["report"]["formats"]:
    _report_out += [
        "results/10_report/report.pdf",
        "results/10_report/exec_summary.pdf",
    ]


def report_inputs(wc):
    ds_ok = datasets_with_samples()
    inp = {
        "scores": "results/08_score/candidate_scores.tsv",
        "ml_metrics": "results/08_score/ml_metrics.json",
        "qmd": "workflow/report/report.qmd",
        "exec_qmd": "workflow/report/exec_summary.qmd",
        "qc": expand("results/02_qc/{ds}/qc_summary.tsv", ds=ds_ok),
        "scib": expand("results/03_integrate/{ds}/scib_metrics.tsv", ds=ds_ok),
        "validation": expand(
            "results/04_annotate/{ds}/compartment_validation.tsv", ds=ds_ok
        ),
        "activity": expand("results/06_pathway/{ds}/activity.tsv", ds=ds_ok),
        "ccc": expand("results/07_ccc/{ds}/differential_ccc.tsv", ds=ds_ok),
    }
    if config["crossdataset"]["enabled"]:
        inp["repro"] = "results/09_crossdataset/repro_scores.tsv"
    return inp


rule render_report:
    input:
        unpack(report_inputs),
    output:
        _report_out,
    params:
        formats=config["report"]["formats"],
        title=config["report"]["title"],
        primary=PRIMARY_DS,
        datasets=lambda wc: datasets_with_samples(),
    log:
        "logs/10_report/render.log",
    conda:
        env("report.yaml")
    script:
        script_py("render_report.py")


rule write_provenance:
    output:
        meta="results/provenance/run_metadata.json",
        cfg="results/provenance/config_resolved.yaml",
    params:
        config=config,
    log:
        "logs/provenance.log",
    conda:
        env("ingest.yaml")
    script:
        script_py("write_provenance.py")

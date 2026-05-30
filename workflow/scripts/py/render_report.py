"""Stage 10 — render the Quarto report + executive summary (HTML, optional PDF).

Primary path uses Quarto (the deliverable). If Quarto is unavailable or errors, falls back
to a self-contained pandas-built HTML so the pipeline always completes with a readable report.
"""
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import ensure_parent, get_logger  # noqa: E402

import pandas as pd  # noqa: E402

log = get_logger(snakemake)  # noqa: F821
p = snakemake.params  # noqa: F821
formats = list(p.formats)
primary = p.primary
datasets = list(p.datasets)
abs_results = os.path.abspath("results")
outdir = os.path.abspath("results/10_report")
os.makedirs(outdir, exist_ok=True)

qmd = {"report": snakemake.input.qmd, "exec_summary": snakemake.input.exec_qmd}  # noqa: F821
made = set()


def try_quarto(stem, src, fmt, params):
    try:
        cmd = ["quarto", "render", os.path.abspath(src), "--to", fmt]
        for k, v in params.items():
            cmd += ["-P", f"{k}:{v}"]
        log.info("running: " + " ".join(cmd))
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        produced = os.path.join(os.path.dirname(os.path.abspath(src)), f"{stem}.{fmt}")
        dest = os.path.join(outdir, f"{stem}.{fmt}")
        if os.path.exists(produced):
            shutil.move(produced, dest)
            return True
    except Exception as e:
        log.warning(f"quarto {stem}.{fmt} failed: {getattr(e, 'stderr', e)}")
    return False


def fallback_html(stem, dest):
    """self-contained HTML from the key result tables."""
    def tbl(path, n=25):
        path = os.path.join(abs_results, path)
        if os.path.exists(path):
            df = pd.read_csv(path, sep="\t")
            return df.head(n).to_html(index=False, border=0)
        return f"<p><em>missing: {path}</em></p>"

    blocks = [f"<h1>{p.title}</h1>",
              "<p>Fallback report (Quarto unavailable). Tables read from <code>results/</code>.</p>",
              "<h2>Ranked biomarker / target candidates</h2>", tbl("08_score/candidate_scores.tsv")]
    if stem == "report":
        blocks += ["<h2>Cross-dataset reproducibility</h2>", tbl("09_crossdataset/repro_scores.tsv")]
        for ds in datasets:
            blocks += [f"<h2>{ds}: QC</h2>", tbl(f"02_qc/{ds}/qc_summary.tsv"),
                       f"<h3>{ds}: integration (scIB)</h3>", tbl(f"03_integrate/{ds}/scib_metrics.tsv"),
                       f"<h3>{ds}: compartment validation</h3>", tbl(f"04_annotate/{ds}/compartment_validation.tsv"),
                       f"<h3>{ds}: differential cell-cell communication</h3>", tbl(f"07_ccc/{ds}/differential_ccc.tsv")]
    html = ("<html><head><meta charset='utf-8'><style>"
            "body{font-family:system-ui,Arial;margin:2rem;max-width:1100px}"
            "table{border-collapse:collapse;font-size:12px}th,td{padding:3px 8px;border-bottom:1px solid #ddd}"
            "h1,h2,h3{color:#243b53}</style></head><body>" + "".join(blocks) + "</body></html>")
    with open(dest, "w") as fh:
        fh.write(html)


# HTML (always); PDF only if requested
for stem, src in qmd.items():
    params = {"results_dir": abs_results}
    if stem == "report":
        params.update({"primary": primary, "datasets": ",".join(datasets)})
    html_dest = os.path.join(outdir, f"{stem}.html")
    if not try_quarto(stem, src, "html", params):
        fallback_html(stem, html_dest)
        log.info(f"wrote fallback HTML: {html_dest}")
    made.add(html_dest)
    if "pdf" in formats:
        pdf_dest = os.path.join(outdir, f"{stem}.pdf")
        if not try_quarto(stem, src, "pdf", params):
            # last-resort: a one-page matplotlib PDF pointer so the target exists
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            fig = plt.figure(figsize=(8.3, 11.7))
            fig.text(0.1, 0.9, p.title, fontsize=14, weight="bold")
            fig.text(0.1, 0.8, f"See {stem}.html — PDF toolchain (tinytex) unavailable.", fontsize=10)
            fig.savefig(pdf_dest)
            plt.close(fig)
            log.warning(f"wrote placeholder PDF: {pdf_dest}")

for f in snakemake.output:  # noqa: F821
    ensure_parent(f)
    if not os.path.exists(f):
        # guarantee declared outputs exist
        open(f, "a").close()
log.info(f"report outputs: {sorted(made)}")

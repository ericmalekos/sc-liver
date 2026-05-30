# Stage 00 — fetch + verify raw data. Only runs for real datasets (the fixture's
# data_path points at committed files, so ingest_input does not request DONE).
rule download_dataset:
    output:
        sentinel="results/00_download/{ds}/DONE",
    params:
        url=lambda wc: config["datasets"][wc.ds]["supp_url"],
        accession=lambda wc: config["datasets"][wc.ds]["geo_accession"],
        md5=lambda wc: config["datasets"][wc.ds].get("md5"),
        outdir=lambda wc: f"results/00_download/{wc.ds}/raw",
    log:
        "logs/00_download/{ds}.log",
    threads: 2
    conda:
        env("ingest.yaml")
    script:
        script_py("download_geo.py")

# Stage 03 — integration over donor/sample (NOT condition) + scIB guard (Q3).
rule integrate:
    input:
        h5ads=lambda wc: expand(
            "results/02_qc/{ds}/{sample}.qc.h5ad",
            ds=wc.ds,
            sample=samples_of(wc.ds),
        ),
        sheet="results/01_ingest/{ds}/harmonized_samplesheet.tsv",
    output:
        h5ad="results/03_integrate/{ds}/integrated.h5ad",
        umap="results/03_integrate/{ds}/figures/umap_overview.png",
    params:
        method=config["integration"]["method"],
        batch_key=config["integration"]["batch_key"],
        n_hvg=config["integration"]["n_hvg"],
        n_latent=config["integration"]["n_latent"],
        n_neighbors=config["integration"]["n_neighbors"],
        gpu=config["gpu"]["enabled"],
        seed=config["seeds"]["global"],
    log:
        "logs/03_integrate/{ds}.log",
    threads: 8
    conda:
        (
            env("integrate_gpu.yaml")
            if config["gpu"]["enabled"]
            and config["integration"]["method"] in ("scvi", "scanvi")
            else env("integrate_cpu.yaml")
        )
    script:
        script_py("integrate.py")


rule scib_benchmark:
    input:
        h5ad="results/03_integrate/{ds}/integrated.h5ad",
    output:
        tsv="results/03_integrate/{ds}/scib_metrics.tsv",
        fig="results/03_integrate/{ds}/figures/scib_summary.png",
    params:
        batch_key=config["integration"]["batch_key"],
        guard=config["integration"]["fibrosis_signal_guard"],
    log:
        "logs/03_integrate/{ds}.scib.log",
    conda:
        env("integrate_cpu.yaml")
    script:
        script_py("scib_benchmark.py")

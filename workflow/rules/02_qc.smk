# Stage 02 — QC (Q2): ambient RNA -> doublets -> adaptive MAD filter.
# Method/language dispatch (decontX=R, scrublet=py, ...) via common.smk helpers.
rule qc_ambient:
    input:
        h5ad="results/01_ingest/{ds}/{sample}.h5ad",
    output:
        h5ad="results/02_qc/{ds}/{sample}.ambient.h5ad",
    params:
        method=config["qc"]["ambient"]["method"],
        modality=lambda wc: SAMPLES.loc[wc.sample, "modality"],
        gpu=config["gpu"]["enabled"],
    log:
        "logs/02_qc/{ds}.{sample}.ambient.log",
    conda:
        ambient_env()
    script:
        ambient_script()


rule qc_doublets:
    input:
        h5ad="results/02_qc/{ds}/{sample}.ambient.h5ad",
    output:
        h5ad="results/02_qc/{ds}/{sample}.dbl.h5ad",
    params:
        method=config["qc"]["doublets"]["method"],
        expected_rate=config["qc"]["doublets"]["expected_rate"],
        seed=config["seeds"]["global"],
    log:
        "logs/02_qc/{ds}.{sample}.doublets.log",
    conda:
        doublet_env()
    script:
        doublet_script()


rule qc_filter:
    input:
        h5ad="results/02_qc/{ds}/{sample}.dbl.h5ad",
    output:
        h5ad="results/02_qc/{ds}/{sample}.qc.h5ad",
        metrics="results/02_qc/{ds}/qc_metrics/{sample}.json",
        fig="results/02_qc/{ds}/figures/{sample}_qc.png",
    params:
        n_mads=config["qc"]["n_mads"],
        mito_pct_max=config["qc"]["mito_pct_max"],
        min_genes=config["qc"]["min_genes"],
        min_counts=config["qc"]["min_counts"],
        min_cells=config["qc"]["min_cells"],
        keep_stressed=config["qc"]["keep_stressed"],
        modality=lambda wc: SAMPLES.loc[wc.sample, "modality"],
    log:
        "logs/02_qc/{ds}.{sample}.filter.log",
    conda:
        env("qc_py.yaml")
    script:
        script_py("qc_filter.py")


rule qc_summary:
    input:
        metrics=lambda wc: expand(
            "results/02_qc/{ds}/qc_metrics/{sample}.json",
            ds=wc.ds,
            sample=samples_of(wc.ds),
        ),
    output:
        tsv="results/02_qc/{ds}/qc_summary.tsv",
    log:
        "logs/02_qc/{ds}.summary.log",
    conda:
        env("qc_py.yaml")
    script:
        script_py("qc_summary.py")

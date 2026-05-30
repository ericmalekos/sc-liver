# Stage 01 — per-sample 10x -> AnnData, and per-dataset metadata harmonization (Q1).
rule ingest_sample:
    input:
        ingest_input,
    output:
        h5ad="results/01_ingest/{ds}/{sample}.h5ad",
    params:
        data_path=sample_data_path,
        meta=lambda wc: SAMPLES.loc[wc.sample].to_dict(),
    log:
        "logs/01_ingest/{ds}.{sample}.log",
    conda:
        env("ingest.yaml")
    script:
        script_py("ingest_10x.py")


rule harmonize_metadata:
    input:
        h5ads=lambda wc: expand(
            "results/01_ingest/{ds}/{sample}.h5ad",
            ds=wc.ds,
            sample=samples_of(wc.ds),
        ),
    output:
        sheet="results/01_ingest/{ds}/harmonized_samplesheet.tsv",
    params:
        ds=lambda wc: wc.ds,
        samples=lambda wc: samples_of(wc.ds),
        fibrosis_mapping=config["metadata"]["fibrosis_mapping"],
        f2plus=config["metadata"]["f2plus_threshold"],
    log:
        "logs/01_ingest/{ds}.harmonize.log",
    conda:
        env("ingest.yaml")
    script:
        script_py("harmonize_metadata.py")

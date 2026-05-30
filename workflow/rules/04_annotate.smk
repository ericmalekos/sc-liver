# Stage 04 — clustering + annotation + compartment validation (Q4).
rule cluster_annotate:
    input:
        h5ad="results/03_integrate/{ds}/integrated.h5ad",
        markers=config["annotation"]["markers_file"],
    output:
        h5ad="results/04_annotate/{ds}/clustered.h5ad",
        dotplot="results/04_annotate/{ds}/figures/marker_dotplot.png",
    params:
        resolution=config["annotation"]["leiden_resolution"],
        seed=config["seeds"]["global"],
    log:
        "logs/04_annotate/{ds}.cluster.log",
    threads: 8
    conda:
        env("annotate.yaml")
    script:
        script_py("cluster_annotate.py")


rule celltypist_decoupler_annotate:
    input:
        h5ad="results/04_annotate/{ds}/clustered.h5ad",
        markers=config["annotation"]["markers_file"],
    output:
        h5ad="results/04_annotate/{ds}/annotated.h5ad",
        validation="results/04_annotate/{ds}/compartment_validation.tsv",
        umap="results/04_annotate/{ds}/figures/umap_celltypes.png",
    params:
        celltypist_model=config["annotation"]["celltypist_model"],
        ora_min_overlap=config["annotation"]["decoupler_ora_min_overlap"],
    log:
        "logs/04_annotate/{ds}.annotate.log",
    conda:
        env("annotate.yaml")
    script:
        script_py("celltypist_annotate.py")

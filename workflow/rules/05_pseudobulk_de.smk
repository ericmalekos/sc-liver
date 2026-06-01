# Stage 05 — DONOR-AWARE pseudobulk differential expression (Q5; highest-risk).
# Python aggregates to donor x condition counts; R runs DESeq2/edgeR with donor as
# the replication unit (avoids cell-level pseudoreplication, Squair 2021).
rule pseudobulk_aggregate:
    input:
        h5ad="results/04_annotate/{ds}/annotated.h5ad",
    output:
        counts="results/05_de/{ds}/{compartment}/pseudobulk_counts.tsv",
        coldata="results/05_de/{ds}/{compartment}/coldata.tsv",
    params:
        compartment=lambda wc: wc.compartment,
        min_cells_per_donor=config["de"]["min_cells_per_donor"],
        aggregate=config["de"]["aggregate"],
    log:
        "logs/05_de/{ds}.{compartment}.aggregate.log",
    conda:
        env("annotate.yaml")
    script:
        script_py("pseudobulk_aggregate.py")


rule deseq2_de:
    input:
        counts="results/05_de/{ds}/{compartment}/pseudobulk_counts.tsv",
        coldata="results/05_de/{ds}/{compartment}/coldata.tsv",
    output:
        de="results/05_de/{ds}/{compartment}/de_results.tsv",
    params:
        engine=config["de"]["engine"],
        design=config["de"]["design"],
        covariates=config["de"]["covariates"],
        contrast=lambda wc: config["de"]["contrasts"][wc.ds],
        min_donors=config["de"]["min_donors_per_group"],
        lfc_threshold=config["de"]["lfc_threshold"],
    log:
        "logs/05_de/{ds}.{compartment}.deseq2.log",
    conda:
        env("de_r.yaml")
    script:
        script_r("deseq2_pseudobulk.R")


# Stage 05c — niche (subpopulation) markers: subcluster within each compartment, flag
# disease-enriched subclusters (donor-aware), and take their one-vs-rest up-markers as the
# subset biomarkers that compartment-level pseudobulk buries (PLVAP, TREM2, ...).
rule niche_markers:
    input:
        h5ad="results/04_annotate/{ds}/annotated.h5ad",
    output:
        summary="results/05_de/{ds}/niche_summary.tsv",
        markers="results/05_de/{ds}/niche_markers.tsv",
    params:
        seed=config["seeds"]["leiden"],
        contrast=lambda wc: config["de"]["contrasts"][wc.ds],
        compartments=COMPARTMENTS,
        resolution=config["niche"]["resolution"],
        min_cells_subcluster=config["niche"]["min_cells_subcluster"],
        enrichment_min_log2fc=config["niche"]["enrichment_min_log2fc"],
        marker_padj=config["niche"]["marker_padj"],
        f2plus_threshold=config["metadata"]["f2plus_threshold"],
    log:
        "logs/05_de/{ds}.niche.log",
    conda:
        env("annotate.yaml")
    script:
        script_py("niche_markers.py")

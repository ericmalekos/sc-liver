# Stage 07 — cell-cell communication (Q7) with overinterpretation guards.
rule ccc_liana:
    input:
        h5ad="results/04_annotate/{ds}/annotated.h5ad",
    output:
        liana="results/07_ccc/{ds}/liana_consensus.tsv",
    params:
        resource=config["ccc"]["resource"],
        expr_prop=config["ccc"]["expr_prop"],
        min_cells=config["ccc"]["min_cells"],
        seed=config["seeds"]["global"],
    log:
        "logs/07_ccc/{ds}.liana.log",
    threads: 8
    conda:
        env("ccc_py.yaml")
    script:
        script_py("ccc_liana.py")


rule ccc_differential:
    input:
        liana="results/07_ccc/{ds}/liana_consensus.tsv",
        de=lambda wc: expand(
            "results/05_de/{ds}/{compartment}/de_results.tsv",
            ds=wc.ds,
            compartment=COMPARTMENTS,
        ),
    output:
        diff="results/07_ccc/{ds}/differential_ccc.tsv",
        nichenet="results/07_ccc/{ds}/nichenet_links.tsv",
    params:
        specificity_cutoff=config["ccc"]["specificity_cutoff"],
        nichenet=config["ccc"]["nichenet_corroborate"],
        min_targets=config["ccc"]["min_downstream_targets"],
    log:
        "logs/07_ccc/{ds}.differential.log",
    conda:
        env("ccc_py.yaml")
    script:
        script_py("ccc_differential.py")

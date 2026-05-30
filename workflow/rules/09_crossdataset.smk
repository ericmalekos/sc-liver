# Stage 09 — cross-dataset reproducibility (Q1 validation arm). Compares DE
# direction/rank between the primary and validation datasets (robust to the
# scRNA vs snRNA modality gap).
def crossdataset_input(wc):
    valid_present = [d for d in VALID_DS if samples_of(d)]
    return {
        "primary": expand(
            "results/05_de/{ds}/{compartment}/de_results.tsv",
            ds=PRIMARY_DS,
            compartment=COMPARTMENTS,
        ),
        "validation": expand(
            "results/05_de/{ds}/{compartment}/de_results.tsv",
            ds=valid_present,
            compartment=COMPARTMENTS,
        ),
    }


rule crossdataset_repro:
    input:
        unpack(crossdataset_input),
    output:
        repro="results/09_crossdataset/repro_scores.tsv",
        fig="results/09_crossdataset/figures/concordance.png",
    params:
        primary=PRIMARY_DS,
        validation=[d for d in VALID_DS if samples_of(d)],
        compartments=COMPARTMENTS,
        rank_metric=config["crossdataset"]["rank_metric"],
        match_on=config["crossdataset"]["match_on"],
    log:
        "logs/09_crossdataset/repro.log",
    conda:
        env("repro.yaml")
    script:
        script_py("crossdataset_repro.py")

# Stage 06 — pathway / mechanism (Q7-mechanism).
# GSEA per compartment on the DE ranking; PROGENy + CollecTRI + per-cell scores per dataset.
rule pathway_gsea:
    input:
        de="results/05_de/{ds}/{compartment}/de_results.tsv",
    output:
        gsea="results/06_pathway/{ds}/{compartment}/gsea.tsv",
    params:
        genesets_dir=config["pathway"]["genesets_dir"],
        gsea_sets=config["pathway"]["gsea_sets"],
        seed=config["seeds"]["global"],
    log:
        "logs/06_pathway/{ds}.{compartment}.gsea.log",
    conda:
        env("pathway.yaml")
    script:
        script_py("gsea.py")


rule pathway_activity:
    input:
        h5ad="results/04_annotate/{ds}/annotated.h5ad",
    output:
        activity="results/06_pathway/{ds}/activity.tsv",
        scores="results/06_pathway/{ds}/geneset_scores.tsv",
        fig="results/06_pathway/{ds}/figures/pathway_activity.png",
    params:
        progeny_top=config["pathway"]["progeny_top"],
        collectri=config["pathway"]["collectri"],
        score_method=config["pathway"]["score_method"],
        genesets_dir=config["pathway"]["genesets_dir"],
    log:
        "logs/06_pathway/{ds}.activity.log",
    conda:
        env("pathway.yaml")
    script:
        script_py("pathway_decoupler.py")

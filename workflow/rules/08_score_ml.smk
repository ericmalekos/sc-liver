# Stage 08 — AI/ML biomarker prioritization (Q6). Operates on the PRIMARY dataset,
# folding in cross-dataset reproducibility, druggability, accessibility, and SHAP.
def build_features_input(wc):
    inp = {
        "de": expand(
            "results/05_de/{ds}/{compartment}/de_results.tsv",
            ds=PRIMARY_DS,
            compartment=COMPARTMENTS,
        ),
        "pathway": f"results/06_pathway/{PRIMARY_DS}/activity.tsv",
        "annotated": f"results/04_annotate/{PRIMARY_DS}/annotated.h5ad",
    }
    if config["crossdataset"]["enabled"]:
        inp["repro"] = "results/09_crossdataset/repro_scores.tsv"
    return inp


rule build_features:
    input:
        unpack(build_features_input),
    output:
        features="results/08_score/feature_matrix.tsv",
    params:
        compartments=COMPARTMENTS,
        specificity_index=config["score"]["specificity_index"],
        druggability_cache=config["score"]["druggability_cache"],
        druggability_sources=config["score"]["druggability_sources"],
        accessibility=config["score"]["accessibility"],
        markers_file=config["annotation"]["markers_file"],
    log:
        "logs/08_score/build_features.log",
    conda:
        env("ml.yaml")
    script:
        script_py("build_features.py")


rule ml_classifier_shap:
    input:
        features="results/08_score/feature_matrix.tsv",
        counts=expand(
            "results/05_de/{ds}/{compartment}/pseudobulk_counts.tsv",
            ds=PRIMARY_DS,
            compartment=COMPARTMENTS,
        ),
        coldata=expand(
            "results/05_de/{ds}/{compartment}/coldata.tsv",
            ds=PRIMARY_DS,
            compartment=COMPARTMENTS,
        ),
    output:
        shap="results/08_score/shap_values.tsv",
        model="results/08_score/ml_model.pkl",
        metrics="results/08_score/ml_metrics.json",
    params:
        model=config["ml"]["model"],
        cv_folds=config["ml"]["cv_folds"],
        group_key=config["ml"]["group_key"],
        contrast=config["de"]["contrasts"][PRIMARY_DS],
        compartments=COMPARTMENTS,
        seed=config["seeds"]["ml"],
    log:
        "logs/08_score/ml_shap.log",
    threads: 8
    conda:
        env("ml.yaml")
    script:
        script_py("ml_classifier_shap.py")


rule compute_score:
    input:
        features="results/08_score/feature_matrix.tsv",
        shap="results/08_score/shap_values.tsv",
    output:
        scores="results/08_score/candidate_scores.tsv",
        fig="results/08_score/figures/top_candidates.png",
    params:
        weights=config["score"]["weights"],
        top_n_per_compartment=config["score"]["top_n_per_compartment"],
        top_n_overall=config["score"]["top_n_overall"],
        required_compartments=config["score"]["required_compartments"],
        tie_breaker=config["score"]["tie_breaker"],
        require_de_significance=config["score"]["require_de_significance"],
        require_home_compartment=config["score"]["require_home_compartment"],
        specificity_ratio_floor=config["score"]["specificity_ratio_floor"],
        padj_threshold=config["de"]["padj_threshold"],
        lfc_threshold=config["de"]["lfc_threshold"],
    log:
        "logs/08_score/compute_score.log",
    conda:
        env("ml.yaml")
    script:
        script_py("compute_score.py")


# Sanity check (Q6): does the ranked list recover known literature positives, and if a
# known positive is missed, WHY (non-significant / filtered / absent in primary)? Each is
# cross-referenced against the validation DE to flag primary-anchored false negatives.
def known_positive_recall_input(wc):
    valid_present = [d for d in VALID_DS if samples_of(d)]
    return {
        "candidates": "results/08_score/candidate_scores.tsv",
        "primary_de": expand(
            "results/05_de/{ds}/{compartment}/de_results.tsv",
            ds=PRIMARY_DS,
            compartment=COMPARTMENTS,
        ),
        "primary_counts": expand(
            "results/05_de/{ds}/{compartment}/pseudobulk_counts.tsv",
            ds=PRIMARY_DS,
            compartment=COMPARTMENTS,
        ),
        "valid_de": expand(
            "results/05_de/{ds}/{compartment}/de_results.tsv",
            ds=valid_present,
            compartment=COMPARTMENTS,
        ),
    }


rule known_positive_recall:
    input:
        unpack(known_positive_recall_input),
    output:
        recall="results/08_score/known_positive_recall.tsv",
        summary="results/08_score/known_positive_recall_summary.tsv",
        fig="results/08_score/figures/known_positive_recall.png",
    params:
        panel=config["score"]["known_positives"],
        padj_threshold=config["de"]["padj_threshold"],
    log:
        "logs/08_score/known_positive_recall.log",
    conda:
        env("ml.yaml")
    script:
        script_py("known_positive_recall.py")

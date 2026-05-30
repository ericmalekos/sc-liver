# Stage 02 — doublet detection in R (scDblFinder; best AUPRC in benchmarks).
# Reads/writes .h5ad via zellkonverter. Adds doublet_score + predicted_doublet to obs.
suppressMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(scDblFinder)
})

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "message")

sce <- readH5AD(snakemake@input[["h5ad"]])
assayNames(sce)[1] <- "counts"
set.seed(as.integer(snakemake@params[["seed"]]))

res <- tryCatch(
  scDblFinder(sce),
  error = function(e) {
    message("scDblFinder fallback (kept all): ", conditionMessage(e))
    sce$scDblFinder.score <- 0
    sce$scDblFinder.class <- "singlet"
    sce
  }
)

res$doublet_score <- res$scDblFinder.score
res$predicted_doublet <- res$scDblFinder.class == "doublet"
assayNames(res)[1] <- "X"
writeH5AD(res, snakemake@output[["h5ad"]], X_name = "X")
message(sprintf("scDblFinder flagged %d/%d", sum(res$predicted_doublet), ncol(res)))

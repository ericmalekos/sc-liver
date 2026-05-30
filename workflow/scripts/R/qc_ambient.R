# Stage 02 — ambient RNA removal in R (decontX; SoupX falls back to decontX since the
# GEO matrices are already filtered and SoupX needs raw droplets). Reads/writes .h5ad
# via zellkonverter. Robust to degenerate inputs (passthrough on error).
suppressMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(decontX)
})

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "message")

infile <- snakemake@input[["h5ad"]]
outfile <- snakemake@output[["h5ad"]]
method <- snakemake@params[["method"]]

if (identical(method, "soupx")) {
  message("SoupX needs raw droplets; GEO ships filtered matrices -> using decontX instead.")
}

sce <- readH5AD(infile)
counts_mat <- assay(sce, 1)

ok <- tryCatch({
  res <- decontX(as(counts_mat, "CsparseMatrix"))
  dec <- round(res$decontXcounts)
  assay(sce, 1) <- dec
  sce$ambient_fraction <- res$contamination
  TRUE
}, error = function(e) {
  message("decontX failed, passing counts through: ", conditionMessage(e))
  sce$ambient_fraction <- 0
  FALSE
})

assayNames(sce)[1] <- "X"
writeH5AD(sce, outfile, X_name = "X")
message(sprintf("ambient(decontX) done=%s, cells=%d", ok, ncol(sce)))

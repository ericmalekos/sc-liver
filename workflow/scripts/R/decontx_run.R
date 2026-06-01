# decontX on a plain MatrixMarket counts file (genes x cells) — no zellkonverter/basilisk.
# Usage: Rscript decontx_run.R <counts.mtx> <out_prefix> [z_clusters.txt]
# Writes <out_prefix>_decontaminated.mtx (rounded corrected counts, genes x cells) and
# <out_prefix>_contamination.txt (per-cell contamination fraction).
# If a cluster-label file is given, it is passed as decontX `z` priors so cross-population
# ambient is removed without stripping a population's genuine marker genes.
suppressMessages({
  library(Matrix)
  library(decontX)
})
args <- commandArgs(trailingOnly = TRUE)
m <- readMM(args[[1]])                 # genes x cells
z <- NULL
if (length(args) >= 3 && nzchar(args[[3]]) && file.exists(args[[3]])) {
  z <- scan(args[[3]], what = integer(), quiet = TRUE)
  if (length(z) != ncol(m) || length(unique(z)) < 2) z <- NULL   # fall back to self-clustering
}
res <- if (is.null(z)) decontX(as(m, "CsparseMatrix")) else decontX(as(m, "CsparseMatrix"), z = z)
writeMM(as(round(res$decontXcounts), "CsparseMatrix"), paste0(args[[2]], "_decontaminated.mtx"))
write.table(res$contamination, paste0(args[[2]], "_contamination.txt"),
            row.names = FALSE, col.names = FALSE)
cat(sprintf("decontX done: %d genes x %d cells, mean contamination %.3f\n",
            nrow(m), ncol(m), mean(res$contamination)))

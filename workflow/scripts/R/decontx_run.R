# decontX on a plain MatrixMarket counts file (genes x cells) — no zellkonverter/basilisk.
# Usage: Rscript decontx_run.R <counts.mtx> <out_prefix>
# Writes <out_prefix>_decontaminated.mtx (rounded corrected counts, genes x cells) and
# <out_prefix>_contamination.txt (per-cell contamination fraction).
suppressMessages({
  library(Matrix)
  library(decontX)
})
args <- commandArgs(trailingOnly = TRUE)
m <- readMM(args[[1]])                 # genes x cells
res <- decontX(as(m, "CsparseMatrix"))
writeMM(as(round(res$decontXcounts), "CsparseMatrix"), paste0(args[[2]], "_decontaminated.mtx"))
write.table(res$contamination, paste0(args[[2]], "_contamination.txt"),
            row.names = FALSE, col.names = FALSE)
cat(sprintf("decontX done: %d genes x %d cells, mean contamination %.3f\n",
            nrow(m), ncol(m), mean(res$contamination)))

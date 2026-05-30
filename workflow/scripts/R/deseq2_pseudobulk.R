# Stage 05b — donor-aware differential expression on pseudobulk counts (Q5).
# Replication unit = donor. Uses DESeq2 (default) / edgeR; covariates added only if they
# vary and are not confounded with the contrast group. Writes an empty (header-only)
# result with a status note when a group has too few donors, so the DAG always completes.
suppressMessages({
  library(jsonlite)
})

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "message")

counts_f <- snakemake@input[["counts"]]
coldata_f <- snakemake@input[["coldata"]]
out_f <- snakemake@output[["de"]]
engine <- snakemake@params[["engine"]]
contrast <- snakemake@params[["contrast"]]
covariates <- snakemake@params[["covariates"]]
min_donors <- as.integer(snakemake@params[["min_donors"]])

grp <- contrast[["group"]]; ref <- contrast[["ref"]]; tst <- contrast[["test"]]

empty_out <- function(note) {
  message("WRITE EMPTY: ", note)
  df <- data.frame(gene = character(), baseMean = numeric(), log2FoldChange = numeric(),
                   lfcSE = numeric(), stat = numeric(), pvalue = numeric(),
                   padj = numeric(), contrast = character(), status = character())
  write.table(df, out_f, sep = "\t", quote = FALSE, row.names = FALSE)
  quit(save = "no", status = 0)
}

counts <- tryCatch(as.matrix(read.delim(counts_f, row.names = 1, check.names = FALSE)),
                   error = function(e) matrix(nrow = 0, ncol = 0))
coldata <- tryCatch(read.delim(coldata_f, row.names = 1, check.names = FALSE),
                    error = function(e) data.frame())

if (nrow(counts) == 0 || ncol(counts) == 0 || nrow(coldata) == 0) empty_out("no pseudobulk data")
if (!grp %in% colnames(coldata)) empty_out(paste("group column missing:", grp))

coldata[[grp]] <- factor(coldata[[grp]])
if (!all(c(ref, tst) %in% levels(coldata[[grp]]))) empty_out("contrast levels absent")
tab <- table(coldata[[grp]])
if (any(tab[c(ref, tst)] < min_donors)) empty_out(paste0("too few donors per group: ",
                                                         paste(names(tab), tab, collapse = ", ")))

# design: varying, non-confounded covariates + group (group last)
use_cov <- c()
for (cv in covariates) {
  if (cv %in% colnames(coldata) && length(unique(coldata[[cv]])) > 1) {
    # skip if perfectly confounded with the group
    if (!any(table(coldata[[cv]], coldata[[grp]]) == 0 &
             rowSums(table(coldata[[cv]], coldata[[grp]])) > 0)) {
      coldata[[cv]] <- factor(coldata[[cv]]); use_cov <- c(use_cov, cv)
    }
  }
}
design <- as.formula(paste0("~ ", paste(c(use_cov, grp), collapse = " + ")))
coldata[[grp]] <- relevel(coldata[[grp]], ref = ref)
message("design: ", deparse(design), " | n=", ncol(counts),
        " | groups: ", paste(names(tab), tab, collapse = ", "))

mode(counts) <- "integer"
keep <- rowSums(counts >= 10) >= min_donors          # expression filter
counts <- counts[keep, , drop = FALSE]
if (nrow(counts) < 2) empty_out("too few genes after filtering")

res_df <- NULL
if (engine == "edger") {
  suppressMessages(library(edgeR))
  y <- DGEList(counts)
  y <- calcNormFactors(y)
  mm <- model.matrix(design, coldata)
  y <- estimateDisp(y, mm)
  fit <- glmQLFit(y, mm)
  coef <- paste0(grp, tst)
  qlf <- glmQLFTest(fit, coef = coef)
  tt <- topTags(qlf, n = Inf)$table
  res_df <- data.frame(gene = rownames(tt), baseMean = 2^(tt$logCPM), log2FoldChange = tt$logFC,
                       lfcSE = NA, stat = tt$F, pvalue = tt$PValue, padj = tt$FDR)
} else {
  suppressMessages(library(DESeq2))
  dds <- DESeqDataSetFromMatrix(counts, colData = coldata, design = design)
  # robust dispersion fit: parametric -> mean -> graceful empty (tiny data can break locfit)
  dds <- tryCatch(
    DESeq(dds, quiet = TRUE),
    error = function(e) tryCatch(DESeq(dds, fitType = "mean", quiet = TRUE),
                                 error = function(e2) NULL)
  )
  if (is.null(dds)) empty_out("DESeq dispersion fit failed on small/degenerate data")
  res <- results(dds, contrast = c(grp, tst, ref))   # keeps the Wald 'stat' for ranking
  # borrow the shrunken log2FC for effect-size ranking, but retain stat/pvalue/padj from `res`
  shr <- tryCatch(lfcShrink(dds, contrast = c(grp, tst, ref), res = res, type = "ashr"),
                  error = function(e) NULL)
  if (!is.null(shr)) res$log2FoldChange <- shr$log2FoldChange
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
}

res_df$contrast <- contrast[["name"]]
res_df$status <- "ok"
res_df <- res_df[order(res_df$padj), ]
write.table(res_df, out_f, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("DE done: %d genes, %d padj<0.05", nrow(res_df),
                sum(res_df$padj < 0.05, na.rm = TRUE)))

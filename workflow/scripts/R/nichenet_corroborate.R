# OPTIONAL — full NicheNet ligand-activity corroboration (alternative to the offline
# python corroboration in ccc_differential.py). Not on the default DAG because nichenetr is
# GitHub-only and downloads a large prior model. To use it, wire a rule that calls this
# script with the differential-CCC ligands and the receiver DE table.
#
# Install (pinned) once into the ccc_r env:
#   remotes::install_github("saeyslab/nichenetr", ref = "v2.0.4")
#
# Usage (standalone):
#   Rscript nichenet_corroborate.R <ligands.txt> <receiver_de.tsv> <out.tsv>
suppressMessages({
  library(nichenetr)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
ligands <- readLines(args[[1]])
de <- read.delim(args[[2]])
out <- args[[3]]

# NicheNet prior models (download once; cache under resources/nichenet/)
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
geneset_oi <- de %>% filter(log2FoldChange > 0.25, padj < 0.25) %>% pull(gene)
background <- de %>% pull(gene)

activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = intersect(ligands, rownames(ligand_target_matrix))
)
write.table(activities, out, sep = "\t", quote = FALSE, row.names = FALSE)

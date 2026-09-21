# ===========================================================================
#  27_scrna_mediator_expression.R
#
#  Purpose
#  -------
#  Final visualisation step: show where the mediator proteins identified by
#  the mediation analysis (step 21) are expressed across the single-cell
#  atlas, and confirm the identity of the main lineages.
#
#  Input
#  -----
#    output/objects/GSE161529_seurat_celltype.rdata   object: scRNA
#    output/objects/mediation_results.rdata           sig_mediation
#
#  Output
#  ------
#    output/figures/27_dotplot_mediators_cancer_vs_normal.pdf
#    output/figures/27_umap_mediator_genes.pdf
#    output/figures/27_umap_lineage_markers.pdf
#    output/tables/27_mediator_genes_tested.csv
#    output/tables/27_mean_expression_mediators.csv
#
#  Note on the genes plotted
#  -------------------------
#  The mediator list is intersected with the genes actually measured in this
#  dataset.  plot1cell::complex_dotplot_multiple() errors out on features
#  that are absent from the assay, and the released object was filtered with
#  min.cells = 3, so a small number of mediators may have been dropped.
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(RColorBrewer)
  library(plot1cell)
  library(scRNAtoolVis)
})

source("config/config.R")

load(file.path(OUT_OBJECTS, "GSE161529_seurat_celltype.rdata"))   # scRNA
load(file.path(OUT_OBJECTS, "mediation_results.rdata"))           # sig_mediation

# ---- 1. Case / control label ---------------------------------------------
scRNA$newident <- ifelse(scRNA$group %in% SC_GROUP_LABELS$cancer,
                         "Cancer", "Normal")
scRNA$newident <- factor(scRNA$newident, levels = c("Normal", "Cancer"))

message("Cells per group:")
print(table(scRNA$newident, scRNA$group))

# ---- 2. Mediator genes present in the dataset -----------------------------
mediators      <- unique(sig_mediation$mediator)
mediators_used <- intersect(mediators, rownames(scRNA))

message(sprintf("Mediators: %d reported, %d measured in this dataset",
                length(mediators), length(mediators_used)))
if (length(setdiff(mediators, mediators_used))) {
  message("  not measured: ", paste(setdiff(mediators, mediators_used), collapse = ", "))
}

write_table(data.frame(mediator = mediators,
                       measured = mediators %in% mediators_used),
            "27_mediator_genes_tested.csv")

# ---- 3. Dot plot, cancer versus normal -----------------------------------
strip_col <- c("#a6cee3", "#1f78b4", "#b2df8a", "#e31a1c", "#fdbf6f",
               "#cab2d6", "#62B197", "#E16E6D", "#9392BE", "#D0E7ED",
               "#D5E4A8", "#6a3d9a", "#4197d8", "#E79397", "#f8c120",
               "#CC88B0", "#413496", "#083356", "#ce1554", "#223e9c",
               "#b12b23", "#aebea6", "#edae11", "#0f6657")
strip_col <- rep_len(strip_col, length(mediators_used))

complex_dotplot_multiple(scRNA,
                         feature      = mediators_used,
                         group        = "newident",
                         color.palette = c("#1F78B4", "#9E0142"),
                         strip.color  = strip_col)

ggsave(file.path(OUT_FIGURES, "27_dotplot_mediators_cancer_vs_normal.pdf"),
       height = 6.5, width = 14)

# ---- 4. Mean expression table (same contrast, numeric form) --------------
mean_expr <- AverageExpression(scRNA, features = mediators_used,
                               group.by = "newident", slot = "data")$RNA
mean_expr <- as.data.frame(mean_expr)
mean_expr <- cbind(gene = rownames(mean_expr), mean_expr)

write_table(mean_expr, "27_mean_expression_mediators.csv")

# ---- 5. Spatial expression of selected mediators -------------------------
feature_umap <- function(features, file, width, height, facet = "newident") {
  p <- featureCornerAxes(object = scRNA, reduction = "umap",
                         groupFacet = facet,
                         relLength = 0.5, relDist = 0.2,
                         features = features)
  ggsave(file, p, width = width, height = height)
}

# FGF2 and OLR1 are the two mediators highlighted in the manuscript
sel <- intersect(c("FGF2", "OLR1"), rownames(scRNA))
if (length(sel)) {
  feature_umap(sel, file.path(OUT_FIGURES, "27_umap_mediator_genes.pdf"),
               width = 12, height = 9)
}

# ---- 6. Lineage identity check ------------------------------------------
lineage <- intersect(c("CD3E", "EPCAM", "PECAM1", "CD68", "LUM"),
                     rownames(scRNA))

feature_umap(lineage, file.path(OUT_FIGURES, "27_umap_lineage_markers.pdf"),
             width = 14, height = 4.5, facet = NULL)

message("Step 27 complete.")

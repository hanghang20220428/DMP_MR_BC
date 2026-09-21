# ===========================================================================
#  02_dmp_heatmap.R
#
#  Purpose
#  -------
#  Heat map of the differentially methylated positions (DMPs) identified in
#  step 01, with samples ordered by group so that the case/control contrast
#  is visible.  Beta values are row-scaled (z-scored across samples).
#
#  Input
#  -----
#    output/objects/GSE104942_norm.rdata       matData, clinical
#    output/objects/GSE104942_DMP_Diff.rdata   dmpDiff
#
#  Output
#  ------
#    output/figures/02_DMP_heatmap.pdf
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

source("config/config.R")

load(file.path(OUT_OBJECTS, "GSE104942_norm.rdata"))      # matData, clinical
load(file.path(OUT_OBJECTS, "GSE104942_DMP_Diff.rdata"))  # dmpDiff

# ---- Subset the beta matrix to the significant CpGs -----------------------
dmp_exp <- matData[rownames(dmpDiff), ]

# Order samples so that all controls come before all cases.
clinical <- clinical %>% dplyr::arrange(desc(Group))

annotation_col    <- data.frame(group = clinical$Group)
rownames(annotation_col) <- clinical$Sample

dmp_exp <- dmp_exp[, rownames(annotation_col)]

# ---- Draw the heat map ----------------------------------------------------
p <- pheatmap(
  dmp_exp,
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  annotation_col = annotation_col,
  annotation_legend = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  scale         = "row",
  color         = colorRampPalette(c("#4197d8", "white", "#E41A1C"))(50),
  fontsize      = 10,
  silent        = TRUE
)

pdf(file.path(OUT_FIGURES, "02_DMP_heatmap.pdf"), width = 10, height = 10)
grid::grid.draw(p$gtable)
dev.off()

message("Heat map written to output/figures/02_DMP_heatmap.pdf")

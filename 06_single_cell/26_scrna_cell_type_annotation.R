# ===========================================================================
#  26_scrna_cell_type_annotation.R
#
#  Purpose
#  -------
#  Cluster marker discovery and cell-type annotation of the GSE161529 atlas
#  produced by step 25.  Cluster identities were assigned manually from the
#  expression of canonical lineage markers, as is standard practice for
#  droplet-based breast-tissue atlases.
#
#  Input
#  -----
#    output/objects/GSE161529_seurat_clustered.rdata     object: scRNA
#
#  Output
#  ------
#    output/objects/GSE161529_seurat_celltype.rdata      object: scRNA
#    output/tables/26_cluster_markers_all.csv
#    output/tables/26_cluster_markers_significant.csv
#    output/figures/26_marker_heatmap.pdf
#    output/figures/26_marker_dotplot_clusters.pdf
#    output/figures/26_umap_celltype.pdf
#    output/figures/26_umap_celltype_by_group.pdf
#
#  Marker statistics (FindAllMarkers)
#      only.pos        = TRUE     only up-regulated markers are kept
#      min.pct         = 0.25     gene must be detected in >= 25% of the
#                                 cells in one of the two groups
#      logfc.threshold = 0.25     minimum log fold change (natural log)
#      test.use        = "wilcox" Seurat default (Wilcoxon rank-sum test)
#      a marker is reported as significant when |avg_log2FC| > 0.5 and
#      Bonferroni-adjusted P < 0.05
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(RColorBrewer)
})

source("config/config.R")

load(file.path(OUT_OBJECTS, "GSE161529_seurat_clustered.rdata"))   # scRNA

Idents(scRNA) <- "seurat_clusters"

# ===========================================================================
# 1. Cluster markers
# ===========================================================================

markers <- FindAllMarkers(scRNA,
                          only.pos        = SC_QC_PARAMS$marker_only_pos,
                          min.pct         = SC_QC_PARAMS$marker_min_pct,
                          logfc.threshold = SC_QC_PARAMS$marker_logfc_thresh)

write_table(markers, "26_cluster_markers_all.csv")

sig_markers <- markers[abs(markers$avg_log2FC) > SC_QC_PARAMS$marker_logfc_keep &
                         markers$p_val_adj < SC_QC_PARAMS$marker_padj_keep, , drop = FALSE]

write_table(sig_markers, "26_cluster_markers_significant.csv")
message(sprintf("Markers: %d detected, %d significant", nrow(markers), nrow(sig_markers)))

# Top 5 markers per cluster, used for the heat map
top5 <- markers %>%
  group_by(cluster) %>%
  top_n(5, wt = avg_log2FC)

aver <- AverageExpression(scRNA, return.seurat = TRUE)

p_heatmap <- DoHeatmap(aver, features = top5$gene,
                       size = 3, draw.lines = FALSE) +
  scale_fill_gradientn(colors = c("#5E4FA2", "white", "#f8c120"))

ggsave(file.path(OUT_FIGURES, "26_marker_heatmap.pdf"),
       p_heatmap, width = 45, height = 45, limitsize = FALSE)

# ===========================================================================
# 2. Canonical lineage markers
# ===========================================================================
#  Reference panel used to inspect every cluster before assigning identities.

ref_gene <- sort(c(
  # epithelial / tumour
  "EPCAM", "KRT18", "KRT19", "MKI67", "CDK1",
  # T cells
  "CD3E", "CD4", "CD8A", "IL7R", "CCL5", "FOXP3", "IL2RA", "CD40LG", "KLRB1",
  # NK cells
  "NCR1", "NKG7",
  # B / plasma cells
  "MS4A1", "CD19", "CD79A", "JCHAIN", "IGHG1", "IGLL5",
  # myeloid
  "LYZ", "CD68", "ITGAX", "FCGR1A", "FCGR3A", "C1QA", "APOC1", "SPP1",
  "S100A12", "FCN1", "S100A9", "CD14", "CD163", "MARCO",
  # dendritic cells
  "FCER1A", "CD1C", "CLEC9A", "LILRA4", "CLEC4C",
  # endothelial
  "PECAM1", "VWF",
  # fibroblasts / smooth muscle
  "COL1A1", "LUM", "PDGFRA", "ACTA2", "MYH11", "MYLK", "FGF7",
  # mast cells
  "CPA3", "KIT",
  # immunotherapy target
  "CD274"))

ref_gene <- ref_gene[ref_gene %in% rownames(scRNA)]

p_dot <- DotPlot(scRNA, features = ref_gene, group.by = "seurat_clusters") +
  coord_flip() +
  theme_bw() +
  scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral")))

ggsave(file.path(OUT_FIGURES, "26_marker_dotplot_clusters.pdf"),
       p_dot, height = 12, width = 8.5)

# ===========================================================================
# 3. Manual cluster -> cell type assignment
# ===========================================================================
#  Cluster numbers refer to FindClusters(resolution = 1) on this dataset.
#  Adjust the mapping below if a different Seurat version changes the cluster
#  numbering, and always re-check the dot plot above before doing so.

cluster_map <- list(
  Epithelial_cells  = c(1, 2, 3, 5, 6, 10, 12, 15, 18, 19, 20, 21, 24, 25, 26, 27, 30, 33),
  Fibroblasts       = c(7, 9, 14, 17, 22, 31),
  T_cells           = c(0, 4, 23),
  Macrophages       = c(11, 13, 28, 32),
  Endothelial_cells = c(8, 34),
  B_cells           = c(16),
  Dendritic_cells   = c(29)
)

scRNA$Celltype <- NA_character_
for (ct in names(cluster_map)) {
  cl <- as.character(cluster_map[[ct]])
  scRNA$Celltype[scRNA$seurat_clusters %in% cl] <- ct
}

unassigned <- sum(is.na(scRNA$Celltype))
if (unassigned > 0) warning(unassigned, " cells were not assigned to a cell type.")

scRNA$Celltype <- factor(scRNA$Celltype, levels = names(cluster_map))
Idents(scRNA) <- scRNA$Celltype

message("Cells per cell type:")
print(table(scRNA$Celltype))

# ===========================================================================
# 4. Plots and save
# ===========================================================================

p_ct <- DimPlot(scRNA, reduction = "umap", label = TRUE,
                raster = FALSE, pt.size = 0.5) + scale_color_nejm()

ggsave(file.path(OUT_FIGURES, "26_umap_celltype.pdf"),
       p_ct, width = 9.5, height = 7)

p_ct_group <- DimPlot(scRNA, reduction = "umap", label = TRUE,
                      split.by = "group", raster = FALSE, pt.size = 0.5) +
  scale_color_nejm()

ggsave(file.path(OUT_FIGURES, "26_umap_celltype_by_group.pdf"),
       p_ct_group, width = 22, height = 6)

# Cell-type composition per group
ct_table <- as.data.frame.matrix(table(scRNA$group, scRNA$Celltype))
ct_table <- cbind(group = rownames(ct_table), ct_table)
write_table(ct_table, "26_celltype_composition_by_group.csv")

save_rdata(scRNA, file = "GSE161529_seurat_celltype.rdata")

message("Step 26 complete.")

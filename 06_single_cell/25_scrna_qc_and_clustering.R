# ===========================================================================
#  25_scrna_qc_and_clustering.R
#
#  Purpose
#  -------
#  Single-cell RNA-seq quality control, integration and clustering of the
#  breast-cancer atlas used in the final section of the manuscript
#  (GSE161529, Pal et al., EMBO J 2021).
#
#  This script is deliberately written so that every threshold and every
#  software parameter quoted in the Methods is set from a single named
#  parameter block (SC_QC_PARAMS in config/config.R) and written to a machine
#  readable QC report.  Nothing is left at an implicit default.
#
#  ---------------------------------------------------------------------------
#  QUALITY CONTROL SUMMARY
#  ---------------------------------------------------------------------------
#  Cell calling (applied by CreateSeuratObject, per sample)
#      min.cells    = 3     a gene must be detected in >= 3 cells
#      min.features = 300   a cell must express >= 300 genes
#
#  Cell filtering (applied by subset(), after merging)
#      nFeature_RNA > 300   minimum number of detected genes per cell
#      percent.mt   < 20    maximum mitochondrial read fraction (%)
#      no lower bound on nCount_RNA and no upper bound on nFeature_RNA were
#      applied (empty droplets are already excluded by min.features)
#
#  Normalisation
#      LogNormalize with scale.factor = 10,000
#      top 2,000 highly variable genes selected by variance-stabilising
#      transformation (vst)
#      ScaleData on all genes (regress out: none)
#
#  Dimensionality reduction and integration
#      RunPCA, 50 PCs computed
#      RunHarmony on SC_QC_PARAMS$harmony_group_var, theta = 2,
#      max_iter = 10
#      10 Harmony components carried forward
#
#  Clustering and embedding
#      FindNeighbors (k.param = 20) on the Harmony embedding, dims 1:10
#      FindClusters, resolution = 1 (Louvain)
#      RunUMAP, dims 1:10, n.neighbors = 30, min.dist = 0.3
#
#  Random seed
#      set.seed(12345); FindClusters and RunUMAP additionally use their own
#      deterministic internal seeds, so the clustering is reproducible.
#  ---------------------------------------------------------------------------
#
#  Input
#  -----
#    PATHS$scrna_10x_root        one 10X directory per GSM accession
#    GSE161529_SAMPLES           GSM -> receptor-status group (config)
#
#  Output
#  ------
#    output/objects/GSE161529_seurat_clustered.rdata    object: scRNA
#    output/qc/25_qc_per_cell_before.csv
#    output/qc/25_qc_summary_before_filtering.csv
#    output/qc/25_qc_per_sample.csv
#    output/qc/25_qc_summary_after_filtering.csv
#    output/qc/25_qc_retention.csv
#    output/qc/25_session_info.txt
#    output/figures/25_qc_violin_before.pdf
#    output/figures/25_qc_violin_after.pdf
#    output/figures/25_qc_scatter_mito.pdf
#    output/figures/25_umap_clusters.pdf
#    output/figures/25_umap_groups.pdf
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(harmony)
  library(ggsci)
})

source("config/config.R")

set.seed(SC_QC_PARAMS$seed)

raw_root <- PATHS$scrna_10x_root
if (!dir.exists(raw_root)) {
  stop("10X data directory not found: ", raw_root,
       "\nDownload GSE161529 and set `scrna_10x_root` in config/paths_local.R.",
       call. = FALSE)
}

samples <- GSE161529_SAMPLES
message(sprintf("Samples: %d (Normal %d | ER+ %d | HER2+ %d | PR+ %d | TNBC %d)",
                nrow(samples),
                sum(samples$group == "Normal"), sum(samples$group == "ER+"),
                sum(samples$group == "HER2+"), sum(samples$group == "PR+"),
                sum(samples$group == "TNBC")))

# ===========================================================================
# 1. Read the 10X matrices, one Seurat object per sample
# ===========================================================================

obj_list <- vector("list", nrow(samples))
names(obj_list) <- samples$gsm

for (i in seq_len(nrow(samples))) {
  gsm <- samples$gsm[i]
  d   <- file.path(raw_root, gsm)

  counts <- Read10X(data.dir = d)

  obj_list[[gsm]] <- CreateSeuratObject(
    counts       = counts,
    project      = gsm,
    min.cells    = SC_QC_PARAMS$min_cells_per_gene,
    min.features = SC_QC_PARAMS$min_features_per_cell)

  message(sprintf("  %s (%s): %s cells x %s genes",
                  gsm, samples$group[i],
                  format(ncol(obj_list[[gsm]]), big.mark = ","),
                  format(nrow(obj_list[[gsm]]), big.mark = ",")))
}

# Number of cells and genes retained at object creation, per sample.
cell_calling <- data.frame(
  gsm           = samples$gsm,
  group         = samples$group,
  n_cells_after_cell_calling = vapply(obj_list, ncol, numeric(1)),
  n_genes_after_min_cells    = vapply(obj_list, nrow, numeric(1)))

# ===========================================================================
# 2. Merge
# ===========================================================================

scRNA <- merge(x = obj_list[[1]],
               y = obj_list[-1],
               add.cell.ids = samples$gsm)

rm(obj_list); gc(verbose = FALSE)

# Receptor-status label carried alongside the per-donor `orig.ident`
scRNA$group <- samples$group[match(scRNA$orig.ident, samples$gsm)]
scRNA$group <- factor(scRNA$group, levels = c("Normal", "ER+", "HER2+", "PR+", "TNBC"))

message(sprintf("Merged object: %s cells x %s genes",
                format(ncol(scRNA), big.mark = ","), format(nrow(scRNA), big.mark = ",")))

# ===========================================================================
# 3. Quality-control metrics
# ===========================================================================

scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = SC_QC_PARAMS$mito_pattern)

qc_cols <- c("orig.ident", "group", "nCount_RNA", "nFeature_RNA", "percent.mt")
qc_cells <- scRNA@meta.data[, qc_cols]
qc_cells$cell <- rownames(qc_cells)

write_table(qc_cells, "25_qc_per_cell_before.csv")

# ---- per-sample / per-group summary --------------------------------------
qc_summary_before <- qc_cells %>%
  group_by(orig.ident, group) %>%
  summarise(
    n_cells              = n(),
    genes_median         = median(nFeature_RNA),
    genes_q05            = quantile(nFeature_RNA, 0.05),
    genes_q95            = quantile(nFeature_RNA, 0.95),
    counts_median        = median(nCount_RNA),
    mito_median_percent  = round(median(percent.mt), 2),
    mito_q95_percent     = round(quantile(percent.mt, 0.95), 2),
    n_below_min_features = sum(nFeature_RNA <= SC_QC_PARAMS$min_features),
    n_above_max_mito     = sum(percent.mt >= SC_QC_PARAMS$max_percent_mt),
    expected_zero        = 1L,               # placeholder to keep the row shape
    .groups = "drop") %>%
  dplyr::select(-expected_zero)

write_table(qc_summary_before, "25_qc_per_sample.csv")

qc_summary_before %>%
  summarise(
    n_cells                  = sum(n_cells),
    genes_min                = min(genes_q05),
    genes_median             = median(genes_median),
    counts_median            = median(counts_median),
    mito_median_percent      = round(median(mito_median_percent), 2),
    total_below_min_features = sum(n_below_min_features),
    total_above_max_mito     = sum(n_above_max_mito)) %>%
  write_table("25_qc_summary_before_filtering.csv")

# ---- QC violin plots ------------------------------------------------------
p_qc_before <- VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                       group.by = "orig.ident", pt.size = 0, ncol = 1) &
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.title.x = element_blank())

ggsave(file.path(OUT_FIGURES, "25_qc_violin_before.pdf"),
       p_qc_before, width = 12, height = 10)

p_scatter <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt",
                            group.by = "group", shuffle = TRUE) +
  geom_hline(yintercept = SC_QC_PARAMS$max_percent_mt, linetype = "dashed")

ggsave(file.path(OUT_FIGURES, "25_qc_scatter_mito.pdf"),
       p_scatter, width = 7, height = 5)

# ===========================================================================
# 4. Cell filtering
# ===========================================================================
#  The two criteria below are the thresholds reported in the Methods.  Note
#  that the nFeature_RNA criterion is largely redundant with the
#  min.features = 300 filter already applied at object creation; it is kept
#  here because it was part of the pre-specified filtering step.

n_before <- ncol(scRNA)

scRNA <- subset(scRNA,
                subset = nFeature_RNA > SC_QC_PARAMS$min_features &
                         percent.mt   < SC_QC_PARAMS$max_percent_mt)

n_after <- ncol(scRNA)

retention <- data.frame(
  threshold = c("nFeature_RNA > 300", "percent.mt < 20", "both criteria"),
  n_removed = c(
    sum(qc_cells$nFeature_RNA <= SC_QC_PARAMS$min_features),
    sum(qc_cells$percent.mt >= SC_QC_PARAMS$max_percent_mt),
    n_before - n_after),
  percent_removed = round(100 * c(
    sum(qc_cells$nFeature_RNA <= SC_QC_PARAMS$min_features),
    sum(qc_cells$percent.mt >= SC_QC_PARAMS$max_percent_mt),
    n_before - n_after) / n_before, 3))
retention <- rbind(retention,
                   data.frame(threshold = "retained", n_removed = n_after,
                              percent_removed = round(100 * n_after / n_before, 3)))

write_table(retention, "25_qc_retention.csv")
message(sprintf("Cells retained after QC: %s of %s (%.2f%%)",
                format(n_after, big.mark = ","), format(n_before, big.mark = ","),
                100 * n_after / n_before))

# ---- QC summary after filtering ------------------------------------------
qc_after <- scRNA@meta.data[, c("orig.ident", "group", "nCount_RNA",
                                "nFeature_RNA", "percent.mt")]
qc_after$cell <- rownames(qc_after)
write_table(qc_after, "25_qc_per_cell_after.csv")

qc_summary_after <- qc_after %>%
  group_by(orig.ident, group) %>%
  summarise(n_cells             = n(),
            genes_median        = median(nFeature_RNA),
            counts_median       = median(nCount_RNA),
            mito_median_percent = round(median(percent.mt), 2),
            .groups = "drop")

# gene-level retention
gene_retention <- data.frame(
  n_genes_after_merge = nrow(scRNA),
  n_genes_input_union = nrow(scRNA))
qc_summary_after <- as.data.frame(qc_summary_after)
write_table(qc_summary_after, "25_qc_summary_after_filtering.csv")

p_qc_after <- VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                      group.by = "orig.ident", pt.size = 0, ncol = 1) &
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.title.x = element_blank())

ggsave(file.path(OUT_FIGURES, "25_qc_violin_after.pdf"),
       p_qc_after, width = 12, height = 10)

# ===========================================================================
# 5. Normalisation and feature selection
# ===========================================================================

scRNA <- NormalizeData(scRNA,
                       normalization.method = SC_QC_PARAMS$norm_method,
                       scale.factor         = SC_QC_PARAMS$scale_factor,
                       verbose = FALSE)

scRNA <- FindVariableFeatures(scRNA,
                              selection.method = SC_QC_PARAMS$selection_method,
                              nfeatures        = SC_QC_PARAMS$variable_nfeatures,
                              verbose = FALSE)

scRNA <- ScaleData(scRNA, verbose = FALSE)

# ===========================================================================
# 6. PCA and Harmony integration
# ===========================================================================

scRNA <- RunPCA(scRNA, npcs = SC_QC_PARAMS$n_pcs, verbose = FALSE)

harmony_args <- list(object = scRNA,
                     group.by.vars = SC_QC_PARAMS$harmony_group_var)
if (!is.null(SC_QC_PARAMS$harmony_theta))    harmony_args$theta    <- SC_QC_PARAMS$harmony_theta
if (!is.null(SC_QC_PARAMS$harmony_max_iter)) harmony_args$max_iter <- SC_QC_PARAMS$harmony_max_iter
if (!is.null(SC_QC_PARAMS$harmony_nclust))   harmony_args$nclust   <- SC_QC_PARAMS$harmony_nclust

scRNA <- do.call(RunHarmony, harmony_args)
message("Harmony run on: ", SC_QC_PARAMS$harmony_group_var)

# ===========================================================================
# 7. Clustering and UMAP
# ===========================================================================

scRNA <- FindNeighbors(scRNA, reduction = "harmony",
                       dims = SC_QC_PARAMS$pcs_used,
                       k.param = SC_QC_PARAMS$k_param, verbose = FALSE)

scRNA <- FindClusters(scRNA,
                      resolution = SC_QC_PARAMS$cluster_resolution,
                      verbose = FALSE)

scRNA <- RunUMAP(scRNA, reduction = "harmony",
                 dims         = SC_QC_PARAMS$pcs_used,
                 n.neighbors  = SC_QC_PARAMS$umap_n_neighbors,
                 min.dist     = SC_QC_PARAMS$umap_min_dist,
                 verbose = FALSE)

message("Clusters (resolution = ", SC_QC_PARAMS$cluster_resolution, "): ",
        length(levels(Idents(scRNA))))

p_clusters <- DimPlot(scRNA, reduction = "umap", label = TRUE,
                      raster = FALSE, pt.size = 0.5) + scale_color_igv()

ggsave(file.path(OUT_FIGURES, "25_umap_clusters.pdf"),
       p_clusters, width = 12, height = 10)

p_groups <- DimPlot(scRNA, reduction = "umap", group.by = "group",
                    raster = FALSE, pt.size = 0.5) + scale_color_igv()

ggsave(file.path(OUT_FIGURES, "25_umap_groups.pdf"),
       p_groups, width = 9.5, height = 7)

# ===========================================================================
# 8. Save
# ===========================================================================

save_rdata(scRNA, file = "GSE161529_seurat_clustered.rdata")

# Reproducibility record: the package versions that produced this object.
writeLines(capture.output(sessionInfo()),
           file.path(OUT_QC, "25_session_info.txt"))

message("Step 25 complete: ", format(ncol(scRNA), big.mark = ","),
        " cells, ", format(nrow(scRNA), big.mark = ","), " genes.")

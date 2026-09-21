# ===========================================================================
#  18_phewas_circos_plot.R
#
#  Purpose
#  -------
#  Circular ("circos") heat map of the phenome-wide MR results: for each
#  causal CpG, the 20 most significant FinnGen endpoints are drawn as
#  concentric tracks coloured by -log10(P), with endpoints grouped by
#  FinnGen disease category.
#
#  Input
#  -----
#    output/objects/PheWAS_FinnGen.rdata    AllOR
#    PATHS$finngen_manifest                 endpoint names and categories
#
#  Output
#  ------
#    output/figures/18_PheWAS_circos.pdf
#    output/tables/18_PheWAS_circos_input.csv
#
#  Layout
#  ------
#  One track per causal CpG (columns 1-3 of the input matrix); the outermost
#  track shows the FinnGen disease category of each endpoint.
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(tibble)
  library(circlize)
  library(ComplexHeatmap)
})

source("config/config.R")

load(file.path(OUT_OBJECTS, "PheWAS_FinnGen.rdata"))   # AllOR

finninfo <- fread(PATHS$finngen_manifest, data.table = FALSE)

AllOR$outcome <- AllOR$id.outcome
data0 <- merge(AllOR, finninfo, by.x = "outcome", by.y = "phenocode")

# ---- 1. Select the 20 most significant endpoints per CpG ------------------
data <- data0 %>%
  dplyr::group_by(exposure) %>%
  dplyr::top_n(20, wt = desc(pval)) %>%
  dplyr::select(exposure, phenotype, pval, category) %>%
  dplyr::ungroup()

data$category <- sapply(strsplit(data$category, split = "[(]"), "[", 1)
data$exposure <- as.factor(data$exposure)
data$category <- as.factor(data$category)

data <- data %>%
  dplyr::arrange(pval) %>%
  dplyr::distinct(phenotype, .keep_all = TRUE) %>%
  pivot_wider(names_from = exposure, values_from = pval) %>%
  column_to_rownames("phenotype") %>%
  dplyr::select(-category, category) %>%
  dplyr::arrange(category)

n_cpg <- ncol(data) - 1L                       # last column is `category`
data[, seq_len(n_cpg)] <- -log10(data[, seq_len(n_cpg)])

write_table(cbind(phenotype = rownames(data), data), "18_PheWAS_circos_input.csv")

# ---- 2. Colour mapping ----------------------------------------------------
col_n <- c("#6a73cf", "#edd064", "#0eb0c8", "#cdb4d7", "#f68b6a",
           "#083356", "#ce1554", "#FFFFBF", "#3288BD", "#e1abbc",
           "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
           "#9E0142", "#fdbf6f", "#ff7f00", "#cab2d6", "#62B197",
           "#E16E6D", "#9392BE", "#D0E7ED", "#D5E4A8", "#6a3d9a",
           "#4197d8", "#E79397", "#f8c120", "#CC88B0", "#a1d5b9",
           "#50b688", "#962d20", "#5b2d90", "#ffda67", "#339900FF",
           "#FF1463FF", "#E7BA52FF", "#990080FF", "#008099FF")
names(col_n) <- unique(as.character(data$category))

col_fun <- c(
  lapply(seq_len(n_cpg), function(i) colorRamp2(c(1, 5), c("#0eb0c8", "#ce1554"))),
  list(col_n)
)

# ---- 3. Draw --------------------------------------------------------------
circos.clear()
pdf(file.path(OUT_FIGURES, "18_PheWAS_circos.pdf"), height = 25, width = 25)
circos.par(gap.degree = 60, start.degree = 30, track.margin = c(0.001, 0.001))

for (i in seq_len(n_cpg + 1L)) {

  data_tmp <- as.matrix(data[, i])
  if (i == 1) rownames(data_tmp) <- rownames(data)
  colnames(data_tmp) <- colnames(data)[i]

  if (i <= n_cpg) {
    circos.heatmap(data_tmp, col = col_fun[[i]],
                   rownames.side = "outside", rownames.cex = 0.8,
                   cluster = FALSE, cell.border = "white",
                   show.sector.labels = TRUE, track.height = 0.025)
  } else {
    circos.heatmap(data_tmp, col = col_fun[[i]],
                   rownames.side = "outside",
                   cluster = TRUE, track.height = 0.02)
  }
}

lgd1 <- Legend(title = "-log10(p)", border = "black",
               grid_height = unit(8.5, "mm"), legend_width = unit(35, "mm"),
               at = c(1, 4), title_position = "topcenter",
               col_fun = col_fun[[1]], direction = "horizontal")
draw(packLegend(lgd1), x = unit(0.62, "npc"), y = unit(0.7, "npc"))

# Track labels for the causal CpGs, placed outside the first sector
m <- -2.5; n <- 1.4; p <- 12
circos.track(track.index = get.current.track.index(),
             panel.fun = function(x, y) {
               if (CELL_META$sector.numeric.index == 1) {
                 circos.rect(CELL_META$cell.xlim[2] + convert_x(p, "mm"), -0.5,
                             CELL_META$cell.xlim[2] + convert_x(p, "mm"), 11,
                             border = NA)
                 for (k in seq_len(n_cpg)) {
                   circos.text(CELL_META$cell.xlim[2] + convert_x(p, "mm"),
                               m + (k + 2) * n, colnames(data)[k],
                               cex = 1, facing = "inside")
                 }
               }
             }, bg.border = NA)

lgd <- Legend(labels = names(col_n), legend_gp = gpar(fill = col_n),
              title = "Category", ncol = 1, row_gap = unit(1, "mm"))
draw(lgd, x = unit(0.515, "npc"), y = unit(0.5, "npc"))

dev.off()
circos.clear()

message("Step 18 complete.")

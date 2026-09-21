# ===========================================================================
#  23_cpg_methylation_boxplot.R
#
#  Purpose
#  -------
#  Box plots of the methylation beta value of each causal CpG in breast-cancer
#  cases versus unaffected controls, in the discovery cohort (GSE104942).
#  This documents the differential methylation that motivated the MR screen.
#
#  Input
#  -----
#    output/objects/GSE104942_norm.rdata       matData, clinical
#    output/objects/mediation_results.rdata    sig_mediation
#
#  Output
#  ------
#    output/figures/23_cpg_methylation_boxplot.pdf
#    output/tables/23_cpg_methylation_by_group.csv
#
#  Statistics
#  ----------
#  Two-sided Student's t-test per CpG, significance symbols via ggsignif
#  (ns P >= 0.05, * P < 0.05, ** P < 0.01, *** P < 0.001, **** P < 0.0001).
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(tibble)
  library(ggpubr)
})

source("config/config.R")

load(file.path(OUT_OBJECTS, "GSE104942_norm.rdata"))      # matData, clinical
load(file.path(OUT_OBJECTS, "mediation_results.rdata"))   # sig_mediation

causal_cpg <- unique(sig_mediation$exposure)
message("Causal CpGs: ", paste(causal_cpg, collapse = ", "))

# ---- 1. Beta values in long format ---------------------------------------
cpg_beta <- matData[causal_cpg, , drop = FALSE] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample")

clinical <- clinical %>% dplyr::arrange(desc(Group))

plot_data <- left_join(clinical, cpg_beta, by = "Sample") %>%
  pivot_longer(cols = all_of(causal_cpg), names_to = "CpG", values_to = "Beta")

plot_data$CpG   <- factor(plot_data$CpG, levels = causal_cpg)
plot_data$Group <- factor(plot_data$Group, levels = c("Control", "Cancer"))

write_table(plot_data %>% dplyr::select(Sample, Group, CpG, Beta),
            "23_cpg_methylation_by_group.csv")

# ---- 2. Pairwise comparisons ---------------------------------------------
comparisons <- combn(levels(plot_data$Group), 2, simplify = FALSE)

p <- ggplot(plot_data, aes(x = Group, y = Beta, fill = Group)) +
  facet_grid(. ~ CpG) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("#E69F00", "#413496")) +
  geom_jitter(size = 0.2, alpha = 0.6) +
  theme_bw() +
  stat_compare_means(comparisons = comparisons,
                     label = "p.signif", method = "t.test")

ggsave(file.path(OUT_FIGURES, "23_cpg_methylation_boxplot.pdf"),
       p, width = 5, height = 3.5)

message("Step 23 complete.")

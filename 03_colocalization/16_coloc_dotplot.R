# ===========================================================================
#  16_coloc_dotplot.R
#
#  Purpose
#  -------
#  Bubble plot summarising the colocalisation results for the causal CpGs
#  against both outcome GWAS.
#
#  Input
#  -----
#    output/objects/coloc_FinnGen.rdata    coloc_summary
#    output/objects/coloc_BCAC.rdata       coloc_summary
#
#  Output
#  ------
#    output/figures/16_coloc_dotplot.pdf
#    output/tables/16_coloc_summary_combined.csv
#
#  The bubble colour and size encode top_snp_PP.H4, the posterior probability
#  that the lead SNP of the CpG is the shared causal variant, for each CpG and
#  each outcome cohort.  This is the quantity reported in the manuscript.
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
})

source("config/config.R")

load(file.path(OUT_OBJECTS, "coloc_FinnGen.rdata"))   # coloc_summary
finn <- coloc_summary
load(file.path(OUT_OBJECTS, "coloc_BCAC.rdata"))      # coloc_summary
bcac <- coloc_summary

data <- rbind(finn, bcac)
data$outcome_short <- ifelse(grepl("FinnGen", data$outcome),
                             "BC (FinnGen R11)", "BC (BCAC 2020)")

write_table(data, "16_coloc_summary_combined.csv")

p <- ggplot(data, aes(x = outcome_short, y = cpg)) +
  geom_point(aes(size = top_snp_PP.H4, fill = top_snp_PP.H4),
             shape = 21, colour = "#4197d8") +
  geom_text(aes(label = sprintf("%.3f", top_snp_PP.H4)), size = 3) +
  scale_size_continuous(range = c(7, 9), guide = "none") +
  scale_fill_gradient(low = "#D0E7ED", high = "#1F78B4", guide = "none") +
  labs(x = NULL, y = NULL,
       title = "Colocalisation: posterior probability of the lead SNP") +
  theme_bw() +
  theme(panel.grid       = element_blank(),
        panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        axis.text        = element_text(size = 11),
        axis.text.y      = element_text(angle = 0, hjust = 1),
        title            = element_text(size = 13))

ggsave(file.path(OUT_FIGURES, "16_coloc_dotplot.pdf"),
       p, width = 8.5, height = 3.5)

message("Step 16 complete.")

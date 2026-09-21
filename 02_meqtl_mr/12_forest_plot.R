# ===========================================================================
#  12_forest_plot.R
#
#  Purpose
#  -------
#  Forest plot comparing the TSMR and GSMR estimates for CpGs that were
#  significant in both the FinnGen and BCAC screens (CpG set A), faceted by
#  method and outcome.
#
#  Input
#  -----
#    output/objects/TSMR_FinnGen_nohet.rdata     AllOR_nohet
#    output/objects/TSMR_BCAC_nohet.rdata        AllOR_nohet
#    output/objects/GSMR_FinnGen.rdata           AllOR
#    output/objects/GSMR_BCAC.rdata              AllOR
#
#  Output
#  ------
#    output/figures/12_MR_forest_plot.pdf
#    output/tables/12_MR_estimates_all_methods.csv
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(cowplot)
})

source("config/config.R")

keep_cols <- c("exposure", "outcome", "method", "pval", "b", "se",
               "or", "or_lci95", "or_uci95")

load(file.path(OUT_OBJECTS, "TSMR_FinnGen_nohet.rdata"))   # AllOR_nohet
Finn_TSMR <- AllOR_nohet[, keep_cols]
Finn_TSMR$method <- "TSMR"

load(file.path(OUT_OBJECTS, "TSMR_BCAC_nohet.rdata"))      # AllOR_nohet
BCAC_TSMR <- AllOR_nohet[, keep_cols]
BCAC_TSMR$method <- "TSMR"

load(file.path(OUT_OBJECTS, "GSMR_FinnGen.rdata"))         # AllOR
Finn_GSMR <- AllOR[AllOR$pval < 0.05, keep_cols]

load(file.path(OUT_OBJECTS, "GSMR_BCAC.rdata"))            # AllOR
BCAC_GSMR <- AllOR[AllOR$pval < 0.05, keep_cols]

# ---- CpGs detected by GSMR in both outcomes -------------------------------
same_cpg <- intersect(Finn_GSMR$exposure, BCAC_GSMR$exposure)
message("CpGs significant in both GSMR screens: ",
        paste(same_cpg, collapse = ", "))

allmr <- rbind(
  Finn_TSMR[Finn_TSMR$exposure %in% same_cpg, ],
  BCAC_TSMR[BCAC_TSMR$exposure %in% same_cpg, ],
  Finn_GSMR[Finn_GSMR$exposure %in% same_cpg, ],
  BCAC_GSMR[BCAC_GSMR$exposure %in% same_cpg, ]
)

allmr$outcome[allmr$outcome == "BC (FinnGen)"] <- "BC_FinnGenR11"

write_table(allmr, "12_MR_estimates_all_methods.csv")

# ---- Forest plot ----------------------------------------------------------
p <- ggplot(allmr, aes(y = exposure, x = or, colour = method)) +
  geom_errorbarh(aes(xmin = or_lci95, xmax = or_uci95), height = 0.2) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#4197d8", "#9E0142")) +
  geom_vline(xintercept = 1, linetype = "longdash") +
  theme_minimal_hgrid(10, rel_small = 1) +
  facet_wrap(~ method + outcome, ncol = 1) +
  labs(y = "", x = "OR (95% CI)",
       title = "Significant MR results between CpG and BC") +
  theme(legend.position = "none",
        strip.text      = element_text(size = 14, face = "bold"),
        axis.title.x    = element_text(size = 16, face = "bold"),
        axis.title.y    = element_text(size = 16, face = "bold"))

ggsave(file.path(OUT_FIGURES, "12_MR_forest_plot.pdf"),
       p, width = 8, height = 11)

message("Step 12 complete.")

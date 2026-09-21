# ===========================================================================
#  24_alluvial_plot.R
#
#  Purpose
#  -------
#  Alluvial (Sankey-style) diagram tracing each significant mediated pathway
#  through the four levels of the analysis:
#
#      instrument SNP  ->  CpG  ->  mediator protein  ->  outcome
#
#  Input
#  -----
#    output/objects/mediation_results.rdata   sig_mediation
#
#  Output
#  ------
#    output/figures/24_alluvial_mediation.pdf
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggalluvial)
})

source("config/config.R")

load(file.path(OUT_OBJECTS, "mediation_results.rdata"))   # sig_mediation

# Columns 1-4 of the mediation table are: SNP, exposure, mediator, outcome
df <- to_lodes_form(sig_mediation[, 1:4], axes = 1:4, id = "value")

col <- c("#7bc4e2", "#acd372", "#fbb05b", "#ed6ca4", "#cdb4d7", "#0eb0c8",
         "#f68b6a", "#083356", "#ce1554", "#FFFFBF", "#3288BD", "#e1abbc",
         "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#008099FF",
         "#9E0142", "#fdbf6f", "#62B197", "#E16E6D", "#9392BE", "#D0E7ED",
         "#D5E4A8", "#6a3d9a", "#4197d8", "#E79397", "#f8c120", "#CC88B0",
         "#a1d5b9", "#50b688", "#962d20", "#5b2d90", "#ffda67", "#339900FF",
         "#E7BA52FF", "#990080FF")

p <- ggplot(df, aes(x = x, fill = stratum, label = stratum,
                    stratum = stratum, alluvium = value)) +
  geom_flow(width = 0.1, curve_type = "sine", alpha = 0.6,
            color = "white", size = 0) +
  geom_stratum(width = 0.35) +
  geom_text(stat = "stratum", size = 2.5, color = "black") +
  scale_fill_manual(values = col) +
  theme_void(1) +
  theme(legend.position = "none")

ggsave(file.path(OUT_FIGURES, "24_alluvial_mediation.pdf"),
       p, width = 6.5, height = 4.5)

message("Step 24 complete.")

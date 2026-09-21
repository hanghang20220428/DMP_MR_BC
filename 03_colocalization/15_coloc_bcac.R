# ===========================================================================
#  15_coloc_bcac.R
#
#  Purpose
#  -------
#  Colocalisation of the causal CpG cis-meQTL signals with breast-cancer
#  susceptibility using the BCAC 2020 meta-analysis as the outcome.
#
#  Input
#  -----
#    output/objects/meta_analysis_significant_cpg.rdata   cpg_meta_significant
#    output/objects/All_clean_MeQTL_GoDMC.rdata           data2
#    PATHS$epic_meqtl                                     EPIC meQTL reference
#    PATHS$bcac_sus_coloc                                 GWASdata (coloc-ready)
#
#  Output
#  ------
#    output/objects/coloc_BCAC.rdata
#    output/tables/15_coloc_BCAC_summary.csv
#    output/tables/15_coloc_BCAC_per_SNP_<cpg>.csv
#
#  See the header of 14_coloc_finngen.R for the definition of the reported
#  statistic (the lead SNP and its conditional posterior, top_snp_PP.H4).
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(coloc)
  library(tidyverse)
  library(data.table)
})

source("config/config.R")
source("helpers/coloc_functions.R")

require_files(PATHS$epic_meqtl, PATHS$bcac_sus_coloc)

load(file.path(OUT_OBJECTS, "meta_analysis_significant_cpg.rdata"))  # cpg_meta_significant
load(file.path(OUT_OBJECTS, "All_clean_MeQTL_GoDMC.rdata"))          # data2
load(PATHS$epic_meqtl)      # -> EPIC_meQTL
load(PATHS$bcac_sus_coloc)  # -> GWASdata (coloc-ready BCAC table)

coloc_summary <- data.frame()
coloc_per_snp <- list()

for (cpg in cpg_meta_significant) {
  r <- run_coloc_one(cpg, data2, GWASdata, EPIC_meQTL,
                     label = "BC (BCAC 2020)")
  if (is.null(r)) { message("  ", cpg, ": not testable"); next }

  coloc_summary <- rbind(coloc_summary, r$summary)
  coloc_per_snp[[cpg]] <- r$per_snp

  write_table(r$per_snp, paste0("15_coloc_BCAC_per_SNP_", cpg, ".csv"))

  message(sprintf("  %s: lead SNP %s, top_snp_PP.H4 = %.4f",
                  cpg, r$summary$top_snp, r$summary$top_snp_PP.H4))
}

save_rdata(coloc_summary, coloc_per_snp, file = "coloc_BCAC.rdata")
write_table(coloc_summary, "15_coloc_BCAC_summary.csv")

message("Step 15 complete.")

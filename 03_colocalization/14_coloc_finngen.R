# ===========================================================================
#  14_coloc_finngen.R
#
#  Purpose
#  -------
#  Colocalisation of the causal CpG cis-meQTL signals with breast-cancer risk
#  using FinnGen R11 as the outcome.  coloc.abf() tests whether the CpG and
#  breast cancer share a single causal variant in the cis window (H4).
#
#  Input
#  -----
#    output/objects/meta_analysis_significant_cpg.rdata   cpg_meta_significant
#    output/objects/All_clean_MeQTL_GoDMC.rdata           data2
#    PATHS$epic_meqtl                                     EPIC meQTL reference
#    PATHS$finngen_bc_coloc                               finncolocdata
#
#  Output
#  ------
#    output/objects/coloc_FinnGen.rdata
#    output/tables/14_coloc_FinnGen_summary.csv     lead SNP + its posterior
#    output/tables/14_coloc_FinnGen_per_SNP_<cpg>.csv
#
#  Reported statistic
#  ------------------
#    top_snp        the lead SNP for the CpG, i.e. the SNP with the largest
#                   per-SNP posterior in the cis window.
#    top_snp_PP.H4  SNP.PP.H4 of that lead SNP - the posterior probability
#                   that this SNP is the shared causal variant, conditional
#                   on H4 being true.  This is the value reported in the
#                   manuscript for each CpG, in both outcome cohorts.
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(coloc)
  library(tidyverse)
  library(data.table)
})

source("config/config.R")
source("helpers/coloc_functions.R")

require_files(PATHS$epic_meqtl, PATHS$finngen_bc_coloc)

load(file.path(OUT_OBJECTS, "meta_analysis_significant_cpg.rdata"))  # cpg_meta_significant
load(file.path(OUT_OBJECTS, "All_clean_MeQTL_GoDMC.rdata"))          # data2
load(PATHS$epic_meqtl)        # -> EPIC_meQTL
load(PATHS$finngen_bc_coloc)  # -> finncolocdata

message("CpGs tested: ", paste(cpg_meta_significant, collapse = ", "))

coloc_summary <- data.frame()
coloc_per_snp <- list()

for (cpg in cpg_meta_significant) {
  r <- run_coloc_one(cpg, data2, finncolocdata, EPIC_meQTL,
                     label = "BC (FinnGen R11)")
  if (is.null(r)) { message("  ", cpg, ": not testable"); next }

  coloc_summary <- rbind(coloc_summary, r$summary)
  coloc_per_snp[[cpg]] <- r$per_snp

  write_table(r$per_snp, paste0("14_coloc_FinnGen_per_SNP_", cpg, ".csv"))

  message(sprintf("  %s: lead SNP %s, top_snp_PP.H4 = %.4f",
                  cpg, r$summary$top_snp, r$summary$top_snp_PP.H4))
}

save_rdata(coloc_summary, coloc_per_snp, file = "coloc_FinnGen.rdata")
write_table(coloc_summary, "14_coloc_FinnGen_summary.csv")

message("Step 14 complete.")

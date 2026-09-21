# ===========================================================================
#  07_tsmr_gcst90011804.R
#
#  Purpose
#  -------
#  Exploratory third-outcome replication of the CpG -> breast-cancer MR screen
#  using GCST90011804 (BCAC/BRIDGES/CIMBA 2020, GRCh37).
#
#  This cohort is used for the three-way intersection in step 08 and for the
#  GSMR sensitivity analysis in step 11.  It is not part of the primary
#  two-outcome analysis.
#
#  Input
#  -----
#    output/objects/GoDMC_as_exp/*.rdata     instrument sets (step 04)
#    PATHS$gcst90011804_outcome              object: outdata
#
#  Output
#  ------
#    output/objects/TSMR_GCST90011804.rdata / _nohet.rdata
#    output/tables/TSMR_GCST90011804.csv   / _nohet.csv
#
#  Difference from steps 05-06
#  ---------------------------
#  The Steiger directionality test was not applied here, so the robustness
#  filter is run with use_direction = FALSE and the result set is defined by
#  heterogeneity / pleiotropy criteria only.
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(tidyverse)
  library(data.table)
})

source("config/config.R")
source("helpers/mr_functions.R")

require_files(PATHS$gcst90011804_outcome)
load(PATHS$gcst90011804_outcome)     # -> outdata (TwoSampleMR outcome)

cpg_list <- sub("\\.rdata$", "",
                list.files(PATHS$godmc_exposure_dir, pattern = "\\.rdata$"))
message("Instrumented CpGs available: ", length(cpg_list))

# ---- 1. CpG-by-CpG MR -----------------------------------------------------
AllOR <- data.frame()

for (k in seq_along(cpg_list)) {
  cpg <- cpg_list[k]
  if (k %% 250 == 0) message("  ... ", k, " / ", length(cpg_list))

  exp <- load_godmc_instruments(cpg)
  if (is.null(exp)) next

  or <- run_tsmr_one(exp, outdata, exposure_id = cpg,
                     outcome_label = "BC_GCST90011804")
  if (is.null(or)) next

  AllOR <- rbind(AllOR, or)
}

message("CpG-outcome tests completed: ", nrow(AllOR))

# ---- 2. Multiple-testing correction --------------------------------------
AllOR <- add_multiple_testing(AllOR)

save_rdata(AllOR, file = "TSMR_GCST90011804.rdata")
write_table(AllOR, "TSMR_GCST90011804.csv")

# ---- 3. Robustness filter (no directionality test) -----------------------
AllOR_nohet <- filter_robust_mr(AllOR, use_direction = FALSE)
save_rdata(AllOR_nohet, file = "TSMR_GCST90011804_nohet.rdata")
write_table(AllOR_nohet, "TSMR_GCST90011804_nohet.csv")

message(sprintf("Nominal (P < 0.05): %d | after robustness filter: %d",
                sum(AllOR$pval < 0.05, na.rm = TRUE), nrow(AllOR_nohet)))
message("Step 07 complete.")

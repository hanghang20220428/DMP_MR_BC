# ===========================================================================
#  17_phewas_mr.R
#
#  Purpose
#  -------
#  Phenome-wide MR screen: each causal CpG (CpG set A, step 08) is tested
#  against every FinnGen release 11 endpoint, to characterise the wider
#  phenotypic consequences of the methylation instruments.
#
#  Input
#  -----
#    output/objects/GoDMC_as_exp/*.rdata              instrument sets (step 04)
#    output/objects/CpG_sets_meqtl_MR.rdata           cpg_intersect_finn_bcac
#    PATHS$finngen_all_outcomes                       <endpoint>_as_outcome.rdata
#    PATHS$finngen_manifest                           endpoint names / categories
#
#  Output
#  ------
#    output/objects/PheWAS_FinnGen.rdata
#    output/objects/PheWAS_FinnGen_nohet.rdata
#    output/objects/PheWAS_FinnGen_FDR_nohet.rdata
#    output/objects/PheWAS_FinnGen_BFR_nohet.rdata
#    output/tables/PheWAS_FinnGen.csv  and the three filtered subsets
#
#  Result sets
#  -----------
#    AllOR                all nominal results
#    *_nohet              P < 0.05 + heterogeneity / pleiotropy / directionality
#    *_FDR_nohet          same, but selected on Benjamini-Hochberg FDR
#    *_BFR_nohet          same, but selected on Bonferroni
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(tidyverse)
  library(data.table)
})

source("config/config.R")
source("helpers/mr_functions.R")

finngen_dir <- PATHS$finngen_all_outcomes
if (!dir.exists(finngen_dir)) {
  stop("FinnGen outcome directory not found: ", finngen_dir,
       "\nSet `finngen_all_outcomes` in config/paths_local.R.", call. = FALSE)
}

load(file.path(OUT_OBJECTS, "CpG_sets_meqtl_MR.rdata"))   # cpg_intersect_finn_bcac

out_files <- list.files(finngen_dir, pattern = "_as_outcome\\.rdata$")
info      <- fread(PATHS$finngen_manifest, data.table = FALSE)

message(sprintf("CpGs: %d | FinnGen endpoints: %d",
                length(cpg_intersect_finn_bcac), length(out_files)))

# ---- 1. CpG x endpoint MR -------------------------------------------------
AllOR <- data.frame()

for (i in cpg_intersect_finn_bcac) {

  exp <- load_godmc_instruments(i)
  if (is.null(exp)) next

  for (j in out_files) {
    shortj <- sub("_as.*$", "", j)

    e <- new.env(parent = emptyenv())
    load(file.path(finngen_dir, j), envir = e)
    if (!exists("finndata", envir = e)) next
    finndata <- get("finndata", envir = e)

    or <- run_tsmr_one(exp, finndata, exposure_id = i,
                       outcome_label = info$phenotype[match(shortj, info$phenocode)])
    if (is.null(or)) next
    or$id.outcome <- shortj

    AllOR <- rbind(AllOR, or)
  }
  message("  ", i, " done (", nrow(AllOR), " tests so far)")
}

message("CpG-endpoint tests completed: ", nrow(AllOR))

# ---- 2. Multiple-testing correction --------------------------------------
AllOR <- add_multiple_testing(AllOR)
save_rdata(AllOR, file = "PheWAS_FinnGen.rdata")
write_table(AllOR, "PheWAS_FinnGen.csv")

# ---- 3. Three robustness filters -----------------------------------------
AllOR_nohet     <- filter_robust_mr(AllOR, signal_col = "pval")
AllOR_FDR_nohet <- filter_robust_mr(AllOR, signal_col = "FDR")
AllOR_BFR_nohet <- filter_robust_mr(AllOR, signal_col = "BFR")

save_rdata(AllOR_nohet,     file = "PheWAS_FinnGen_nohet.rdata")
save_rdata(AllOR_FDR_nohet, file = "PheWAS_FinnGen_FDR_nohet.rdata")
save_rdata(AllOR_BFR_nohet, file = "PheWAS_FinnGen_BFR_nohet.rdata")

write_table(AllOR_nohet,     "PheWAS_FinnGen_nohet.csv")
write_table(AllOR_FDR_nohet, "PheWAS_FinnGen_FDR_nohet.csv")
write_table(AllOR_BFR_nohet, "PheWAS_FinnGen_BFR_nohet.csv")

message(sprintf("Nominal %d | FDR %d | Bonferroni %d",
                nrow(AllOR_nohet), nrow(AllOR_FDR_nohet), nrow(AllOR_BFR_nohet)))
message("Step 17 complete.")

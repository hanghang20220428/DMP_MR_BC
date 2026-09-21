# ===========================================================================
#  11_gsmr_gcst90011804.R
#
#  Purpose
#  -------
#  GSMR sensitivity analysis for CpG set C (significant in FinnGen, BCAC and
#  GCST90011804), using GCST90011804 as the outcome.
#
#  Input
#  -----
#    output/objects/CpG_sets_meqtl_MR.rdata      cpg_intersect_finn_bcac_gcst
#    output/objects/All_clean_MeQTL_GoDMC.rdata  data2
#    PATHS$gcst90011804_raw                      GCST90011804 summary statistics
#    PATHS$snp_eaf_eur                           ALFA European allele frequencies
#    PATHS$pancancer_info                        study metadata (samplesize)
#
#  Output
#  ------
#    output/objects/GSMR_GCST90011804.rdata
#    output/tables/GSMR_GCST90011804.csv
#
#  Outcome harmonisation
#  ---------------------
#  GCST90011804 is published without beta, standard error or effect-allele
#  frequency, so these are reconstructed exactly as in step 00b:
#      beta <- log(odds_ratio) ;  se <- get_se(beta, p) ; eaf <- ALFA EUR
#
#  !! See docs/04_known_issues.md for the GSMR argument-order issue !!
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(gsmr2)
  library(tidyverse)
  library(data.table)
})

source("config/config.R")
source("helpers/gsmr_functions.R")

require_files(PATHS$gcst90011804_raw, PATHS$snp_eaf_eur, PATHS$pancancer_info)

# ---- 1. CpG set C ---------------------------------------------------------
load(file.path(OUT_OBJECTS, "CpG_sets_meqtl_MR.rdata"))   # cpg_intersect_finn_bcac_gcst
load(file.path(OUT_OBJECTS, "All_clean_MeQTL_GoDMC.rdata"))
load(PATHS$snp_eaf_eur)        # -> snp_eaf_eur (SNP, eaf_eur)
load(PATHS$pancancer_info)     # -> info (accessionId, samplesize)

message("CpGs in set C (FinnGen + BCAC + GCST90011804): ",
        length(cpg_intersect_finn_bcac_gcst))

# ---- 2. Outcome table -----------------------------------------------------
data0 <- fread(PATHS$gcst90011804_raw, data.table = FALSE)

data0$beta           <- log(data0$odds_ratio)
data0$standard_error <- TwoSampleMR::get_se(data0$beta, data0$p_value)
data0$effect_allele_frequency <-
  snp_eaf_eur$eaf_eur[data.table::chmatch(data0$variant_id, snp_eaf_eur$SNP)]
data0$samplesize <- info$samplesize[match("GCST90011804", info$accessionId)]

save_rdata(data0, file = "GCST90011804_data0.rdata")

outdata <- data0 %>%
  dplyr::select(SNP           = "variant_id",
                effect_allele = "effect_allele",
                other_allele  = "other_allele",
                beta          = "beta",
                se            = "standard_error",
                eaf           = "effect_allele_frequency",
                p             = "p_value",
                chr           = "chromosome",
                pos           = "base_pair_location",
                samplesize    = "samplesize") %>%
  na.omit()

outdata <- gsmr_outcome_table(outdata)

message(sprintf("GCST90011804 outcome: %s SNPs, N = %s",
                format(nrow(outdata), big.mark = ","),
                format(unique(outdata$samplesize), big.mark = ",")))

# ---- 3. GSMR scan ---------------------------------------------------------
AllOR <- run_gsmr_scan(cpg_intersect_finn_bcac_gcst, data2, outdata,
                       outcome_label = "BC_GCST90011804")

save_rdata(AllOR, file = "GSMR_GCST90011804.rdata")
write_table(AllOR, "GSMR_GCST90011804.csv")

message(sprintf("GSMR associations at P < 0.05: %d of %d CpGs tested",
                sum(AllOR$pval < 0.05, na.rm = TRUE), nrow(AllOR)))
message("Step 11 complete.")

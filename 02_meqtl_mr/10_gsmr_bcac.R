# ===========================================================================
#  10_gsmr_bcac.R
#
#  Purpose
#  -------
#  GSMR sensitivity analysis for CpG set A (FinnGen + BCAC intersection),
#  using the BCAC 2020 susceptibility meta-analysis as the outcome.
#
#  Input
#  -----
#    output/objects/CpG_sets_meqtl_MR.rdata      cpg_intersect_finn_bcac
#    output/objects/All_clean_MeQTL_GoDMC.rdata  data2
#    PATHS$bcac_sus_raw                          ICOS/OncoArray summary statistics
#
#  Output
#  ------
#    output/objects/GSMR_BCAC.rdata
#    output/tables/GSMR_BCAC.csv
#
#  Outcome harmonisation
#  ---------------------
#  The BCAC file stores SNP identifiers as "rsID:allele1:allele2", so the
#  first field is kept; SNP rows that do not carry an rsID are dropped.
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

require_files(PATHS$bcac_sus_raw)

# ---- 1. CpG set A ---------------------------------------------------------
load(file.path(OUT_OBJECTS, "CpG_sets_meqtl_MR.rdata"))   # cpg_intersect_finn_bcac
load(file.path(OUT_OBJECTS, "All_clean_MeQTL_GoDMC.rdata"))

message("CpGs in set A (FinnGen + BCAC intersection): ",
        length(cpg_intersect_finn_bcac))

# ---- 2. Outcome table -----------------------------------------------------
BCAC_sus <- fread(PATHS$bcac_sus_raw, data.table = FALSE)

outdata <- BCAC_sus %>%
  dplyr::select("SNP.Onco", "Effect.Onco", "Baseline.Onco", "Beta.meta",
                "sdE.meta", "EAFcases.Onco", "p.meta", "chr.Onco", "Position.Onco")

colnames(outdata) <- c("SNP", "effect_allele", "other_allele", "beta", "se",
                       "eaf", "p", "chr", "pos")

outdata$SNP <- sapply(strsplit(outdata$SNP, split = ":"), "[", 1)
outdata     <- outdata[grepl("rs", outdata$SNP), ]

outdata$samplesize <- 247173
outdata <- gsmr_outcome_table(outdata)

message(sprintf("BCAC outcome: %s SNPs, N = %s",
                format(nrow(outdata), big.mark = ","),
                format(unique(outdata$samplesize), big.mark = ",")))

# ---- 3. GSMR scan ---------------------------------------------------------
AllOR <- run_gsmr_scan(cpg_intersect_finn_bcac, data2, outdata,
                       outcome_label = "BC_BCAC")

save_rdata(AllOR, file = "GSMR_BCAC.rdata")
write_table(AllOR, "GSMR_BCAC.csv")

message(sprintf("GSMR associations at P < 0.05: %d of %d CpGs tested",
                sum(AllOR$pval < 0.05, na.rm = TRUE), nrow(AllOR)))
message("Step 10 complete.")

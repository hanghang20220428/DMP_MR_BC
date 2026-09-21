# ===========================================================================
#  09_gsmr_finngen.R
#
#  Purpose
#  -------
#  GSMR sensitivity analysis for the CpGs that were significant in both the
#  FinnGen and BCAC TSMR screens (CpG set A, step 08), using FinnGen R11
#  breast cancer as the outcome.
#
#  GSMR is used as an independent causal-inference implementation that adds
#  HEIDI-outlier filtering, i.e. it tests whether the exposure and outcome
#  associations are driven by the same underlying variant rather than by
#  linkage disequilibrium.
#
#  Input
#  -----
#    output/objects/CpG_sets_meqtl_MR.rdata   cpg_intersect_finn_bcac
#    output/objects/All_clean_MeQTL_GoDMC.rdata   data2
#    PATHS$finngen_bc_raw, PATHS$finngen_manifest
#
#  Output
#  ------
#    output/objects/GSMR_FinnGen.rdata
#    output/tables/GSMR_FinnGen.csv
#
#  !! See docs/04_known_issues.md !!
#  The published GSMR runs used a positionally mis-aligned argument list
#  against gsmr2 1.1.1, which effectively disabled HEIDI-outlier filtering and
#  the LD-FDR pruning step.  GSMR_PARAMS$legacy_published_arguments = TRUE
#  reproduces those runs; set it to FALSE to use the intended parameters.
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(gsmr2)
  library(tidyverse)
  library(data.table)
})

source("config/config.R")
source("helpers/gsmr_functions.R")

require_files(PATHS$finngen_bc_raw, PATHS$finngen_manifest)

# ---- 1. CpG set A ---------------------------------------------------------
load(file.path(OUT_OBJECTS, "CpG_sets_meqtl_MR.rdata"))   # cpg_intersect_finn_bcac
load(file.path(OUT_OBJECTS, "All_clean_MeQTL_GoDMC.rdata"))

message("CpGs in set A (FinnGen + BCAC intersection): ",
        length(cpg_intersect_finn_bcac))

# ---- 2. Outcome table -----------------------------------------------------
data0 <- fread(PATHS$finngen_bc_raw, data.table = FALSE)

finn_info <- fread(PATHS$finngen_manifest, data.table = FALSE)
trait_row <- finn_info[grepl("finngen_R11_C3_BREAST_EXALLC.gz",
                             finn_info$path_https), ]

outdata <- data0 %>%
  dplyr::select(SNP           = "rsids",
                effect_allele = "alt",
                other_allele  = "ref",
                beta          = "beta",
                se            = "sebeta",
                eaf           = "af_alt",
                p             = "pval",
                chr           = "#chrom",
                pos           = "pos") %>%
  na.omit() %>%
  separate_rows(SNP, sep = ",")          # FinnGen lists multiple rsIDs per row

outdata$samplesize <- trait_row$num_cases + trait_row$num_controls
outdata <- gsmr_outcome_table(outdata)

message(sprintf("FinnGen outcome: %s SNPs, N = %s",
                format(nrow(outdata), big.mark = ","),
                format(unique(outdata$samplesize), big.mark = ",")))

# ---- 3. GSMR scan ---------------------------------------------------------
AllOR <- run_gsmr_scan(cpg_intersect_finn_bcac, data2, outdata,
                       outcome_label = "BC_FinnGenR11")

save_rdata(AllOR, file = "GSMR_FinnGen.rdata")
write_table(AllOR, "GSMR_FinnGen.csv")

message(sprintf("GSMR associations at P < 0.05: %d of %d CpGs tested",
                sum(AllOR$pval < 0.05, na.rm = TRUE), nrow(AllOR)))
message("Step 09 complete.")

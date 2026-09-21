# ===========================================================================
#  00a_finngen_r11_bc_to_coloc.R
#
#  Purpose
#  -------
#  Convert the FinnGen release 11 breast-cancer summary statistics
#  (endpoint C3_BREAST_EXALLC) into the "coloc-ready" table used by the
#  colocalisation step (step 13).  The conversion adds the case/control
#  counts required by coloc.abf() and derives varbeta, MAF, s and z.
#
#  Input
#  -----
#    PATHS$finngen_bc_raw    finngen_R11_C3_BREAST_EXALLC.gz   (FinnGen R11)
#    PATHS$finngen_manifest  finngen_R11_manifest.csv          (endpoint metadata)
#
#  Output
#  ------
#    OUTPUT_DIR/objects/finndataR11_BC_as_coloc.rdata   object: finncolocdata
#      SNP  chr  pos  effect_allele  other_allele  eaf  beta  se  P
#      samplesize  ncase  varbeta  MAF  s  z
#
#  Notes
#  -----
#    * FinnGen R11 breast cancer: 20,586 cases and 201,494 controls
#      (N = 222,080), GRCh38.
#    * `s` is the proportion of cases, as required by coloc for type = "cc".
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

source("config/config.R")
require_files(PATHS$finngen_bc_raw, PATHS$finngen_manifest)

# ---- 1. Read the raw FinnGen endpoint and its manifest entry --------------
data0     <- fread(PATHS$finngen_bc_raw, data.table = FALSE)
finn_info <- fread(PATHS$finngen_manifest, data.table = FALSE)

# Match the manifest row describing this endpoint so that the exact
# case/control counts released by FinnGen are used.
trait_row <- finn_info[grepl("finngen_R11_C3_BREAST_EXALLC.gz", finn_info$path_https), ]
stopifnot(nrow(trait_row) == 1)

data0$ncase      <- trait_row$num_cases
data0$ncontrol   <- trait_row$num_controls
data0$samplesize <- trait_row$num_cases + trait_row$num_controls

message(sprintf("FinnGen R11 C3_BREAST_EXALLC: %s cases / %s controls (N = %s)",
                format(trait_row$num_cases, big.mark = ","),
                format(trait_row$num_controls, big.mark = ","),
                format(trait_row$samplesize, big.mark = ",")))

# ---- 2. Harmonise column names -------------------------------------------
data1 <- data0 %>%
  dplyr::select(SNP            = "rsids",
                chr            = "#chrom",
                pos            = "pos",
                effect_allele  = "alt",
                other_allele   = "ref",
                eaf            = "af_alt",
                beta           = "beta",
                se             = "sebeta",
                P              = "pval",
                samplesize     = "samplesize",
                ncase          = "ncase")

# ---- 3. Derive the quantities coloc.abf() needs ---------------------------
data1$varbeta <- data1$se^2
data1$MAF     <- ifelse(data1$eaf < 0.5, data1$eaf, 1 - data1$eaf)
data1$z       <- data1$beta / data1$se
data1$s       <- data1$ncase / data1$samplesize          # case proportion
data1         <- data1[!duplicated(data1$SNP), ]

# Drop SNPs with missing values in any required field.  SNPs missing an EAF
# are back-filled upstream (see step 00b) when necessary.
finncolocdata <- na.omit(data1)

# ---- 4. Save --------------------------------------------------------------
save_rdata(finncolocdata, file = "finndataR11_BC_as_coloc.rdata")
message("Saved ", nrow(finncolocdata), " SNPs to output/objects/finndataR11_BC_as_coloc.rdata")

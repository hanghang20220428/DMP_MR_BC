# ===========================================================================
#  03_godmc_meqtl_cleaning.R
#
#  Purpose
#  -------
#  Clean the GoDMC whole-blood cis-meQTL meta-analysis so that it can be used
#  as the exposure resource for the cis-meQTL Mendelian randomisation screen.
#
#  Input
#  -----
#    PATHS$godmc_assoc  assoc_meta_all.csv.gz   GoDMC mQTL meta-analysis (GRCh37)
#    PATHS$godmc_snps   snps.csv.gz             SNP annotation: chr:pos <-> rsID
#
#  Output
#  ------
#    output/objects/All_clean_MeQTL_GoDMC.rdata  object: data2
#        cpg  chr_pos  effect_allele  other_allele  beta  se  eaf  p
#        samplesize  SNP  chr  pos
#    output/objects/all_cpg_list_GoDMC.rdata     object: all_GoDM_cpg_list
#    output/tables/03_GoDMC_meQTL_summary.csv
#
#  Cleaning steps
#  --------------
#    * keep only the columns required for MR
#    * map the GoDMC "chr:pos" SNP identifier to an rsID using snps.csv.gz
#    * drop pairs without an rsID or with missing statistics
#    * split chr:pos into separate chromosome and position columns
#
#  Reference build
#  ---------------
#    GoDMC is distributed on GRCh37 / hg19; all downstream tools (PLINK EUR
#    reference panel, FinnGen, BCAC) are therefore matched on GRCh37.
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

source("config/config.R")
require_files(PATHS$godmc_assoc, PATHS$godmc_snps)

# ---- 1. Read the meQTL meta-analysis -------------------------------------
GoDMC_meQTL <- fread(PATHS$godmc_assoc, data.table = FALSE)
message("GoDMC pairs read: ", format(nrow(GoDMC_meQTL), big.mark = ","))

data1 <- GoDMC_meQTL %>%
  dplyr::select("cpg", "snp", "allele1", "allele2", "beta_a1",
                "se", "freq_a1", "pval", "samplesize")

colnames(data1) <- c("cpg", "chr_pos", "effect_allele", "other_allele",
                     "beta", "se", "eaf", "p", "samplesize")

# ---- 2. Map chr:pos identifiers to rsIDs ---------------------------------
snps <- fread(PATHS$godmc_snps, data.table = FALSE)

data1$SNP <- snps$rsid[match(data1$chr_pos, snps$name)]

n_before <- nrow(data1)
data2    <- na.omit(data1)
message(sprintf("Pairs with a resolvable rsID and complete statistics: %s of %s",
                format(nrow(data2), big.mark = ","),
                format(n_before, big.mark = ",")))

# ---- 3. Split chr:pos into chromosome and position ----------------------
data2$chr <- sub("^chr", "", sapply(strsplit(data2$chr_pos, ":"), "[", 1))
data2$pos <- as.integer(sapply(strsplit(data2$chr_pos, ":"), "[", 2))

# ---- 4. Save ------------------------------------------------------------
save_rdata(data2, file = "All_clean_MeQTL_GoDMC.rdata")

all_GoDM_cpg_list <- unique(data2$cpg)
save_rdata(all_GoDM_cpg_list, file = "all_cpg_list_GoDMC.rdata")

write_table(data.frame(
  n_pairs_input         = n_before,
  n_pairs_after_clean   = nrow(data2),
  n_unique_cpg          = length(all_GoDM_cpg_list),
  median_samplesize     = median(data2$samplesize)
), "03_GoDMC_meQTL_summary.csv")

message("Unique CpGs in GoDMC: ", format(length(all_GoDM_cpg_list), big.mark = ","))

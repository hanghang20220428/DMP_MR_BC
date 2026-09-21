# ===========================================================================
#  00c_pancancer17_to_mrdata.R
#
#  Purpose
#  -------
#  Convert public GWAS Catalog summary statistics (one .tsv.gz per trait) into
#  TwoSampleMR exposure / outcome objects, so that they can be used directly
#  in the phenome-wide MR screens (steps 16 and 04).
#
#  Instrument selection (exposure side)
#  ------------------------------------
#    * P < 5 x 10^-8                     (genome-wide significant)
#    * LD clumping r^2 < 0.001, 10,000 kb window, 1000 Genomes EUR reference
#    * MHC region (chr6:28-34 Mb) removed
#
#  Input
#  -----
#    <trait_dir>/*.h.tsv.gz      GWAS Catalog harmonised summary statistics
#    <trait_dir>/info.rdata      study metadata (accessionId, trait, samplesize)
#
#  Output
#  ------
#    OUTPUT_OBJECTS/<accession>_as_outcome.rdata     object: outdata
#    OUTPUT_OBJECTS/<accession>_as_exposure.rdata    object: exp
#
#  Note
#  ----
#  Reference build: GRCh38 for the GWAS Catalog harmonised files.
#  Set `gwas_catalog_dir` in config/paths_local.R.
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(tidyverse)
  library(data.table)
})

source("config/config.R")

gwas_catalog_dir <- .resolve_path("gwas_catalog_dir",
                                  file.path(DATA_ROOT, "GwasCatalog"))
if (!dir.exists(gwas_catalog_dir)) {
  stop("GWAS Catalog directory not found: ", gwas_catalog_dir,
       "\nSet `gwas_catalog_dir` in config/paths_local.R.", call. = FALSE)
}

allf <- list.files(gwas_catalog_dir, pattern = "tsv.gz", full.names = TRUE)

# Study metadata keyed by GWAS Catalog accession ID.
load(file.path(gwas_catalog_dir, "info.rdata"))   # -> info

allMRinfo <- data.frame()

for (i in allf) {
  shorti <- gsub(".h.tsv.gz", "", basename(i))
  data0  <- fread(i, data.table = FALSE, integer64 = "numeric")

  data1 <- data0 %>%
    dplyr::select(SNP           = "rsid",
                  effect_allele = "effect_allele",
                  other_allele  = "other_allele",
                  beta          = "beta",
                  se            = "standard_error",
                  eaf           = "effect_allele_frequency",
                  p             = "p_value",
                  chr           = "chromosome",
                  pos           = "base_pair_location") %>%
    na.omit()

  data1$samplesize <- info$samplesize[match(shorti, info$accessionId)]

  # ---- (a) outcome object ------------------------------------------------
  outdata <- format_data(data1,
                         type              = "outcome",
                         snp_col           = "SNP",
                         beta_col          = "beta",
                         se_col            = "se",
                         eaf_col           = "eaf",
                         effect_allele_col = "effect_allele",
                         other_allele_col  = "other_allele",
                         pval_col          = "p",
                         chr_col           = "chr",
                         pos_col           = "pos",
                         samplesize_col    = "samplesize")
  outdata$outcome    <- info$trait[match(shorti, info$accessionId)]
  outdata$id.outcome <- shorti
  save_rdata(outdata, file = paste0(shorti, "_as_outcome.rdata"))

  # ---- (b) exposure object ------------------------------------------------
  # Remove the MHC region before instrument selection.
  mhc   <- data1 %>% dplyr::filter(chr == "6" & pos > 28000000 & pos < 34000000)
  data2 <- data1[!data1$SNP %in% mhc$SNP, ]

  exp_dat1 <- data2 %>%
    dplyr::select("SNP", "p") %>%
    dplyr::rename(rsid = SNP, pval = p) %>%
    dplyr::filter(pval < MR_PARAMS$clump_p) %>%
    dplyr::distinct(rsid, .keep_all = TRUE)

  if (nrow(exp_dat1) == 0) {
    allMRinfo <- rbind(allMRinfo, data.frame(
      trait = info$trait[match(shorti, info$accessionId)],
      traitID = shorti, error = "no SNP below 5e-8"))
    next
  }

  exp_dat2 <- try(ieugwasr::ld_clump_local(exp_dat1,
                                           clump_p  = 1,
                                           clump_r2 = MR_PARAMS$clump_r2,
                                           clump_kb = MR_PARAMS$clump_kb,
                                           bfile    = PLINK_BFILE,
                                           plink_bin = PLINK_BIN),
                  silent = TRUE)

  if (inherits(exp_dat2, "try-error")) {
    allMRinfo <- rbind(allMRinfo, data.frame(
      trait = info$trait[match(shorti, info$accessionId)],
      traitID = shorti, error = "no SNP left after LD clumping"))
    next
  }

  exp <- data2 %>%
    dplyr::filter(SNP %in% exp_dat2$rsid) %>%
    format_data(type              = "exposure",
                snp_col           = "SNP",
                beta_col          = "beta",
                se_col            = "se",
                eaf_col           = "eaf",
                effect_allele_col = "effect_allele",
                other_allele_col  = "other_allele",
                pval_col          = "p",
                chr_col           = "chr",
                pos_col           = "pos",
                samplesize_col    = "samplesize")
  exp$exposure    <- info$trait[match(shorti, info$accessionId)]
  exp$id.exposure <- shorti
  save_rdata(exp, file = paste0(shorti, "_as_exposure.rdata"))
}

if (nrow(allMRinfo) > 0) {
  write_table(allMRinfo, "pancancer_traits_without_instruments.csv")
}
message("Exposure / outcome objects written to ", OUT_OBJECTS)

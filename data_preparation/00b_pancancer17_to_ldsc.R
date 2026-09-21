# ===========================================================================
#  00b_pancancer17_to_ldsc.R
#
#  Purpose
#  -------
#  Build GWAS summary-statistic files for the 17 cancer types reported in the
#  pan-cancer analysis of Rashkin et al. (PMID 32887889).  The published
#  files (GRCh37, GRCh38 for some endpoints) do not contain beta, standard
#  error or effect-allele frequency, so these are reconstructed:
#
#     beta      <- log(odds_ratio)
#     se        <- TwoSampleMR::get_se(beta, p)
#     eaf       <- European allele frequency from the ALFA reference panel
#     N         <- cohort size parsed from the study manifest
#
#  Input
#  -----
#    PATHS$pancancer_dir/*.tsv.gz    GWAS Catalog files, e.g.
#                                    32887889-GCST90011804-EFO_0000305-Build37.f.tsv.gz
#    <study manifest>                PMID32887889_studies_export.tsv
#    PATHS$snp_eaf_eur               0snp_eaf_eur_ALFA.rdata  (object: snp_eaf_eur)
#
#  Output
#  ------
#    OUTPUT_OBJECTS/<accession>_as_LDSC.rdata     object: LDSCdata
#      SNP A1 A2 Z N
#
#  Note
#  ----
#  This step is only needed if the pan-cancer / LDSC branch of the project is
#  re-run; the breast-cancer endpoints GCST90011804 (step 06) and the FinnGen
#  phenome-wide screen (step 16) do not depend on it.
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(tidyverse)
  library(data.table)
})

source("config/config.R")

# ---- 1. Configure the pan-cancer input directory --------------------------
# Set the local path in config/paths_local.R (variable: pancancer_dir).
pancancer_dir <- .resolve_path("pancancer_dir",
                               file.path(DATA_ROOT, "PanCancer17"))
if (!dir.exists(pancancer_dir)) {
  stop("Pan-cancer directory not found: ", pancancer_dir,
       "\nSet `pancancer_dir` in config/paths_local.R.", call. = FALSE)
}

allf <- list.files(pancancer_dir, pattern = "tsv.gz", full.names = TRUE)
stopifnot(length(allf) > 0)

require_files(PATHS$snp_eaf_eur)
load(PATHS$snp_eaf_eur)   # -> snp_eaf_eur  (SNP, eaf_eur)

# ---- 2. Parse the study manifest -----------------------------------------
manifest_file <- file.path(pancancer_dir, "PMID32887889_studies_export.tsv")
require_files(manifest_file)

info <- fread(manifest_file, data.table = FALSE)

info$samplesize <- gsub(" European", "", info$discoverySampleAncestry) %>% as.numeric()
info$ncases     <- sapply(strsplit(info$initialSampleDescription,
                                   split = " European ancestry cases"), "[", 1)
info$ncontrols  <- sapply(strsplit(info$initialSampleDescription,
                                   split = " European ancestry controls"), "[", 1)
info$ncontrols  <- sapply(strsplit(info$ncontrols,
                                   split = " European ancestry cases,"), "[", 2)
info$trait      <- info$reportedTrait

save_rdata(info, file = "pancancer17_info.rdata")

# ---- 3. Convert every endpoint to LDSC format -----------------------------
for (i in allf) {
  shorti <- basename(i)
  shorti <- sapply(strsplit(shorti, split = "-EFO"), "[", 1)
  shorti <- sapply(strsplit(shorti, split = "-"),    "[", 2)

  data0 <- fread(i, data.table = FALSE, integer64 = "numeric")

  # beta and se are not published; reconstruct them from the odds ratio and
  # p-value (get_se inverts the two-sided normal test).
  data0$beta           <- log(data0$odds_ratio)
  data0$standard_error <- TwoSampleMR::get_se(data0$beta, data0$p_value)

  # Effect-allele frequency is not published either; back-fill with the
  # European ALFA reference panel.
  data0$effect_allele_frequency <-
    snp_eaf_eur$eaf_eur[data.table::chmatch(data0$variant_id, snp_eaf_eur$SNP)]

  data0$samplesize <- info$samplesize[match(shorti, info$accessionId)]

  LDSCdata <- data0 %>%
    dplyr::select("variant_id", "effect_allele", "other_allele", "beta",
                  "standard_error", "samplesize") %>%
    dplyr::mutate(Z = beta / standard_error) %>%
    dplyr::select(SNP = variant_id, A1 = effect_allele, A2 = other_allele,
                  Z = Z, N = samplesize) %>%
    na.omit()

  save_rdata(LDSCdata, file = paste0(shorti, "_as_LDSC.rdata"))
  message(sprintf("  %-14s %8d SNPs  (N = %s)", shorti, nrow(LDSCdata),
                  format(unique(LDSCdata$N)[1], big.mark = ",")))
}

message("LDSC-formatted pan-cancer files written to ", OUT_OBJECTS)

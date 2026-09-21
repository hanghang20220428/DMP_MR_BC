# ===========================================================================
#  helpers/gsmr_functions.R
#
#  Shared functions for the GSMR (Generalised Summary-data-based Mendelian
#  Randomisation) sensitivity analyses (steps 09-11).
#
#  GSMR is used here as an independent causal-inference implementation that
#  additionally applies HEIDI-outlier filtering, i.e. it tests whether the
#  exposure and outcome associations are driven by the same underlying
#  variant rather than by linkage.
#
#  All thresholds are taken from GSMR_PARAMS in config/config.R.
#
#      source("config/config.R")
#      source("helpers/gsmr_functions.R")
# ===========================================================================

#' Standardise an outcome GWAS table into the columns GSMR needs.
#'
#' @param dat  Data frame with columns SNP, effect_allele, other_allele,
#'             beta, se, eaf, p, chr, pos, samplesize.
gsmr_outcome_table <- function(dat) {
  keep <- c("SNP", "effect_allele", "other_allele", "beta", "se",
            "eaf", "p", "chr", "pos", "samplesize")
  missing <- setdiff(keep, colnames(dat))
  if (length(missing)) {
    stop("Outcome table is missing column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  dat <- dat[, keep, drop = FALSE]
  dat$chr <- as.character(dat$chr)
  dat$pos <- as.integer(dat$pos)
  na.omit(dat)
}


#' Run GSMR for a single CpG.
#'
#' @param cpg             CpG probe identifier.
#' @param meqtl           Cleaned GoDMC table (object `data2`, step 03).
#' @param outdata         Standardised outcome table (see gsmr_outcome_table).
#' @param outcome_label   Label stored in the `outcome` column of the result.
#'
#' @return One-row data frame with the GSMR estimate, or NULL when the CpG
#'         does not yield a valid GSMR run.
#'
#' Instrument handling, matching the published analysis:
#'   * exposure SNPs at P < 5e-8 for the CpG
#'   * the MHC region (chr6:28,477,797-33,448,354) is removed
#'   * the LD correlation matrix is estimated from 1000 Genomes EUR
#'   * SNPs absent from the LD matrix are dropped from both datasets
#'   * GSMR defaults from GSMR_PARAMS, including HEIDI-outlier filtering
run_gsmr_one <- function(cpg, meqtl, outdata, outcome_label) {

  # (a) exposure instruments at genome-wide significance
  data0 <- meqtl[meqtl$cpg == cpg & meqtl$p < GSMR_PARAMS$gwas_thresh, , drop = FALSE]
  if (nrow(data0) == 0) return(NULL)

  data1 <- data0[, c("SNP", "effect_allele", "other_allele", "beta", "se",
                     "eaf", "p", "chr", "pos", "samplesize"), drop = FALSE]
  data1$chr <- as.character(data1$chr)
  data1$pos <- as.integer(data1$pos)
  data1 <- na.omit(data1)
  if (nrow(data1) == 0) return(NULL)

  # (b) remove the MHC region
  mhc <- data1[data1$chr == as.character(GSMR_PARAMS$mhc_chr) &
                 data1$pos > GSMR_PARAMS$mhc_start &
                 data1$pos < GSMR_PARAMS$mhc_end, , drop = FALSE]
  exp_dat <- data1[!data1$SNP %in% mhc$SNP, , drop = FALSE]
  if (nrow(exp_dat) == 0) return(NULL)

  # (c) shared SNPs with the outcome
  same_snp <- intersect(exp_dat$SNP, outdata$SNP)
  if (length(same_snp) == 0) return(NULL)

  # (d) LD correlation matrix (1000 Genomes EUR)
  ldrho <- try(
    ieugwasr::ld_matrix_local(list(same_snp),
                              bfile       = PLINK_BFILE,
                              plink_bin   = PLINK_BIN,
                              with_alleles = FALSE),
    silent = TRUE)
  if (inherits(ldrho, "try-error") || is.null(ldrho)) return(NULL)

  ld_snp <- colnames(ldrho)
  if (length(ld_snp) < 2) return(NULL)

  # (e) align both datasets to the SNP order of the LD matrix
  exp2 <- exp_dat[match(ld_snp, exp_dat$SNP), , drop = FALSE]
  out2 <- outdata[match(ld_snp, outdata$SNP), , drop = FALSE]

  ok   <- !is.na(exp2$SNP) & !is.na(out2$SNP)
  if (sum(ok) < 2) return(NULL)
  exp2 <- exp2[ok, , drop = FALSE]
  out2 <- out2[ok, , drop = FALSE]
  ldrho <- ldrho[exp2$SNP, exp2$SNP, drop = FALSE]

  # (f) run GSMR
  #
  #  See docs/04_known_issues.md.  When legacy_published_arguments = TRUE the
  #  values are passed to the arguments they actually reached in the
  #  published run, so that the published GSMR estimates are reproduced
  #  bit-for-bit.  Set it to FALSE in config.R for the intended parameters.
  if (isTRUE(GSMR_PARAMS$legacy_published_arguments)) {
    gp <- GSMR_PARAMS$legacy_effective_values
    warning("GSMR is running in LEGACY mode: HEIDI-outlier filtering is ",
            "effectively disabled (multi_snps_heidi_thresh = 5) and the ",
            "minimum instrument count is 0.05. See docs/04_known_issues.md.",
            call. = FALSE)
  } else {
    gp <- GSMR_PARAMS
  }

  res <- try(
    gsmr(bzx            = exp2$beta,
         bzx_se         = exp2$se,
         bzx_pval       = exp2$p,
         bzy            = out2$beta,
         bzy_se         = out2$se,
         bzy_pval       = out2$p,
         ldrho          = ldrho,
         snpid          = exp2$SNP,
         n_ref          = GSMR_PARAMS$n_ref,
         heidi_outlier_flag      = GSMR_PARAMS$heidi_outlier_flag,
         gwas_thresh             = GSMR_PARAMS$gwas_thresh,
         single_snp_heidi_thresh = gp$single_snp_heidi_thresh,
         multi_snps_heidi_thresh = gp$multi_snps_heidi_thresh,
         nsnps_thresh            = gp$nsnps_thresh,
         ld_r2_thresh            = gp$ld_r2_thresh,
         ld_fdr_thresh           = gp$ld_fdr_thresh,
         gsmr2_beta              = gp$gsmr2_beta),
    silent = TRUE)
  if (inherits(res, "try-error")) return(NULL)

  out <- data.frame(exposure = cpg,
                    outcome  = outcome_label,
                    method   = "GSMR",
                    b        = res$bxy,
                    se       = res$bxy_se,
                    pval     = res$bxy_pval,
                    nsnp     = if (is.null(res$used_snp)) NA_integer_
                               else length(res$used_snp),
                    n_instruments_input = nrow(exp2))
  TwoSampleMR::generate_odds_ratios(out)
}


#' Convenience wrapper: loop GSMR over a set of CpGs.
run_gsmr_scan <- function(cpg_set, meqtl, outdata, outcome_label, verbose = TRUE) {
  AllOR <- data.frame()
  for (k in seq_along(cpg_set)) {
    cpg <- cpg_set[k]
    if (verbose && k %% 50 == 0) message("  ... ", k, " / ", length(cpg_set))
    r <- try(run_gsmr_one(cpg, meqtl, outdata, outcome_label), silent = TRUE)
    if (inherits(r, "try-error") || is.null(r)) next
    AllOR <- rbind(AllOR, r)
  }
  rownames(AllOR) <- NULL
  AllOR
}

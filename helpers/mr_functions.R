# ===========================================================================
#  helpers/mr_functions.R
#
#  Shared functions for the cis-meQTL Mendelian randomisation steps
#  (steps 04-13).  They are factored out here so that the instrument
#  selection, harmonisation and sensitivity-analysis logic is defined exactly
#  once and is applied identically to every outcome GWAS.
#
#  Source this file after config/config.R:
#
#      source("config/config.R")
#      source("helpers/mr_functions.R")
# ===========================================================================

# ---------------------------------------------------------------------------
# 1. Instrument construction
# ---------------------------------------------------------------------------

#' Build a clumped instrument set for one CpG from the GoDMC cis-meQTL table.
#'
#' @param cpg     CpG probe identifier (e.g. "cg15146122").
#' @param meqtl   Data frame of GoDMC SNP-CpG pairs, one row per pair, with
#'                columns SNP, chr, pos, cpg, effect_allele, other_allele,
#'                beta, se, eaf, p, samplesize (see step 03).
#' @param ...     Passed on to TwoSampleMR::format_data().
#'
#' @return A TwoSampleMR exposure data frame, or NULL if the CpG yields no
#'         usable instrument.
#'
#' Steps
#'   1. LD clumping with the 1000 Genomes EUR reference panel
#'      (P < 5e-8, r^2 < 0.001, 10,000 kb)
#'   2. formatting as a TwoSampleMR exposure
#'   3. weak-instrument filter: F = (beta / se)^2 > 10
build_godmc_instruments <- function(cpg, meqtl, quiet = FALSE) {

  exp_data <- meqtl[meqtl$cpg == cpg, , drop = FALSE]
  if (nrow(exp_data) == 0) return(NULL)

  exp_dat1 <- data.frame(rsid = exp_data$SNP, pval = exp_data$p)

  exp_dat2 <- try(
    ieugwasr::ld_clump_local(exp_dat1,
                             clump_p   = MR_PARAMS$clump_p,
                             clump_r2  = MR_PARAMS$clump_r2,
                             clump_kb  = MR_PARAMS$clump_kb,
                             bfile     = PLINK_BFILE,
                             plink_bin = PLINK_BIN),
    silent = TRUE)

  if (inherits(exp_dat2, "try-error") || is.null(exp_dat2) || nrow(exp_dat2) == 0) {
    return(NULL)
  }

  keep <- intersect(exp_data$SNP, exp_dat2$rsid)
  if (length(keep) == 0) return(NULL)

  data <- exp_data[exp_data$SNP %in% keep, , drop = FALSE]
  data <- data[order(data$p), , drop = FALSE]
  data <- data[!duplicated(data$SNP), , drop = FALSE]

  exp <- TwoSampleMR::format_data(
    data,
    type              = "exposure",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    eaf_col           = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    pval_col          = "p",
    samplesize_col    = "samplesize",
    chr_col           = "chr",
    pos_col           = "pos")

  # Weak-instrument filter (per-SNP F statistic)
  exp$fval <- (exp$beta.exposure / exp$se.exposure)^2
  exp <- exp[exp$fval > MR_PARAMS$f_stat_min, , drop = FALSE]

  if (nrow(exp) == 0) return(NULL)
  if (!quiet) message("    ", cpg, ": ", nrow(exp), " instruments (F > ",
                      MR_PARAMS$f_stat_min, ")")
  exp
}


#' Load a previously saved instrument set for one CpG.
load_godmc_instruments <- function(cpg, dir = PATHS$godmc_exposure_dir) {
  f <- file.path(dir, paste0(cpg, ".rdata"))
  if (!file.exists(f)) return(NULL)
  e <- new.env(parent = emptyenv())
  load(f, envir = e)
  if (!exists("exp", envir = e)) return(NULL)
  get("exp", envir = e)
}

# ---------------------------------------------------------------------------
# 2. Harmonisation and primary MR for a single CpG-outcome pair
# ---------------------------------------------------------------------------

#' Harmonise one exposure/outcome pair and run the primary MR estimate.
#'
#' @param exp            TwoSampleMR exposure data frame (instrument set).
#' @param outcome_dat    TwoSampleMR outcome data frame.
#' @param exposure_id    Label written to the `exposure` column.
#' @param outcome_label  Label written to the `outcome` column.
#' @param sensitivity_all  If FALSE (default) heterogeneity / pleiotropy /
#'                       directionality statistics are computed only for
#'                       results that reach P < 0.05, which is all the
#'                       robustness filter needs.  Set to TRUE to compute
#'                       them for every CpG, as done in the published
#'                       protein -> breast-cancer screen (step 20).
#'
#' @return One-row data frame (odds ratios plus sensitivity statistics), or
#'         NULL if the pair yields no usable instruments after harmonisation.
#'
#' Instrument filtering applied here
#'   * SNPs associated with the outcome at P < 1e-5 are removed, to limit
#'     reverse-causation / pleiotropic bias through the outcome.
#'   * palindromic SNPs are resolved with action = 2 (infer from allele
#'     frequencies); SNPs that cannot be resolved are dropped by mr_keep.
run_tsmr_one <- function(exp, outcome_dat, exposure_id, outcome_label,
                         sensitivity_all = FALSE) {

  stopifnot(inherits(exp, "data.frame"), nrow(exp) > 0)

  # (a) remove instruments that associate with the outcome
  out_snp <- outcome_dat$SNP[outcome_dat$pval.outcome < MR_PARAMS$outcome_p_remove]
  exp <- exp[!exp$SNP %in% out_snp, , drop = FALSE]
  if (nrow(exp) == 0) return(NULL)

  # (b) restrict the outcome to the instrument SNPs
  out <- outcome_dat[outcome_dat$SNP %in% exp$SNP, , drop = FALSE]
  if (nrow(out) == 0) return(NULL)

  # (c) harmonise
  dat <- try(
    TwoSampleMR::harmonise_data(exposure_dat = exp, outcome_dat = out,
                                action = MR_PARAMS$harmonise_action),
    silent = TRUE)
  if (inherits(dat, "try-error")) return(NULL)
  dat <- subset(dat, dat$mr_keep == TRUE)
  if (nrow(dat) == 0) return(NULL)

  # (d) primary estimate.  GagnonMR::primary_MR_analysis() returns the
  #     method appropriate to the number of instruments: Wald ratio for a
  #     single instrument, inverse-variance weighted (IVW) otherwise.
  res <- GagnonMR::primary_MR_analysis(dat = dat)
  res$exposure <- exposure_id
  res$outcome  <- outcome_label
  or <- TwoSampleMR::generate_odds_ratios(res)

  # (e) sensitivity analyses for the IVW estimate
  het_q  <- NA_real_
  pleio  <- NA_real_
  direct <- NA

  if (sensitivity_all || (is.finite(res$pval) && res$pval < 0.05)) {
    if (identical(as.character(or$method[1]), "Inverse variance weighted")) {
      het <- try(TwoSampleMR::mr_heterogeneity(dat, method_list = "mr_ivw"), silent = TRUE)
      if (!inherits(het, "try-error")) het_q <- het$Q_pval[1]
      pl <- try(TwoSampleMR::mr_pleiotropy_test(dat), silent = TRUE)
      if (!inherits(pl, "try-error")) pleio <- pl$pval[1]
    }
    dt <- try(TwoSampleMR::directionality_test(dat), silent = TRUE)
    if (!inherits(dt, "try-error")) direct <- dt$correct_causal_direction[1]
  }

  cbind(or,
        het_Qval    = het_q,
        pleio_pval  = pleio,
        direction   = direct,
        n_instruments = nrow(dat))
}

# ---------------------------------------------------------------------------
# 3. Multiple-testing correction and robustness filtering
# ---------------------------------------------------------------------------

add_multiple_testing <- function(AllOR) {
  AllOR$FDR <- p.adjust(AllOR$pval, method = MR_PARAMS$fdr_method)
  AllOR$BFR <- p.adjust(AllOR$pval, method = MR_PARAMS$bonferroni_method)
  AllOR
}


#' Keep only results that survive the pre-specified robustness criteria.
#'
#' A result is retained when
#'   * P < 0.05, AND
#'   * either  - Cochran's Q P > 0.05 and MR-Egger intercept P > 0.05
#'              (multi-instrument / IVW results), or
#'              - no heterogeneity test was applicable
#'              (single-instrument / Wald-ratio results), AND
#'   * the Steiger directionality test indicates the correct causal direction.
#'
#' @param use_direction  Set FALSE for the exploratory GCST90011804 screen,
#'                       where the directionality test was not applied.
#' @param signal_col     Column used as the significance criterion
#'                       ("pval", "FDR" or "BFR").
#' @param signal_cut     Cut-off applied to `signal_col` (default 0.05).
filter_robust_mr <- function(AllOR, use_direction = TRUE,
                             signal_col = "pval", signal_cut = 0.05) {

  if (!signal_col %in% colnames(AllOR)) {
    stop("Column '", signal_col, "' not found in the results table.", call. = FALSE)
  }

  sig <- !is.na(AllOR[[signal_col]]) & AllOR[[signal_col]] < signal_cut
  if (use_direction) {
    sig <- sig & !is.na(AllOR$direction) & as.logical(AllOR$direction)
  }

  pass_het <- sig &
    !is.na(AllOR$het_Qval) &
    AllOR$het_Qval  > MR_PARAMS$heterogeneity_p &
    AllOR$pleio_pval > MR_PARAMS$pleiotropy_p

  no_het <- sig & is.na(AllOR$het_Qval)

  out <- AllOR[pass_het | no_het, , drop = FALSE]
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# 4. CpG probe coordinates
# ---------------------------------------------------------------------------

#' Read an Illumina probe manifest (a GEO platform file) into a data frame.
#'
#' GEO platform files such as GPL13534 (450K) and GPL21145 (EPIC) open with a
#' block of comment lines written as `#COLUMN = ...` (for example
#' `#ID = IlmnID`) that precedes the real tab-delimited header.  These files
#' are not read with `read.delim()`; `data.table::fread()` skips the comment
#' block and lands on the true header row, so the returned frame has the
#' expected `ID`, `CHR` and `MAPINFO` columns.  The exact header layout differs
#' slightly between platform releases, so the required columns are validated
#' here and a informative error is raised up-front rather than letting a
#' missing column propagate silently as `NULL` further down the pipeline.
#'
#' @param path Path to the manifest file.
#'
#' @return A data frame with at least the columns `ID`, `CHR` and `MAPINFO`.
read_probe_manifest <- function(path) {

  require_files(path)

  df <- data.table::fread(path, data.table = FALSE, check.names = FALSE)

  needed  <- c("ID", "CHR", "MAPINFO")
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0) {
    stop("Probe manifest '", path, "' is missing column(s): ",
         paste(missing, collapse = ", "),
         ". Expected a GEO platform file (e.g. GPL13534 or GPL21145) whose ",
         "header contains ID, CHR and MAPINFO.", call. = FALSE)
  }

  df
}

#' Attach chromosome and position to a table of CpG probes.
#'
#' Coordinates come from the Illumina probe manifests (GEO platform files
#' GPL13534 for 450K and GPL21145 for EPIC, both GRCh37/hg19).  The 450K
#' manifest takes precedence because the discovery cohort (GSE104942) and the
#' GoDMC meQTL resource are both 450K-based; the EPIC manifest is used as a
#' fall-back for probes absent from the 450K array.
#'
#' @param dat        Data frame containing a CpG identifier column.
#' @param probe_col  Name of that column (default "exposure").
#'
#' @return `dat` with `CHR` and `BP` columns added.
add_probe_coordinates <- function(dat, probe_col = "exposure") {

  anno_450k <- read_probe_manifest(PATHS$anno_450k)
  anno_epic <- read_probe_manifest(PATHS$anno_epic)

  p <- dat[[probe_col]]

  chr <- anno_epic$CHR[match(p, anno_epic$ID)]
  bp  <- anno_epic$MAPINFO[match(p, anno_epic$ID)]

  chr450 <- anno_450k$CHR[match(p, anno_450k$ID)]
  bp450  <- anno_450k$MAPINFO[match(p, anno_450k$ID)]

  dat$CHR <- ifelse(is.na(chr450), chr, chr450)
  dat$BP  <- ifelse(is.na(bp450),  bp,  bp450)
  dat
}

# ===========================================================================
#  helpers/coloc_functions.R
#
#  Shared functions for the colocalisation analysis (steps 14-16).
#
#  For each CpG, a +/- 1 Mb cis window is extracted from the GoDMC meQTL
#  summary statistics and tested against the breast-cancer GWAS with
#  coloc.abf(), conditioning on a single shared causal variant (H4).
#
#  Reported statistic
#  ------------------
#  The per-SNP posterior table returned by coloc.abf() is ordered by
#  SNP.PP.H4, and for each CpG the reported value is taken from its first row:
#
#    top_snp        the lead SNP, i.e. the SNP with the largest per-SNP
#                   posterior in the cis window.
#    top_snp_PP.H4  SNP.PP.H4 of that lead SNP - the posterior probability
#                   that this SNP is the shared causal variant, conditional
#                   on H4 being true.
#
#  Every number in the summary table can therefore be traced back to a single
#  SNP, and the full per-SNP table is written out alongside it.
#
#      source("config/config.R")
#      source("helpers/coloc_functions.R")
# ===========================================================================

#' Extract the CpG coordinate from the EPIC meQTL reference table.
cpg_position <- function(cpg, epic_meqtl) {
  a <- epic_meqtl[epic_meqtl$CpG == cpg, c("CpG", "chr_cpg", "pos_cpg")]
  if (nrow(a) == 0) return(NULL)
  list(chr = as.character(a$chr_cpg[1]), pos = as.integer(a$pos_cpg[1]))
}


#' Run coloc.abf for one CpG against one outcome GWAS.
#'
#' @param cpg         CpG probe identifier.
#' @param meqtl       Cleaned GoDMC table (object `data2`, step 03).
#' @param gwas        coloc-ready outcome table with columns
#'                    SNP, beta, se, varbeta, samplesize, MAF, z, P (and s for cc);
#'                    `data_preparation/00a_*.R` builds this table.
#' @param epic_meqtl  EPIC meQTL reference table giving CpG coordinates.
#' @param label       Outcome label.
#'
#' @return A list with `summary` (one row: the lead SNP and its conditional
#'         posterior probability), `label` and `per_snp` (the full per-SNP
#'         table), or NULL if the region cannot be tested.
run_coloc_one <- function(cpg, meqtl, gwas, epic_meqtl, label) {

  pos <- cpg_position(cpg, epic_meqtl)
  if (is.null(pos)) return(NULL)

  # ---- QTL side -----------------------------------------------------------
  d <- meqtl[meqtl$cpg == cpg, , drop = FALSE]
  if (nrow(d) == 0) return(NULL)

  d$MAF     <- ifelse(d$eaf < 0.5, d$eaf, 1 - d$eaf)
  d$varbeta <- d$se^2
  d$z       <- d$beta / d$se

  window <- COLOC_PARAMS$window_kb * 1000
  qtl <- d[as.character(d$chr) == pos$chr &
             d$pos > pos$pos - window & d$pos < pos$pos + window, , drop = FALSE]
  qtl <- na.omit(qtl)
  qtl <- qtl[!duplicated(qtl$SNP), ]

  # ---- shared SNPs --------------------------------------------------------
  same_snp <- intersect(qtl$SNP, gwas$SNP)
  if (length(same_snp) == 0) return(NULL)

  qtl <- qtl[qtl$SNP %in% same_snp, , drop = FALSE]
  qtl <- qtl[order(qtl$SNP), , drop = FALSE]
  qtl <- qtl[!duplicated(qtl$SNP), , drop = FALSE]

  gw <- gwas[gwas$SNP %in% same_snp, , drop = FALSE]
  gw <- gw[order(gw$SNP), , drop = FALSE]
  gw <- gw[!duplicated(gw$SNP), , drop = FALSE]

  # ---- coloc.abf ----------------------------------------------------------
  # Two API details of coloc >= 5, both of which fail *silently* in a pipeline:
  #  * there is no `silent` argument (the published script never passed one);
  #  * dataset2$s must be a single number - check_dataset() tests
  #    `!is.numeric(d$s) || d$s <= 0 || d$s >= 1`, and a length > 1 vector only
  #    survived on R 4.2.3, where `||` used the first element with a warning.
  #    R >= 4.3 errors, so the case proportion is passed as a scalar. It is
  #    constant within an outcome study, and the two forms were verified to give
  #    bit-identical posteriors (max |diff| = 0 over the whole per-SNP vector).
  s_cc <- unique(gw$s)
  if (length(s_cc) != 1) s_cc <- s_cc[1]

  fit <- try(coloc::coloc.abf(
    dataset1 = list(snp = qtl$SNP, beta = qtl$beta, varbeta = qtl$varbeta,
                    N = qtl$samplesize, MAF = qtl$MAF, z = qtl$z,
                    pvalues = qtl$p, type = COLOC_PARAMS$type1),
    dataset2 = list(snp = gw$SNP, beta = gw$beta, varbeta = gw$varbeta,
                    N = gw$samplesize, MAF = gw$MAF, z = gw$z,
                    pvalues = gw$P, type = COLOC_PARAMS$type2, s = s_cc)))

  if (inherits(fit, "try-error")) return(NULL)

  per_snp <- fit$results
  per_snp <- per_snp[order(-per_snp$SNP.PP.H4), , drop = FALSE]
  top     <- per_snp[1, , drop = FALSE]

  list(
    label   = label,
    summary = data.frame(
      cpg           = cpg,
      outcome       = label,
      n_snp         = length(same_snp),
      top_snp       = top$snp,
      top_snp_PP.H4 = top$SNP.PP.H4),
    per_snp = per_snp)
}

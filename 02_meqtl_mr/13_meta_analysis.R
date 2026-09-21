# ===========================================================================
#  13_meta_analysis.R
#
#  Purpose
#  -------
#  Random-effects meta-analysis of the CpG -> breast-cancer MR estimates
#  across the four method x outcome combinations (TSMR/GSMR x FinnGen/BCAC),
#  and selection of the CpGs that are significant in the pooled estimate.
#
#  The meta-analysis is what defines the final causal CpG set that enters the
#  colocalisation (steps 14-16) and protein mediation (steps 19-24) analyses.
#
#  Input
#  -----
#    output/objects/12_MR_estimates_all_methods.csv   (written by step 12)
#
#  Output
#  ------
#    output/objects/meta_analysis_significant_cpg.rdata   cpg_meta_significant
#    output/tables/13_meta_analysis_results.csv
#    output/figures/13_meta_<CpG>.pdf                     one forest plot per CpG
#
#  Decision rule (as pre-specified in the manuscript)
#  --------------------------------------------------
#    If the estimates are homogeneous (Cochran's Q P > 0.05 and I^2 < 50%),
#    the common-effect estimate decides significance; otherwise the
#    random-effects estimate is used.  A CpG is called significant when the
#    relevant pooled P value is < 0.05.
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(meta)
  library(tidyverse)
  library(data.table)
})

source("config/config.R")

allmr_file <- file.path(OUT_TABLES, "12_MR_estimates_all_methods.csv")
require_files(allmr_file)

allmr <- fread(allmr_file, data.table = FALSE)

# One study label per method x outcome combination.
allmr$study <- paste0(allmr$outcome, allmr$method)

same_cpg <- unique(allmr$exposure)
message("CpGs entering the meta-analysis: ", length(same_cpg))

settings.meta("JAMA")

meta_rows <- list()
cpg_meta_significant <- character(0)

for (cpg in same_cpg) {

  cpgmr <- allmr[allmr$exposure == cpg, , drop = FALSE]

  fit <- meta::metagen(TE = b, seTE = se, data = cpgmr, studlab = study,
                       sm = "OR", common = TRUE, random = TRUE)

  # Forest plot (columns of the plot are unchanged from the analysis script)
  pdf(file.path(OUT_FIGURES, paste0("13_meta_", cpg, ".pdf")),
      width = 8, height = 6)
  meta::forest(fit)
  dev.off()

  homogeneous <- isTRUE(fit$pval.Q > 0.05 & fit$I2 < 0.5)
  p_used      <- if (homogeneous) fit$pval.common else fit$pval.random
  significant <- isTRUE(p_used < 0.05)
  if (significant) cpg_meta_significant <- c(cpg_meta_significant, cpg)

  meta_rows[[cpg]] <- data.frame(
    cpg            = cpg,
    n_studies      = fit$k,
    b_common       = fit$TE.common,
    se_common      = fit$seTE.common,
    p_common       = fit$pval.common,
    b_random       = fit$TE.random,
    se_random      = fit$seTE.random,
    p_random       = fit$pval.random,
    OR_random      = exp(fit$TE.random),
    OR_lci95       = exp(fit$lower.random),
    OR_uci95       = exp(fit$upper.random),
    Q              = fit$Q,
    p_Q            = fit$pval.Q,
    I2_percent     = 100 * fit$I2,
    model_used     = if (homogeneous) "common" else "random",
    p_used         = p_used,
    significant    = significant)
}

meta_results <- do.call(rbind, meta_rows)
write_table(meta_results, "13_meta_analysis_results.csv")

save_rdata(cpg_meta_significant, meta_results,
           file = "meta_analysis_significant_cpg.rdata")

message("Significant CpGs after meta-analysis: ",
        paste(cpg_meta_significant, collapse = ", "))
message("Step 13 complete.")

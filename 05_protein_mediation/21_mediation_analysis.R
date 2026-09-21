# ===========================================================================
#  21_mediation_analysis.R
#
#  Purpose
#  -------
#  Product-of-coefficients mediation analysis for the CpG -> protein -> breast
#  cancer pathway.
#
#      a  = CpG -> protein effect        (step 19)
#      b  = protein -> breast cancer     (step 20)
#      a*b = indirect (mediated) effect
#
#  Only proteins that are significant in BOTH arms are eligible mediators, so
#  the CpG -> protein and protein -> breast cancer result sets are
#  intersected before the indirect effect is estimated.
#
#  Input
#  -----
#    output/objects/CpG_protein_TSMR_nohet.rdata    CpG -> protein
#    output/objects/protein_BC_TSMR_nohet.rdata     protein -> breast cancer
#
#  Output
#  ------
#    output/objects/mediation_results.rdata
#        all_cpg_protein_mediation   all tested CpG-protein pairs
#        sig_cpg_protein_mediation   pairs with a 95% CI excluding zero
#    output/tables/21_mediation_<cpg>.csv          per-CpG results
#    output/tables/21_mediation_all.csv
#    output/tables/21_mediation_significant.csv
#
#  Confidence interval
#  -------------------
#  RMediation::medci with type = "prodclin" (distribution of the product),
#  which is appropriate for the product of two regression coefficients and
#  does not assume normality of a*b.
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(RMediation)
})

source("config/config.R")

load(file.path(OUT_OBJECTS, "CpG_protein_TSMR_nohet.rdata"))    # AllOR_nohet
cpg_protein <- AllOR_nohet

load(file.path(OUT_OBJECTS, "protein_BC_TSMR_nohet.rdata"))     # AllOR_nohet
protein_bc <- AllOR_nohet

all_mediation <- data.frame()

for (cpg in unique(cpg_protein$exposure)) {

  arm1 <- cpg_protein[cpg_protein$exposure == cpg, , drop = FALSE]

  # proteins significant in both arms
  shared_protein <- intersect(arm1$outcome, protein_bc$exposure)
  if (length(shared_protein) == 0) {
    message("  ", cpg, ": no protein significant in both arms")
    next
  }

  arm1_shared <- arm1[arm1$outcome %in% shared_protein, , drop = FALSE]
  arm2_shared <- protein_bc[protein_bc$exposure %in% shared_protein, , drop = FALSE]

  rows <- vector("list", length(shared_protein))

  for (k in seq_along(shared_protein)) {

    mediator <- shared_protein[k]
    d1 <- arm1_shared[arm1_shared$outcome == mediator, , drop = FALSE]
    d2 <- arm2_shared[arm2_shared$exposure == mediator, , drop = FALSE]
    if (nrow(d1) == 0 || nrow(d2) == 0) next

    beta1 <- d1$b[1]; se_beta1 <- d1$se[1]     # CpG -> protein
    beta2 <- d2$b[1]; se_beta2 <- d2$se[1]     # protein -> BC

    ci <- medci(mu.x = beta1, mu.y = beta2,
                se.x = se_beta1, se.y = se_beta2,
                alpha = MEDIATION_PARAMS$ci_alpha,
                type  = MEDIATION_PARAMS$ci_type)

    rows[[k]] <- data.frame(
      SNP           = d1$nsnp[1],
      exposure      = cpg,
      mediator      = mediator,
      outcome       = "BC",
      beta1         = beta1,
      beta2         = beta2,
      Mbeta         = beta1 * beta2,
      lo_ci         = ci$`95% CI`[1],
      up_ci         = ci$`95% CI`[2])
  }

  allME <- do.call(rbind, rows)
  if (is.null(allME) || nrow(allME) == 0) next

  save_rdata(allME, file = paste0("mediation_", cpg, ".rdata"))
  write_table(allME, paste0("21_mediation_", cpg, ".csv"))

  all_mediation <- rbind(all_mediation, allME)
  message(sprintf("  %s: %d candidate mediators", cpg, nrow(allME)))
}

# A mediator is significant when the 95% CI of the indirect effect excludes
# zero, i.e. the lower and upper limits share the same sign.
sig_mediation <- all_mediation[all_mediation$lo_ci * all_mediation$up_ci > 0, , drop = FALSE]
rownames(sig_mediation) <- NULL

save_rdata(all_mediation, sig_mediation, file = "mediation_results.rdata")
write_table(all_mediation, "21_mediation_all.csv")
write_table(sig_mediation, "21_mediation_significant.csv")

message(sprintf("Significant mediated CpG-protein pairs: %d", nrow(sig_mediation)))
if (nrow(sig_mediation)) print(sig_mediation[, c("exposure", "mediator", "Mbeta")])
message("Step 21 complete.")

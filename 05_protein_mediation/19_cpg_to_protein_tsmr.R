# ===========================================================================
#  19_cpg_to_protein_tsmr.R
#
#  Purpose
#  -------
#  First arm of the mediation analysis: test whether the causal CpGs are
#  associated with the abundance of plasma proteins (UK Biobank Pharma
#  Proteomics Project, 2,940 SomaScan aptamers).
#
#      exposure : GoDMC cis-meQTL instruments for each causal CpG
#      outcome  : UKB-PPP plasma protein levels (one file per protein)
#
#  Input
#  -----
#    output/objects/GoDMC_as_exp/*.rdata                instrument sets (step 04)
#    output/objects/meta_analysis_significant_cpg.rdata cpg_meta_significant
#    PATHS$ukbppp_outcome_dir                           <protein>_as_outcome.rdata
#    PATHS$generef_grch37                               gene coordinates (GRCh37)
#
#  Output
#  ------
#    output/objects/CpG_protein_TSMR.rdata / _nohet.rdata
#    output/tables/CpG_protein_TSMR.csv   / _nohet.csv
#    output/figures/19_manhattan_<cpg>.pdf             one plot per causal CpG
#    output/figures/19_enrichment_<cpg>.pdf
#
#  Robustness filter
#  -----------------
#    P < 0.05 AND (Cochran's Q P > 0.05 AND MR-Egger intercept P > 0.05, or
#    no heterogeneity test applicable).  The Steiger directionality test was
#    not applied in this arm, so use_direction = FALSE.
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(tidyverse)
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggprism)
  library(gground)
})

source("config/config.R")
source("helpers/mr_functions.R")
source("helpers/enrichment_plot.R")

prot_dir <- PATHS$ukbppp_outcome_dir
if (!dir.exists(prot_dir)) {
  stop("UKB-PPP outcome directory not found: ", prot_dir,
       "\nSet `ukbppp_outcome_dir` in config/paths_local.R.", call. = FALSE)
}

load(file.path(OUT_OBJECTS, "meta_analysis_significant_cpg.rdata"))  # cpg_meta_significant

prot_files <- list.files(prot_dir, pattern = "\\.rdata$")
message(sprintf("Causal CpGs: %d | UKB-PPP proteins: %d",
                length(cpg_meta_significant), length(prot_files)))

# ---- 1. CpG x protein MR --------------------------------------------------
AllOR <- data.frame()

for (i in cpg_meta_significant) {

  exp <- load_godmc_instruments(i)
  if (is.null(exp)) next

  for (j in prot_files) {

    shortj <- sapply(strsplit(j, "_"), "[", 1)      # protein name

    e <- new.env(parent = emptyenv())
    load(file.path(prot_dir, j), envir = e)
    if (!exists("outdata", envir = e)) next

    or <- run_tsmr_one(exp, get("outdata", envir = e),
                       exposure_id = i, outcome_label = shortj)
    if (is.null(or)) next

    AllOR <- rbind(AllOR, or)
  }
  message("  ", i, " done (", nrow(AllOR), " tests so far)")
}

message("CpG-protein tests completed: ", nrow(AllOR))

# ---- 2. Multiple-testing correction --------------------------------------
AllOR <- add_multiple_testing(AllOR)
save_rdata(AllOR, file = "CpG_protein_TSMR.rdata")
write_table(AllOR, "CpG_protein_TSMR.csv")

AllOR_nohet <- filter_robust_mr(AllOR, use_direction = FALSE)
save_rdata(AllOR_nohet, file = "CpG_protein_TSMR_nohet.rdata")
write_table(AllOR_nohet, "CpG_protein_TSMR_nohet.csv")

message(sprintf("CpG-protein associations at P < 0.05 after robustness filter: %d",
                nrow(AllOR_nohet)))

# ---- 3. Manhattan plots per CpG ------------------------------------------
load(PATHS$generef_grch37)   # -> generef (gene, chromosome, start)

AllOR_nohet$CHR <- generef$chromosome[match(AllOR_nohet$outcome, generef$gene)]
AllOR_nohet$BP  <- generef$start[match(AllOR_nohet$outcome, generef$gene)]

for (i in unique(AllOR_nohet$exposure)) {

  dat <- AllOR_nohet %>%
    dplyr::filter(exposure == i, !is.na(CHR), !is.na(BP)) %>%
    dplyr::select(Protein = outcome, CHR, BP, P = pval)

  if (nrow(dat) == 0) next

  prosig <- unique(dat$Protein[dat$P < 0.005])

  pdf(file.path(OUT_FIGURES, paste0("19_manhattan_", i, ".pdf")), height = 5, width = 8)
  CMplot::CMplot(dat, plot.type = "m", threshold = 0.05, cex = 0.5,
                 amplify = TRUE,
                 highlight = prosig, highlight.text = prosig,
                 highlight.text.col = "black",
                 chr.den.col = c("#1F78B4", "#f8c120"),
                 file.output = FALSE,
                 main = paste0("Proteins causally associated with ", i))
  dev.off()
}

# ---- 4. Enrichment per CpG -----------------------------------------------
for (i in unique(AllOR_nohet$exposure)) {

  genes <- AllOR_nohet$outcome[AllOR_nohet$exposure == i]
  pathway <- enrich_go_kegg(genes)
  plot_enrichment(pathway,
                  file.path(OUT_FIGURES, paste0("19_enrichment_", i, ".pdf")),
                  width = 12, height = 11)
}

message("Step 19 complete.")

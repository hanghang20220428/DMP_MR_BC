# ===========================================================================
#  20_protein_to_bc_tsmr.R
#
#  Purpose
#  -------
#  Second arm of the mediation analysis: test whether plasma protein
#  abundance is causally associated with breast-cancer risk.
#
#      exposure : UKB-PPP cis-pQTL instruments (up to 6 SNPs per protein)
#      outcome  : BCAC 2020 breast-cancer susceptibility
#
#  Together with step 19 (CpG -> protein) this provides the two coefficients
#  combined in the mediation analysis of step 21.
#
#  Input
#  -----
#    PATHS$ukbppp_exposure_dir   <protein>_as_exposure.rdata
#    PATHS$bcac_sus_meta         object: BCAC_sus_2020
#    PATHS$generef_grch38        gene coordinates for the Manhattan plot
#
#  Output
#  ------
#    output/objects/protein_BC_TSMR.rdata / _nohet.rdata
#    output/tables/protein_BC_TSMR.csv   / _nohet.csv
#    output/figures/20_manhattan_protein_BC.pdf
#    output/figures/20_enrichment_protein_BC.pdf
#
#  Robustness filter
#  -----------------
#    P < 0.05 AND (Cochran's Q P > 0.05 AND MR-Egger intercept P > 0.05, or
#    a single-instrument Wald ratio) AND correct Steiger direction.
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

exp_dir <- PATHS$ukbppp_exposure_dir
if (!dir.exists(exp_dir)) {
  stop("UKB-PPP exposure directory not found: ", exp_dir,
       "\nSet `ukbppp_exposure_dir` in config/paths_local.R.", call. = FALSE)
}

require_files(PATHS$bcac_sus_meta)
load(PATHS$bcac_sus_meta)          # -> BCAC_sus_2020

prot_files <- list.files(exp_dir, pattern = "\\.rdata$")
message("UKB-PPP proteins with instruments: ", length(prot_files))

# ---- 1. Protein x breast cancer MR ---------------------------------------
AllOR <- data.frame()

for (k in seq_along(prot_files)) {

  if (k %% 250 == 0) message("  ... ", k, " / ", length(prot_files))

  e <- new.env(parent = emptyenv())
  load(file.path(exp_dir, prot_files[k]), envir = e)
  if (!exists("exp", envir = e)) next

  exp    <- get("exp", envir = e)
  shorti <- sapply(strsplit(prot_files[k], "_"), "[", 1)

  # Weak-instrument filter (per-SNP F statistic), as in step 04
  exp$fval <- (exp$beta.exposure / exp$se.exposure)^2
  exp <- exp[exp$fval > MR_PARAMS$f_stat_min, , drop = FALSE]
  if (nrow(exp) == 0) next

  or <- run_tsmr_one(exp, BCAC_sus_2020,
                     exposure_id = shorti,
                     outcome_label = "Breast_cancer_BCAC",
                     sensitivity_all = TRUE)
  if (is.null(or)) next

  AllOR <- rbind(AllOR, or)
}

message("Protein-breast cancer tests completed: ", nrow(AllOR))

# ---- 2. Multiple-testing correction --------------------------------------
AllOR <- add_multiple_testing(AllOR)
save_rdata(AllOR, file = "protein_BC_TSMR.rdata")
write_table(AllOR, "protein_BC_TSMR.csv")

AllOR_nohet <- filter_robust_mr(AllOR, use_direction = TRUE)
save_rdata(AllOR, AllOR_nohet, file = "protein_BC_TSMR_nohet.rdata")
write_table(AllOR_nohet, "protein_BC_TSMR_nohet.csv")

message(sprintf("Proteins causally associated with BC risk after robustness filter: %d",
                nrow(AllOR_nohet)))

# ---- 3. Manhattan plot ---------------------------------------------------
load(PATHS$generef_grch38)   # -> generef (gene, chromosome, start)

AllOR_nohet$CHR <- generef$chromosome[match(AllOR_nohet$exposure, generef$gene)]
AllOR_nohet$BP  <- generef$start[match(AllOR_nohet$exposure, generef$gene)]

dat <- AllOR_nohet %>%
  dplyr::filter(!is.na(CHR), !is.na(BP)) %>%
  dplyr::select(Protein = exposure, CHR, BP, P = pval)

prosig <- unique(c(dat$Protein[dat$P < 0.0005], AllOR_nohet$exposure))
prosig <- prosig[!is.na(prosig)]

pdf(file.path(OUT_FIGURES, "20_manhattan_protein_BC.pdf"), height = 5, width = 8)
CMplot::CMplot(dat, plot.type = "m", threshold = 0.05, cex = 0.5,
               amplify = TRUE,
               highlight = prosig, highlight.text = prosig,
               highlight.text.col = "black",
               chr.den.col = c("#1F78B4", "#f8c120"),
               file.output = FALSE,
               main = "Proteins causally associated with BC risk")
dev.off()

# ---- 4. Enrichment -------------------------------------------------------
pathway <- enrich_go_kegg(unique(AllOR_nohet$exposure))
plot_enrichment(pathway,
                file.path(OUT_FIGURES, "20_enrichment_protein_BC.pdf"),
                pal = c("#4197d8", "#acd372", "#f8c120", "#BC80BD"),
                width = 9.5, height = 9.5)

message("Step 20 complete.")

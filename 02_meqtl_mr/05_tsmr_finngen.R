# ===========================================================================
#  05_tsmr_finngen.R
#
#  Purpose
#  -------
#  Two-sample Mendelian randomisation of blood DNA methylation (exposure)
#  against breast-cancer risk (outcome) using FinnGen release 11 as the
#  outcome GWAS.
#
#      exposure : GoDMC whole-blood cis-meQTL (one instrument set per CpG)
#      outcome  : FinnGen R11 C3_BREAST_EXALLC, 20,586 cases / 201,494 controls
#
#  Input
#  -----
#    output/objects/GoDMC_as_exp/*.rdata        instrument sets (step 04)
#    PATHS$finngen_bc_outcome                   object: finndata
#
#  Output
#  ------
#    output/objects/TSMR_FinnGen.rdata                   AllOR
#    output/objects/TSMR_FinnGen_nohet.rdata             AllOR_nohet
#    output/tables/TSMR_FinnGen.csv
#    output/tables/TSMR_FinnGen_nohet.csv
#    output/figures/05_volcano_FinnGen.pdf
#    output/figures/05_manhattan_FinnGen.pdf
#
#  Methods
#  -------
#    * Primary estimate: inverse-variance weighted (IVW) for CpGs with >= 2
#      instruments, Wald ratio for a single instrument
#      (GagnonMR::primary_MR_analysis).
#    * Sensitivity: Cochran's Q for heterogeneity, MR-Egger intercept for
#      horizontal pleiotropy, Steiger directionality test.
#    * Multiple testing: Benjamini-Hochberg FDR and Bonferroni.
#    * Robustness filter applied to define the "no-het" result set:
#      P < 0.05 AND (Q P > 0.05 AND Egger intercept P > 0.05, or no
#      heterogeneity test applicable) AND correct causal direction.
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(tidyverse)
  library(data.table)
})

source("config/config.R")
source("helpers/mr_functions.R")

require_files(PATHS$finngen_bc_outcome)
load(PATHS$finngen_bc_outcome)      # -> finndata (TwoSampleMR outcome)

cpg_list <- sub("\\.rdata$", "",
                list.files(PATHS$godmc_exposure_dir, pattern = "\\.rdata$"))
message("Instrumented CpGs available: ", length(cpg_list))

# ---- 1. CpG-by-CpG MR -----------------------------------------------------
AllOR <- data.frame()

for (k in seq_along(cpg_list)) {
  cpg <- cpg_list[k]
  if (k %% 250 == 0) message("  ... ", k, " / ", length(cpg_list))

  exp <- load_godmc_instruments(cpg)
  if (is.null(exp)) next

  or <- run_tsmr_one(exp, finndata, exposure_id = cpg,
                     outcome_label = "BC (FinnGen)")
  if (is.null(or)) next

  AllOR <- rbind(AllOR, or)
}

message("CpG-outcome tests completed: ", nrow(AllOR))

# ---- 2. Multiple-testing correction --------------------------------------
AllOR <- add_multiple_testing(AllOR)

save_rdata(AllOR, file = "TSMR_FinnGen.rdata")
write_table(AllOR, "TSMR_FinnGen.csv")

# ---- 3. Robustness filter ------------------------------------------------
AllOR_nohet <- filter_robust_mr(AllOR, use_direction = TRUE)
save_rdata(AllOR_nohet, file = "TSMR_FinnGen_nohet.rdata")
write_table(AllOR_nohet, "TSMR_FinnGen_nohet.csv")

message(sprintf("Nominal (P < 0.05): %d | after heterogeneity / pleiotropy / directionality filter: %d",
                sum(AllOR$pval < 0.05, na.rm = TRUE), nrow(AllOR_nohet)))

# ---- 4. Volcano plot -----------------------------------------------------
plot_dat <- AllOR
plot_dat$Effect <- ifelse(plot_dat$FDR < 0.05 & plot_dat$b > 0, "Positive",
                   ifelse(plot_dat$FDR < 0.05 & plot_dat$b < 0, "Negative",
                          "Insignificant"))

p_volcano <- ggplot(plot_dat, aes(x = b, y = -log10(FDR), colour = Effect)) +
  geom_point(aes(size = abs(b)), alpha = 0.8) +
  scale_color_manual(values = c("Insignificant" = "#B3B3B3",
                                "Positive"      = "#1F78B4",
                                "Negative"      = "#f8c120")) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  theme_classic() +
  theme(panel.grid   = element_blank(),
        legend.title = element_text(size = 6.5),
        legend.text  = element_text(size = 6.5)) +
  labs(x = "Beta (effect size)",
       y = parse(text = "-log[10]*(FDR)"),
       title = "All MR results between CpG sites and BC risk (FinnGen)")

ggsave(file.path(OUT_FIGURES, "05_volcano_FinnGen.pdf"),
       p_volcano, width = 6.5, height = 4.5)

# ---- 5. Manhattan plot ---------------------------------------------------
plot_dat <- add_probe_coordinates(plot_dat)
manh <- plot_dat %>%
  dplyr::transmute(CpG = exposure, CHR, BP, P = pval) %>%
  dplyr::filter(!is.na(CHR), !is.na(BP))

prosig <- manh$CpG[manh$P < 1e-4]
prosig <- unique(c(prosig, AllOR_nohet$exposure))
prosig <- prosig[!is.na(prosig)]

pdf(file.path(OUT_FIGURES, "05_manhattan_FinnGen.pdf"), height = 5, width = 8)
CMplot::CMplot(manh, plot.type = "m", threshold = 0.05, cex = 0.5,
               amplify = TRUE,
               highlight = prosig, highlight.text = prosig,
               highlight.text.col = "black",
               chr.den.col = c("#1F78B4", "#f8c120"),
               file.output = FALSE,
               main = "CpGs causally associated with BC risk (FinnGen)")
dev.off()

message("Step 05 complete.")

# ===========================================================================
#  06_tsmr_bcac.R
#
#  Purpose
#  -------
#  Two-sample MR of blood DNA methylation against breast-cancer risk using the
#  BCAC 2020 susceptibility meta-analysis as the outcome GWAS.
#
#      exposure : GoDMC whole-blood cis-meQTL (one instrument set per CpG)
#      outcome  : BCAC 2020 overall breast cancer susceptibility
#                 (iCOGS + OncoArray + UK Biobank, N = 247,173)
#
#  Input
#  -----
#    output/objects/GoDMC_as_exp/*.rdata     instrument sets (step 04)
#    PATHS$bcac_sus_meta                     object: BCAC_sus_2020
#
#  Output
#  ------
#    output/objects/TSMR_BCAC.rdata / TSMR_BCAC_nohet.rdata
#    output/tables/TSMR_BCAC.csv  / TSMR_BCAC_nohet.csv
#    output/figures/06_volcano_BCAC.pdf
#    output/figures/06_manhattan_BCAC.pdf
#
#  Methods, parameters and the robustness filter are identical to step 05;
#  see the header of 05_tsmr_finngen.R.
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(tidyverse)
  library(data.table)
})

source("config/config.R")
source("helpers/mr_functions.R")

require_files(PATHS$bcac_sus_meta)
load(PATHS$bcac_sus_meta)           # -> BCAC_sus_2020 (TwoSampleMR outcome)

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

  or <- run_tsmr_one(exp, BCAC_sus_2020, exposure_id = cpg,
                     outcome_label = "BC_BCAC")
  if (is.null(or)) next

  AllOR <- rbind(AllOR, or)
}

message("CpG-outcome tests completed: ", nrow(AllOR))

# ---- 2. Multiple-testing correction --------------------------------------
AllOR <- add_multiple_testing(AllOR)

save_rdata(AllOR, file = "TSMR_BCAC.rdata")
write_table(AllOR, "TSMR_BCAC.csv")

# ---- 3. Robustness filter ------------------------------------------------
AllOR_nohet <- filter_robust_mr(AllOR, use_direction = TRUE)
save_rdata(AllOR_nohet, file = "TSMR_BCAC_nohet.rdata")
write_table(AllOR_nohet, "TSMR_BCAC_nohet.csv")

message(sprintf("Nominal (P < 0.05): %d | after robustness filter: %d",
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
       title = "All MR results between CpG sites and BC risk (BCAC)")

ggsave(file.path(OUT_FIGURES, "06_volcano_BCAC.pdf"),
       p_volcano, width = 6.5, height = 4.5)

# ---- 5. Manhattan plot ---------------------------------------------------
plot_dat <- add_probe_coordinates(plot_dat)
manh <- plot_dat %>%
  dplyr::transmute(CpG = exposure, CHR, BP, P = pval) %>%
  dplyr::filter(!is.na(CHR), !is.na(BP))

prosig <- unique(c(manh$CpG[manh$P < 2e-5], AllOR_nohet$exposure))
prosig <- prosig[!is.na(prosig)]

pdf(file.path(OUT_FIGURES, "06_manhattan_BCAC.pdf"), height = 5, width = 8)
CMplot::CMplot(manh, plot.type = "m", threshold = 0.05, cex = 0.5,
               amplify = TRUE,
               highlight = prosig, highlight.text = prosig,
               highlight.text.col = "black",
               chr.den.col = c("#1F78B4", "#f8c120"),
               file.output = FALSE,
               main = "CpGs causally associated with BC risk (BCAC)")
dev.off()

message("Step 06 complete.")

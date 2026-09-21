# ===========================================================================
#  04_build_godmc_instruments.R
#
#  Purpose
#  -------
#  Build the genetic instrument set for every differentially methylated
#  position (DMP) that is present in GoDMC.  One exposure object is written
#  per CpG; all three TSMR steps (05, 06, 07) and the phenome-wide screen
#  (step 17) re-use these files.
#
#  Input
#  -----
#    output/objects/All_clean_MeQTL_GoDMC.rdata   object: data2
#    output/objects/GSE104942_DMP_Diff.rdata      object: dmpDiff
#    PLINK_BFILE, PLINK_BIN                       1000 Genomes EUR reference
#
#  Output
#  ------
#    output/objects/GoDMC_as_exp/<cpg>.rdata      object: exp  (exposure)
#    output/tables/04_instrument_summary.csv      per-CpG instrument counts
#
#  Instrument definition
#  ---------------------
#    * cis-meQTL SNPs at P < 5e-8 for the CpG
#    * LD clumped (r^2 < 0.001, 10,000 kb window, 1000 Genomes EUR)
#    * per-SNP F statistic > 10
#
#  Note
#  ----
#  CpGs for which no instrument survives these filters are recorded as
#  "not instrumentable" in the summary table and are skipped downstream.
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(tidyverse)
  library(data.table)
})

source("config/config.R")
source("helpers/mr_functions.R")

load(file.path(OUT_OBJECTS, "All_clean_MeQTL_GoDMC.rdata"))   # data2
load(file.path(OUT_OBJECTS, "GSE104942_DMP_Diff.rdata"))      # dmpDiff

instrument_dir <- file.path(OUT_OBJECTS, "GoDMC_as_exp")
dir.create(instrument_dir, showWarnings = FALSE, recursive = TRUE)

# Only CpGs that are (a) significant DMPs and (b) present in GoDMC can be
# instrumented.
dmp_cpg    <- rownames(dmpDiff)
testable   <- intersect(dmp_cpg, unique(data2$cpg))
message(sprintf("DMPs: %d | present in GoDMC: %d", length(dmp_cpg), length(testable)))

summary_rows <- vector("list", length(testable))
n_failed     <- 0L

for (k in seq_along(testable)) {
  cpg <- testable[k]
  if (k %% 500 == 0) message("  ... ", k, " / ", length(testable))

  exp <- try(build_godmc_instruments(cpg, data2, quiet = TRUE), silent = TRUE)

  if (inherits(exp, "try-error") || is.null(exp)) {
    n_failed <- n_failed + 1L
    summary_rows[[k]] <- data.frame(cpg = cpg, n_instruments = 0L,
                                    min_f = NA_real_, instrumented = FALSE)
    next
  }

  save(exp, file = file.path(instrument_dir, paste0(cpg, ".rdata")))
  summary_rows[[k]] <- data.frame(
    cpg = cpg,
    n_instruments = nrow(exp),
    min_f = min(exp$fval),
    instrumented = TRUE)
}

inst_summary <- do.call(rbind, summary_rows)
write_table(inst_summary, "04_instrument_summary.csv")

message(sprintf("Instrumented CpGs: %d | not instrumentable: %d",
                sum(inst_summary$instrumented), n_failed))
message(sprintf("Median number of instruments per CpG: %.0f",
                median(inst_summary$n_instruments[inst_summary$instrumented])))

# ===========================================================================
#  01_dmp_discovery_GSE104942.R
#
#  Purpose
#  -------
#  Stage 1 of the pipeline: identify differentially methylated positions
#  (DMPs) between breast-cancer cases and unaffected controls in blood, using
#  Illumina 450K data from GEO series GSE104942.  The significant CpGs become
#  the exposure set for the cis-meQTL Mendelian randomisation screen.
#
#  Input
#  -----
#    PATHS$gse104942_matrix  GSE104942_series_matrix.txt.gz
#
#  Output
#  ------
#    output/objects/GSE104942_norm.rdata       matData (beta matrix), clinical
#    output/objects/GSE104942_DMP_Diff.rdata   dmpDiff  (FDR < 0.05)
#    output/tables/GSE104942_DMP_Diff.csv
#    output/figures/01_*Box.pdf, 01_densityBeanPlot.pdf, 01_mdsPlot.pdf
#    output/qc/01_GSE104942_QC_summary.csv
#
#  Processing steps
#  ----------------
#    1. k-nearest-neighbour imputation of missing beta values (impute.knn)
#    2. removal of probes with a mean beta below 0.005 (unreliable signal)
#    3. between-array quantile normalisation (wateRmelon::betaqn)
#    4. array-level QC plots (boxplot, density beanplot, MDS)
#    5. differential methylation testing on M-values (minfi::dmpFinder,
#       categorical comparison, FDR controlled by Benjamini-Hochberg)
#
#  Probe annotation
#  ----------------
#    makeGenomicRatioSetFromMatrix() defaults to the IlluminaHumanMethylation450k
#    / ilmn12.hg19 annotation, which matches the platform used by GSE104942.
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(GEOquery)
  library(minfi)
  library(impute)
  library(wateRmelon)
})

source("config/config.R")
require_files(PATHS$gse104942_matrix)

set.seed(20240101)

# ---- 1. Read the series matrix and build the phenotype table --------------
rawfile <- getGEO(filename = PATHS$gse104942_matrix, getGPL = FALSE)
info    <- rawfile@phenoData@data

clinical <- info %>%
  dplyr::select(characteristics_ch1) %>%
  mutate(Sample = rownames(info),
         Group  = factor(ifelse(characteristics_ch1 == "status: unaffected",
                                "Control", "Cancer"),
                         levels = c("Control", "Cancer")))

message("Samples per group:")
print(table(clinical$Group))

# phenoData row order is identical to the assay column order, so `clinical`
# can be used directly as the sample annotation for the beta matrix.
stopifnot(identical(rownames(clinical), colnames(rawfile@assayData[["exprs"]])))

# ---- 2. k-NN imputation of missing beta values ----------------------------
exp <- rawfile@assayData[["exprs"]]

mat <- impute.knn(exp, k = DMP_PARAMS$impute_knn_k)
matData <- mat$data
matData <- matData + 1e-5      # avoid exact zeros before log-type transforms

# ---- 3. Filter and normalise ---------------------------------------------
# Probes with a very low mean beta carry little signal and behave poorly
# during quantile normalisation.
n_before <- nrow(matData)
matData  <- matData[rowMeans(matData) > DMP_PARAMS$detection_p, ]
message(sprintf("Probes retained: %d of %d (mean beta > %.3f)",
                nrow(matData), n_before, DMP_PARAMS$detection_p))

# Array-level QC before normalisation
pdf(file.path(OUT_FIGURES, "01_rawBox.pdf"))
boxplot(matData, col = "#4197d8", xaxt = "n", outline = FALSE)
dev.off()

matData <- betaqn(matData)     # between-array quantile normalisation

# Array-level QC after normalisation
pdf(file.path(OUT_FIGURES, "01_normalBox.pdf"))
boxplot(matData, col = "#f8c120", xaxt = "n", outline = FALSE)
dev.off()

pdf(file.path(OUT_FIGURES, "01_densityBeanPlot.pdf"))
par(oma = c(2, 10, 2, 2))
densityBeanPlot(matData, sampGroups = clinical$Group, sampNames = clinical$Sample)
dev.off()

pdf(file.path(OUT_FIGURES, "01_mdsPlot.pdf"))
mdsPlot(matData, numPositions = 1000,
        sampGroups = clinical$Group, sampNames = clinical$Sample)
dev.off()

# ---- 4. Differential methylation analysis ---------------------------------
# dmpFinder is run on M-values (logit of beta), which are better behaved for
# linear modelling than beta values.
grset <- makeGenomicRatioSetFromMatrix(matData, what = "Beta")
M     <- getM(grset)

dmp     <- dmpFinder(M, pheno = clinical$Group, type = "categorical")
dmpDiff <- dmp[!is.na(dmp$qval) & dmp$qval < DMP_PARAMS$dmp_qvalue_cutoff, ]

message(sprintf("DMPs at FDR < %.2f: %d", DMP_PARAMS$dmp_qvalue_cutoff, nrow(dmpDiff)))

# ---- 5. Save --------------------------------------------------------------
save_rdata(matData, clinical, file = "GSE104942_norm.rdata")
save_rdata(dmpDiff,           file = "GSE104942_DMP_Diff.rdata")

fwrite(dmpDiff, file.path(OUT_TABLES, "GSE104942_DMP_Diff.csv"), row.names = TRUE)

qc_summary <- data.frame(
  n_samples       = ncol(matData),
  n_control       = sum(clinical$Group == "Control"),
  n_cancer        = sum(clinical$Group == "Cancer"),
  n_probes_input  = n_before,
  n_probes_kept   = nrow(matData),
  n_dmp_fdr005    = nrow(dmpDiff)
)
write_table(qc_summary, "01_GSE104942_QC_summary.csv")
print(qc_summary)

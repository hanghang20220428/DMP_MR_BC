# ===========================================================================
#  22_phenogram_input.R
#
#  Purpose
#  -------
#  Build the input tables for PhenoGram (https://visualization.ritchielab.org/
#  phenograms/plot), which draws each causal CpG and its significant mediator
#  proteins on a chromosome ideogram.
#
#  Input
#  -----
#    output/objects/mediation_results.rdata   sig_mediation
#    PATHS$anno_450k, PATHS$anno_epic         probe coordinates (GRCh37)
#    PATHS$generef_grch37                     gene coordinates (GRCh37)
#
#  Output
#  ------
#    output/tables/22_phenogram_<cpg>.txt     tab-separated, one file per CpG
#        cpg_pro   chr   pos   phenotype
#
#  Notes
#  -----
#    * The first column holds the feature identifier (a CpG probe or a gene
#      symbol) and is named `cpg_pro` for compatibility with the files
#      originally submitted to PhenoGram; PhenoGram maps columns by name or
#      position at upload time.
#    * Output files are written with LF line endings; PhenoGram expects
#      position and chromosome in GRCh37/hg19 coordinates.
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

source("config/config.R")
source("helpers/mr_functions.R")   # read_probe_manifest()

load(file.path(OUT_OBJECTS, "mediation_results.rdata"))   # all_mediation, sig_mediation
load(PATHS$generef_grch37)                                # -> generef

anno_450k <- read_probe_manifest(PATHS$anno_450k)
anno_epic <- read_probe_manifest(PATHS$anno_epic)

cpg_position <- function(cpg) {
  chr <- anno_epic$CHR[match(cpg, anno_epic$ID)]
  pos <- anno_epic$MAPINFO[match(cpg, anno_epic$ID)]
  chr450 <- anno_450k$CHR[match(cpg, anno_450k$ID)]
  pos450 <- anno_450k$MAPINFO[match(cpg, anno_450k$ID)]
  list(chr = ifelse(is.na(chr450), chr, chr450),
       pos = ifelse(is.na(pos450), pos, pos450))
}

for (i in unique(sig_mediation$exposure)) {

  mediators <- unique(sig_mediation$mediator[sig_mediation$exposure == i])

  # (a) the CpG itself
  p <- cpg_position(i)
  data1 <- data.frame(cpg_pro = i, chr = p$chr, pos = p$pos, phenotype = i)

  # (b) the mediator proteins, placed at their gene start coordinate
  data2 <- data.frame(
    cpg_pro   = mediators,
    chr       = generef$chromosome[match(mediators, generef$gene)],
    pos       = generef$start[match(mediators, generef$gene)],
    phenotype = mediators)

  out <- rbind(data1, data2)

  write.table(out, file.path(OUT_TABLES, paste0("22_phenogram_", i, ".txt")),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

  message(sprintf("  %s: %d features (1 CpG + %d proteins)",
                  i, nrow(out), length(mediators)))
}

message("Step 22 complete.")

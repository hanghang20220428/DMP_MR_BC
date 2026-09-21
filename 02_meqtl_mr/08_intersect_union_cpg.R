# ===========================================================================
#  08_intersect_union_cpg.R
#
#  Purpose
#  -------
#  Define the CpG sets that are carried forward to the GSMR sensitivity
#  analysis and to the downstream colocalisation / mediation work.
#
#      set A : CpGs significant in BOTH FinnGen and BCAC          (intersection)
#      set B : CpGs significant in EITHER FinnGen or BCAC         (union)
#      set C : CpGs significant in FinnGen AND BCAC AND GCST90011804
#
#  Set A (two-outcome intersection) is the primary causal set: it defines the
#  CpGs taken into GSMR (steps 09-10) and, after meta-analysis (step 13),
#  into colocalisation (steps 14-16) and the protein mediation analysis
#  (steps 19-24).  Sets B and C are reported as sensitivity checks.
#
#  Input
#  -----
#    output/objects/TSMR_FinnGen_nohet.rdata          AllOR_nohet
#    output/objects/TSMR_BCAC_nohet.rdata             AllOR_nohet
#    output/objects/TSMR_GCST90011804_nohet.rdata     AllOR_nohet
#
#  Output
#  ------
#    output/objects/CpG_sets_meqtl_MR.rdata
#        cpg_intersect_finn_bcac        (set A)
#        cpg_union_finn_bcac            (set B)
#        cpg_intersect_finn_bcac_gcst   (set C)
#    output/tables/08_CpG_sets_summary.csv
# ===========================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
})

source("config/config.R")

load(file.path(OUT_OBJECTS, "TSMR_FinnGen_nohet.rdata"))        # AllOR_nohet
finn <- AllOR_nohet

load(file.path(OUT_OBJECTS, "TSMR_BCAC_nohet.rdata"))           # AllOR_nohet
bcac <- AllOR_nohet

load(file.path(OUT_OBJECTS, "TSMR_GCST90011804_nohet.rdata"))   # AllOR_nohet
gcst <- AllOR_nohet

# ---- CpG sets -------------------------------------------------------------
cpg_intersect_finn_bcac      <- Reduce(intersect, list(finn$exposure, bcac$exposure))
cpg_union_finn_bcac          <- Reduce(union,     list(finn$exposure, bcac$exposure))
cpg_intersect_finn_bcac_gcst <- Reduce(intersect, list(finn$exposure, bcac$exposure,
                                                       gcst$exposure))

save_rdata(cpg_intersect_finn_bcac,
           cpg_union_finn_bcac,
           cpg_intersect_finn_bcac_gcst,
           file = "CpG_sets_meqtl_MR.rdata")

summary_tab <- data.frame(
  set = c("intersection_FinnGen_BCAC",
          "union_FinnGen_BCAC",
          "intersection_FinnGen_BCAC_GCST90011804"),
  n_cpg = c(length(cpg_intersect_finn_bcac),
            length(cpg_union_finn_bcac),
            length(cpg_intersect_finn_bcac_gcst)),
  cpg_list = c(paste(cpg_intersect_finn_bcac, collapse = "; "),
               paste(cpg_union_finn_bcac, collapse = "; "),
               paste(cpg_intersect_finn_bcac_gcst, collapse = "; "))
)
write_table(summary_tab, "08_CpG_sets_summary.csv")
print(summary_tab[, 1:2])

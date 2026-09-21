# ===========================================================================
#  config/paths_local.R   --   MACHINE-SPECIFIC, NOT TRACKED BY GIT
#
#  Local data locations for this workstation.  Values here override the
#  defaults in config/config.R.
# ===========================================================================

# ---- Reference / software ------------------------------------------------
plink_bfile <- "D:/1kg.v3/EUR"
plink_bin   <- "D:/Program files/plink_win64_20230116/plink.exe"

# ---- Discovery cohort (GSE104942) ---------------------------------------
gse104942_matrix <- "K:/4DMP_MR_BC/GSE104942_series_matrix.txt.gz"

# ---- GoDMC cis-meQTL -----------------------------------------------------
godmc_assoc        <- "K:/4DMP_MR_BC/assoc_meta_all.csv.gz"
godmc_snps         <- "K:/4DMP_MR_BC/snps.csv.gz"
godmc_exposure_dir <- "K:/4DMP_MR_BC/GoDMC_as_exp"

# ---- FinnGen R11 ---------------------------------------------------------
finngen_bc_raw       <- "K:/4DMP_MR_BC/Finn_R11_BC/finngen_R11_C3_BREAST_EXALLC.gz"
finngen_manifest     <- "K:/4DMP_MR_BC/Finn_R11_BC/finngen_R11_manifest.csv"
finngen_bc_outcome   <- "K:/6Ovarian_BC_Comorbidity/Finn_R11_BC/C3_BREAST_EXALLC_finndataR11_as_outcome.rdata"
finngen_bc_coloc     <- "K:/4DMP_MR_BC/Finn_R11_BC/finndataR11_BC_as_coloc.rdata"
finngen_all_outcomes <- "H:/2自己数据/01GWAS数据/FinnR11_all_outcome_rdata"

# ---- BCAC ----------------------------------------------------------------
bcac_sus_meta  <- "M:/0BC_data_BCAC/1BCAC_sus_2020_meta.rdata"
bcac_sus_raw   <- "M:/0BC_data_BCAC/susceptibility_2020/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt"
bcac_sus_coloc <- "M:/0BC_data_BCAC/2BCAC_2020_suscep_as_coloc.rdata"

# ---- Third breast-cancer outcome (GCST90011804) -------------------------
gcst90011804_raw     <- "K:/4DMP_MR_BC/泛癌17NC2020/32887889-GCST90011804-EFO_0000305-Build37.f.tsv.gz"
gcst90011804_outcome <- "K:/4DMP_MR_BC/泛癌17NC2020/GCST90011804_as_outcome.rdata"
pancancer_info       <- "K:/4DMP_MR_BC/泛癌17NC2020/17泛癌信息.rdata"
pancancer_dir        <- "K:/4DMP_MR_BC/泛癌17NC2020"
gwas_catalog_dir     <- "K:/4DMP_MR_BC/泛癌17NC2020"
snp_eaf_eur          <- "M:/0SNP_eaf_ref/0snp_eaf_eur_ALFA.rdata"

# ---- UKB-PPP plasma proteomics ------------------------------------------
ukbppp_outcome_dir  <- "H:/2自己数据/pQTL数据/2940UK_pQTL_Outcome"
ukbppp_exposure_dir <- "H:/2自己数据/pQTL数据/2940UK_pQTL_Exposure6SNP"

# ---- Single-cell (GSE161529) --------------------------------------------
# NOTE: the raw 10X matrices are not present on this workstation any more.
# Re-download them from GEO into this directory before re-running step 25,
# or point this at the existing location.
scrna_10x_root <- "K:/4DMP_MR_BC/GSE161529_10X"

# ---- Annotation ----------------------------------------------------------
anno_450k      <- "K:/4DMP_MR_BC/GPL13534-11288.txt"
anno_epic      <- "K:/4DMP_MR_BC/GPL21145-48548.txt"
generef_grch37 <- "K:/4DMP_MR_BC/GRCh37.87.geneloc_ref.rdata"
generef_grch38 <- "K:/6Ovarian_BC_Comorbidity/GRCh38.100.geneloc_ref.rdata"
epic_meqtl     <- "K:/4DMP_MR_BC/2.5EPIC_meQTL.rdata"

# ===========================================================================
#  config/config.R
#
#  Central configuration for the 4DMP-MR-BC analysis pipeline
#  (DNA-methylation-driven causal and druggable discovery in breast cancer).
#
#  Every analysis script in this repository starts with:
#
#      source("config/config.R")
#
#  and then refers to objects defined here instead of hard-coded absolute
#  paths.  This keeps the public code portable while still allowing the
#  authors to run it unchanged on their own machine.
#
#  Path resolution order (first match wins)
#  ----------------------------------------
#    1. an environment variable (e.g. DMPMR_BCAC_DIR)
#    2. a variable of the same name in config/paths_local.R   (git-ignored)
#    3. the default shipped with the repository
#
#  To run the pipeline on your own machine:
#    cp config/paths_local.R.example config/paths_local.R   and edit it.
# ===========================================================================

# ---------------------------------------------------------------------------
# 0. Locate the repository root
# ---------------------------------------------------------------------------
#  The scripts are meant to be run with the repository root as the working
#  directory (open 4DMP-MR-BC.Rproj in RStudio, or setwd() to the clone).
#  DMPMR_BC_ROOT can be used to point at the root explicitly.

.find_repo_root <- function() {
  root <- Sys.getenv("DMPMR_BC_ROOT", unset = "")
  if (nzchar(root) && dir.exists(root)) {
    return(normalizePath(root, winslash = "/", mustWork = FALSE))
  }
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, "config", "config.R"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

REPO_ROOT <- .find_repo_root()

# ---------------------------------------------------------------------------
# 1. User-supplied local paths (optional, not tracked by git)
# ---------------------------------------------------------------------------

# The parent must be baseenv(), not emptyenv(): the assignments inside
# paths_local.R are ordinary calls to `<-`, which R looks up like any other
# function, so an emptyenv() parent makes every assignment fail with
# "could not find function \"<-\"".
.dmpmr_local <- new.env(parent = baseenv())

#' Load config/paths_local.R into an environment, robustly.
#'
#' `source(encoding = "UTF-8")` is the natural way to read a file that may
#' contain non-ASCII directory names (several data roots do), but on Windows
#' builds whose native locale is not UTF-8 -- e.g. `LC_CTYPE=C` or `CP936` --
#' it attempts an iconv translation into the native charset, fails on the
#' first non-ASCII byte and leaves the string literal unterminated, which
#' raises a parse error.  The fallback below reads the file as UTF-8 text and
#' evaluates the parsed expressions directly, which preserves the bytes.
.load_local_paths <- function(path, env) {
  ok <- tryCatch({
    source(path, local = env, encoding = "UTF-8")
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)

  if (!ok) {
    # `source()` may have assigned part of the file before failing; clear the
    # environment so the fallback starts from a clean state.
    rm(list = ls(env, all.names = TRUE), envir = env)
    lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
    eval(parse(text = lines, encoding = "UTF-8"), envir = env)
  }
  invisible(TRUE)
}

local_paths_file <- file.path(REPO_ROOT, "config", "paths_local.R")
if (file.exists(local_paths_file)) {
  .load_local_paths(local_paths_file, .dmpmr_local)
  message("[config] local paths loaded from config/paths_local.R")
}

# Resolve one path: env var -> paths_local.R -> default
.resolve_path <- function(name, default) {
  env_name <- paste0("DMPMR_", toupper(name))
  v <- Sys.getenv(env_name, unset = "")
  if (nzchar(v)) return(v)
  if (exists(name, envir = .dmpmr_local, inherits = FALSE)) {
    v <- get(name, envir = .dmpmr_local, inherits = FALSE)
    if (!is.null(v) && length(v) == 1L && nzchar(v)) return(v)
  }
  default
}

# ---------------------------------------------------------------------------
# 2. Output directories
# ---------------------------------------------------------------------------

OUTPUT_DIR   <- file.path(REPO_ROOT, "output")
OUT_TABLES   <- file.path(OUTPUT_DIR, "tables")
OUT_FIGURES  <- file.path(OUTPUT_DIR, "figures")
OUT_OBJECTS  <- file.path(OUTPUT_DIR, "objects")
OUT_QC       <- file.path(OUTPUT_DIR, "qc")
for (.d in c(OUTPUT_DIR, OUT_TABLES, OUT_FIGURES, OUT_OBJECTS, OUT_QC)) {
  dir.create(.d, showWarnings = FALSE, recursive = TRUE)
}
rm(.d)

# Convenience writers -------------------------------------------------------
#  save_rdata(): keep intermediate R objects out of the repository root
save_rdata <- function(..., file) {
  save(..., file = file.path(OUT_OBJECTS, file))
}
#  write_table(): always write a CSV plus an RDS copy of the same table
write_table <- function(x, file, row.names = FALSE) {
  utils::write.csv(x, file.path(OUT_TABLES, file), row.names = row.names)
  saveRDS(x, file.path(OUT_TABLES, sub("\\.csv$", ".rds", file)))
}

# ---------------------------------------------------------------------------
# 3. External reference / software locations
# ---------------------------------------------------------------------------
#  These are large third-party resources that are not redistributed with the
#  code.  See docs/03_data_sources.md for download instructions.

PLINK_BFILE <- .resolve_path("plink_bfile", "data/1kg.v3/EUR")          # 1000 Genomes EUR, PLINK 1 binary
PLINK_BIN   <- .resolve_path("plink_bin",   "plink")                    # PLINK 1.9 executable

# ---------------------------------------------------------------------------
# 4. Source datasets
# ---------------------------------------------------------------------------
#  Every entry documents (a) what the file is and (b) where it came from.
#  Defaults are relative to DATA_ROOT; set paths_local.R for local copies.

DATA_ROOT <- .resolve_path("data_root", file.path(REPO_ROOT, "data"))

PATHS <- list(

  # ---- 4.1 Discovery cohort: blood DNA methylation (Illumina 450K) --------
  # GEO GSE104942 series matrix (breast-cancer cases vs unaffected controls)
  gse104942_matrix = .resolve_path("gse104942_matrix",
                                   file.path(DATA_ROOT, "GSE104942",
                                             "GSE104942_series_matrix.txt.gz")),

  # ---- 4.2 GoDMC blood cis-meQTL meta-analysis (GRCh37) -------------------
  # GoDMC consortium, mQTL meta-analysis of whole-blood samples.
  godmc_assoc      = .resolve_path("godmc_assoc",
                                   file.path(DATA_ROOT, "GoDMC", "assoc_meta_all.csv.gz")),
  # SNP annotation table used to map GoDMC "chr:pos" identifiers to rsIDs
  godmc_snps       = .resolve_path("godmc_snps",
                                   file.path(DATA_ROOT, "GoDMC", "snps.csv.gz")),
  # Directory holding one <CpG>.rdata exposure file per instrumented CpG
  godmc_exposure_dir = .resolve_path("godmc_exposure_dir",
                                     file.path(DATA_ROOT, "GoDMC", "GoDMC_as_exp")),

  # ---- 4.3 Outcome GWAS ---------------------------------------------------
  # FinnGen release 11, C3_BREAST_EXALLC (20,586 cases / 201,494 controls)
  finngen_bc_raw      = .resolve_path("finngen_bc_raw",
                                      file.path(DATA_ROOT, "FinnGen_R11", "finngen_R11_C3_BREAST_EXALLC.gz")),
  finngen_manifest    = .resolve_path("finngen_manifest",
                                      file.path(DATA_ROOT, "FinnGen_R11", "finngen_R11_manifest.csv")),
  finngen_bc_outcome  = .resolve_path("finngen_bc_outcome",
                                      file.path(DATA_ROOT, "FinnGen_R11", "C3_BREAST_EXALLC_finndataR11_as_outcome.rdata")),
  finngen_bc_coloc    = .resolve_path("finngen_bc_coloc",
                                      file.path(DATA_ROOT, "FinnGen_R11", "finndataR11_BC_as_coloc.rdata")),
  # Directory of FinnGen R11 endpoints already formatted as TwoSampleMR outcomes
  # (used by the phenome-wide MR screen, step 16)
  finngen_all_outcomes = .resolve_path("finngen_all_outcomes",
                                       file.path(DATA_ROOT, "FinnGen_R11", "all_outcome_rdata")),

  # BCAC 2020 susceptibility meta-analysis (247,173 cases/controls)
  bcac_sus_meta       = .resolve_path("bcac_sus_meta",
                                      file.path(DATA_ROOT, "BCAC", "1BCAC_sus_2020_meta.rdata")),
  bcac_sus_raw        = .resolve_path("bcac_sus_raw",
                                      file.path(DATA_ROOT, "BCAC",
                                                "icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt")),
  bcac_sus_coloc      = .resolve_path("bcac_sus_coloc",
                                      file.path(DATA_ROOT, "BCAC", "2BCAC_2020_suscep_as_coloc.rdata")),

  # Third breast-cancer GWAS used as an independent replication outcome
  # (GCST90011804, BCAC/BRIDGES/CIMBA, GRCh37, hg19)
  gcst90011804_raw    = .resolve_path("gcst90011804_raw",
                                      file.path(DATA_ROOT, "PanCancer17",
                                                "32887889-GCST90011804-EFO_0000305-Build37.f.tsv.gz")),
  gcst90011804_outcome = .resolve_path("gcst90011804_outcome",
                                       file.path(DATA_ROOT, "PanCancer17", "GCST90011804_as_outcome.rdata")),
  pancancer_info      = .resolve_path("pancancer_info",
                                      file.path(DATA_ROOT, "PanCancer17", "pancancer17_info.rdata")),

  # Ethnic-matched allele frequencies (ALFA) used to back-fill missing EAFs
  snp_eaf_eur         = .resolve_path("snp_eaf_eur",
                                      file.path(DATA_ROOT, "reference", "0snp_eaf_eur_ALFA.rdata")),

  # ---- 4.4 UKB-PPP plasma proteomics (2,940 aptamers) ---------------------
  # <protein>_as_outcome.rdata  : UKB-PPP proteins used as outcome
  ukbppp_outcome_dir  = .resolve_path("ukbppp_outcome_dir",
                                      file.path(DATA_ROOT, "UKB_PPP", "2940UK_pQTL_Outcome")),
  # <protein>_as_exposure.rdata : UKB-PPP proteins used as exposure (6-SNP panel)
  ukbppp_exposure_dir = .resolve_path("ukbppp_exposure_dir",
                                      file.path(DATA_ROOT, "UKB_PPP", "2940UK_pQTL_Exposure6SNP")),

  # ---- 4.5 Single-cell RNA-seq: GSE161529 ---------------------------------
  # 10X Genomics Chromium 3' libraries from normal breast and breast tumours
  # (Pal et al., EMBO J 2021).  The array below lists the 47 samples retained
  # here, grouped by receptor status.
  scrna_10x_root      = .resolve_path("scrna_10x_root",
                                      file.path(DATA_ROOT, "GSE161529")),
  # Seurat object produced by step 24/25 and reused by step 26
  scrna_seurat_celltype = .resolve_path("scrna_seurat_celltype",
                                        file.path(OUT_OBJECTS, "GSE161529_seurat_celltype.rdata")),

  # ---- 4.6 Genome annotation ---------------------------------------------
  # Illumina probe manifests, downloaded from GEO (GPL13534 = 450K,
  # GPL21145 = EPIC, GPL33022 = EPIC on GRCh38).  GRCh37/hg19 coordinates.
  anno_450k           = .resolve_path("anno_450k",
                                      file.path(DATA_ROOT, "manifest", "GPL13534-11288.txt")),
  anno_epic           = .resolve_path("anno_epic",
                                      file.path(DATA_ROOT, "manifest", "GPL21145-48548.txt")),
  # Gene coordinates (Ensembl release 87 / GRCh37, and GRCh38 for UKB-PPP)
  generef_grch37      = .resolve_path("generef_grch37",
                                      file.path(DATA_ROOT, "reference", "GRCh37.87.geneloc_ref.rdata")),
  generef_grch38      = .resolve_path("generef_grch38",
                                      file.path(DATA_ROOT, "reference", "GRCh38.100.geneloc_ref.rdata")),
  # EPIC meQTL summary statistics (probe -> chromosome/position lookup)
  epic_meqtl          = .resolve_path("epic_meqtl",
                                      file.path(DATA_ROOT, "reference", "EPIC_meQTL.rdata"))
)

# ---------------------------------------------------------------------------
# 5. GSE161529 sample definition
# ---------------------------------------------------------------------------
#  GEO series GSE161529 (Pal et al., EMBO J 2021): 10X Genomics Chromium 3'
#  libraries from normal breast tissue and breast tumours.  The 47 samples
#  retained here are listed in the manuscript; `group` is the receptor-status
#  label used throughout the analysis.

GSE161529_SAMPLES <- data.frame(
  gsm = c(
    # ---- Normal breast tissue (n = 13) -----------------------------------
    "GSM4909253", "GSM4909254", "GSM4909257", "GSM4909261", "GSM4909263",
    "GSM4909265", "GSM4909266", "GSM4909268", "GSM4909270", "GSM4909271",
    "GSM4909272", "GSM4909274", "GSM4909276",
    # ---- ER+ tumours (n = 19) --------------------------------------------
    "GSM4909296", "GSM4909297", "GSM4909298", "GSM4909299", "GSM4909300",
    "GSM4909301", "GSM4909302", "GSM4909303", "GSM4909304", "GSM4909305",
    "GSM4909306", "GSM4909307", "GSM4909309", "GSM4909311", "GSM4909313",
    "GSM4909315", "GSM4909317", "GSM4909319", "GSM4909320",
    # ---- HER2+ tumours (n = 6) -------------------------------------------
    "GSM4909289", "GSM4909290", "GSM4909291", "GSM4909292", "GSM4909293",
    "GSM4909294",
    # ---- PR+ tumour (n = 1) ----------------------------------------------
    "GSM4909295",
    # ---- TNBC tumours (n = 8) --------------------------------------------
    "GSM4909281", "GSM4909282", "GSM4909283", "GSM4909284", "GSM4909285",
    "GSM4909286", "GSM4909287", "GSM4909288"),
  group = c(
    rep("Normal", 13), rep("ER+", 19), rep("HER2+", 6), rep("PR+", 1),
    rep("TNBC", 8)),
  stringsAsFactors = FALSE
)

# Expected 10X directory layout under PATHS$scrna_10x_root:
#     <scrna_10x_root>/<GSM>/barcodes.tsv.gz
#     <scrna_10x_root>/<GSM>/features.tsv.gz
#     <scrna_10x_root>/<GSM>/matrix.mtx.gz

# ---------------------------------------------------------------------------
# 6. Analysis parameters
# ---------------------------------------------------------------------------
#  Named here so that every threshold quoted in the manuscript can be traced
#  back to a single line of code.

## ---- 6.1 Methylation preprocessing (step 01) ------------------------------
DMP_PARAMS <- list(
  impute_knn_k      = 10,        # impute.knn(): neighbours used for k-NN imputation
  detection_p       = 0.005,     # drop probes with mean beta < 0.005 before BMIQ
  normalize_method  = "betaqn",  # wateRmelon::betaqn quantile normalisation
  dmp_qvalue_cutoff = 0.05       # limma/dmpFinder FDR (Benjamini-Hochberg) cutoff
)

## ---- 6.2 Instrument selection for cis-meQTL MR (steps 04-06) --------------
MR_PARAMS <- list(
  clump_p            = 5e-8,     # genome-wide significance for instrument selection
  clump_r2           = 0.001,    # LD clumping r^2
  clump_kb           = 10000,    # LD clumping window (kb)
  f_stat_min         = 10,       # minimum per-SNP F statistic (weak-instrument filter)
  outcome_p_remove   = 1e-5,     # drop instruments associated with the outcome (P < 1e-5)
  harmonise_action   = 2,        # TwoSampleMR: infer palindromic SNP strand
  pleiotropy_p       = 0.05,     # MR-Egger intercept / MR-PRESSO global test
  heterogeneity_p    = 0.05,     # Cochran's Q
  fdr_method         = "BH",
  bonferroni_method  = "bonferroni"
)

## ---- 6.3 GSMR parameters (steps 08-10) ------------------------------------
##  IMPORTANT - read docs/04_known_issues.md before changing legacy mode.
##
##  The published GSMR runs called gsmr2::gsmr() with POSITIONAL arguments
##  written against the argument list of gsmr v1 (gwas_thresh immediately
##  followed by multi_snps_heidi_thresh).  gsmr2 1.1.1 inserts
##  `single_snp_heidi_thresh` between them, so every value from position 12
##  onwards shifted by one and `gsmr2_beta` fell back to its default.
##
##  `legacy_published_arguments = TRUE` reproduces the published runs exactly
##  by passing the intended values to the arguments they actually reached.
##  Set it to FALSE to use the intended parameter values - this WILL change
##  the GSMR results and requires re-running steps 09-11.
GSMR_PARAMS <- list(

  # --- arguments that were passed correctly in the published run -----------
  n_ref                 = 503,   # reference sample size for the LD matrix (1kg EUR)
  gwas_thresh           = 5e-8,  # instrument selection threshold
  heidi_outlier_flag    = TRUE,
  mhc_chr               = "6",   # MHC region excluded from instruments
  mhc_start             = 28477797,
  mhc_end               = 33448354,

  # --- intended values (used when legacy_published_arguments = FALSE) ------
  single_snp_heidi_thresh = 0.01,
  multi_snps_heidi_thresh = 0.01,
  nsnps_thresh            = 5,     # minimum number of instruments per CpG
  ld_r2_thresh            = 0.05,  # prune SNP pairs in high LD before GSMR
  ld_fdr_thresh           = 0.05,
  gsmr2_beta              = 1,     # 1 = gsmr2 HEIDI-outlier implementation

  # --- what the published run actually delivered --------------------------
  legacy_published_arguments = TRUE,
  legacy_effective_values = list(
    single_snp_heidi_thresh = 0.01,  # script's "multi_snps_heidi_thresh"
    multi_snps_heidi_thresh = 5,     # script's "nsnps_thresh"  -> HEIDI-outlier never fires
    nsnps_thresh            = 0.05,  # script's "ld_r2_thresh" -> no minimum instrument count
    ld_r2_thresh            = 0.05,  # script's "ld_fdr_thresh"
    ld_fdr_thresh           = 1,     # script's "gsmr2_beta"    -> LD-FDR pruning disabled
    gsmr2_beta              = 0      # gsmr2 default            -> legacy HEIDI-outlier method
  )
)

## ---- 6.4 Colocalisation (steps 13-14) -------------------------------------
COLOC_PARAMS <- list(
  window_kb = 1000,              # cis window (+/- 1 Mb) around the CpG
  type1     = "quant",           # meQTL is a quantitative trait
  type2     = "cc"               # breast cancer GWAS is case-control
)

## ---- 6.5 Mediation analysis (step 20) -------------------------------------
MEDIATION_PARAMS <- list(
  ci_alpha = 0.05,               # 95% confidence interval
  ci_type  = "prodclin",         # RMediation::medci distribution-of-product method
  n_boot   = NULL                # Monte-Carlo bootstrap not used; product-of-coefficients CI
)

## ---- 6.6 Single-cell RNA-seq: quality control and clustering --------------
##  These values are the ones reported in the manuscript Methods.  They are
##  collected here so that the exact thresholds used can be re-used verbatim.
##  See docs/02_scRNA_seq_QC_parameters.md for the rationale of each choice.
SC_QC_PARAMS <- list(

  # ---- Cell / gene calling -------------------------------------------------
  min_cells_per_gene  = 3,       # a gene must be detected in >= 3 cells to be kept
  min_features_per_cell = 300,   # a cell must express >= 300 genes

  # ---- Empty-droplet and doublet-associated filtering ----------------------
  min_features        = 300,     # lower bound on nFeature_RNA (unique genes)
  max_percent_mt      = 20,      # upper bound on mitochondrial read fraction (%)
  mito_pattern        = "^MT-",  # human mitochondrial gene symbol prefix

  # ---- Normalisation -------------------------------------------------------
  norm_method         = "LogNormalize",
  scale_factor        = 10000,
  variable_nfeatures  = 2000,
  selection_method    = "vst",

  # ---- Dimensionality reduction and batch correction -----------------------
  n_pcs               = 50,      # PCs computed by RunPCA
  pcs_used            = 1:10,    # PCs carried forward into the graph / UMAP

  ## Harmony integration variable.
  ##   "group"      -> 5 receptor-status batches (Normal / ER+ / HER2+ / PR+ /
  ##                   TNBC).  This reproduces the published run, in which the
  ##                   Seurat `orig.ident` field held the receptor-status
  ##                   label.  See docs/04_known_issues.md.
  ##   "orig.ident" -> per-donor (GSM) integration.  This is the recommended
  ##                   choice for any re-analysis, because correcting on the
  ##                   group label only removes between-group differences and
  ##                   leaves within-group donor effects in the embedding.
  harmony_group_var   = "group",
  harmony_theta       = 2,
  harmony_max_iter    = 10,
  harmony_nclust      = NULL,    # NULL -> harmony's own default

  # ---- Graph-based clustering ---------------------------------------------
  k_param             = 20,      # FindNeighbors k.param
  cluster_resolution  = 1,       # Leiden/Louvain resolution
  umap_min_dist       = 0.3,
  umap_n_neighbors    = 30,

  # ---- Marker discovery ----------------------------------------------------
  marker_only_pos     = TRUE,
  marker_min_pct      = 0.25,
  marker_logfc_thresh = 0.25,
  marker_logfc_keep   = 0.5,     # additional filtering of the marker table
  marker_padj_keep    = 0.05,

  # ---- Reproducibility -----------------------------------------------------
  seed                = 12345
)

## ---- 6.7 Single-cell group label used for the case/control contrast -------
SC_GROUP_LABELS <- list(
  cancer = c("ER+", "HER+", "PR+", "TNBC"),
  normal = "Normal"
)

## ---- 6.8 Plotting theme ---------------------------------------------------
PALETTE <- c("#1F78B4", "#9E0142", "#f8c120", "#4197d8", "#E16E6D",
             "#cab2d6", "#62B197", "#B3B3B3", "#6a3d9a", "#D0E7ED")
THEME_BASE_SIZE <- 10

# ---------------------------------------------------------------------------
# 7. Sanity check helper
# ---------------------------------------------------------------------------
#  Call require_files(PATHS$x, PATHS$y) at the top of a script that depends on
#  large inputs, so that a missing file fails with a readable message instead
#  of a cryptic error deep inside a loop.
require_files <- function(...) {
  p <- unlist(list(...), use.names = FALSE)
  missing <- p[!file.exists(p)]
  if (length(missing)) {
    stop("Missing input file(s):\n  ", paste(missing, collapse = "\n  "),
         "\nSee docs/03_data_sources.md and config/paths_local.R.example.",
         call. = FALSE)
  }
  invisible(TRUE)
}

message("[config] repository root : ", REPO_ROOT)
message("[config] output directory: ", OUTPUT_DIR)

# 4DMP-MR-BC

**Causal and druggable discovery in breast cancer from blood DNA methylation:
an integrated cis-meQTL Mendelian randomisation, colocalisation, proteomic
mediation and single-cell analysis.**

This repository contains the complete analysis code for the study, organised as
a six-stage pipeline from differentially methylated positions through causal
inference to the single-cell context of the resulting protein mediators.

```
GSE104942 blood 450K          GoDMC blood cis-meQTL
        │                             │
        └── 19,169 DMPs ──────────────┘
                    │
     TSMR vs FinnGen R11 / BCAC 2020 / GCST90011804
                    │
     intersection + GSMR + meta-analysis
                    │
              3 causal CpGs
                    │
     ┌──────────────┼──────────────────┐
 coloc.abf     FinnGen PheWAS     UKB-PPP mediation
                                      │
                          15 CpG-protein pathways
                                      │
                        GSE161529 single-cell atlas
                            (271,532 cells)
```

---

## Quick start

```bash
git clone <REPOSITORY URL>
cd 4DMP-MR-BC
cp config/paths_local.R.example config/paths_local.R
$EDITOR config/paths_local.R        # point at your local copies of the data
```

Keep `config/paths_local.R` ASCII-only where possible; if a data directory
genuinely contains non-ASCII characters, save the file as UTF-8 without a BOM
(`config.R` loads it with a UTF-8-aware reader that tolerates a non-UTF-8 R
locale, but a file re-saved in a legacy codepage cannot be recovered).

Then run the stages in order from the repository root (or open
`4DMP-MR-BC.Rproj` in RStudio):

```r
source("config/config.R")                       # verify the paths resolve
Rscript("01_dmp/01_dmp_discovery_GSE104942.R")  # or: Rscript 01_dmp/01_...R
```

All inputs, outputs and analysis parameters are declared in
`config/config.R`; no absolute path appears anywhere else in the code.

---

## Repository layout

```
config/
    config.R                     central paths, parameters, helpers
    paths_local.R.example        template for machine-specific data locations
helpers/
    mr_functions.R               instrument building, TSMR, robustness filters
    gsmr_functions.R             GSMR scan with explicit argument mapping
    coloc_functions.R            coloc.abf wrapper, H0-H4 posteriors
    enrichment_plot.R            GO / KEGG rounded-bar plots
data_preparation/                00a-00c   format outcome and reference GWAS
01_dmp/                          01-02     differential methylation discovery
02_meqtl_mr/                     03-13     cis-meQTL MR, GSMR, meta-analysis
03_colocalization/               14-16     coloc.abf and summary figure
04_phewas/                       17-18     phenome-wide MR screen
05_protein_mediation/            19-24     CpG -> protein -> BC mediation
06_single_cell/                  25-27     scRNA-seq QC, clustering, annotation
docs/
    01_workflow_overview.md      stage-by-stage description
    02_scRNA_seq_QC_parameters.md  ** every scRNA-seq QC parameter **
    03_data_sources.md           accessions and download instructions
    04_known_issues.md           ** read before publishing the code **
```

---

## Analysis parameters at a glance

| Stage | Key thresholds |
|---|---|
| DMP discovery | k-NN imputation k = 10; mean beta > 0.005; `betaqn` normalisation; `dmpFinder` FDR < 0.05 |
| Instrument selection | P < 5e-8; LD r² < 0.001 within 10,000 kb (1000G EUR); F > 10; palindromic SNPs resolved with `action = 2` |
| Instrument exclusion | SNPs associated with the outcome at P < 1e-5 removed |
| Primary MR | IVW (≥ 2 instruments) / Wald ratio (1 instrument); Cochran's Q P > 0.05; MR-Egger intercept P > 0.05; correct Steiger direction |
| Multiple testing | Benjamini-Hochberg FDR and Bonferroni on all CpG-outcome tests |
| Meta-analysis | random-effects `metagen` across 4 method × outcome strata |
| Colocalisation | ± 1 Mb cis window, `coloc.abf`, `type = "quant"` vs `"cc"` |
| Mediation | product of coefficients, 95 % CI by `RMediation::medci(type = "prodclin")` |
| **scRNA-seq QC** | **`min.cells = 3`; `min.features = 300`; `nFeature_RNA > 300`; `percent.mt < 20`; LogNormalize (10,000); 2,000 HVGs (vst); 50 PCs; Harmony `theta = 2`, `max_iter = 10`; 10 components; Louvain resolution 1; UMAP `n.neighbors = 30`, `min.dist = 0.3`; seed 12345** |

The single-cell parameters are documented in full — with rationale, the
observed QC distributions and a ready-to-paste Methods paragraph — in
[`docs/02_scRNA_seq_QC_parameters.md`](docs/02_scRNA_seq_QC_parameters.md).

---

## Environment

Developed and run under **R 4.2.3 (2023-03-15, ucrt)** on Windows 10/11 x64.

Core packages:

| Package | Version | Used in |
|---|---|---|
| Seurat | 4.4.0 | steps 25-27 |
| harmony | 1.2.1 | step 25 |
| TwoSampleMR | 0.6.17 | steps 05-07, 17, 19-20 |
| gsmr2 | 1.1.1 | steps 09-11 |
| coloc | 5.2.3 | steps 14-15 |
| GagnonMR | 0.0.0.9000 | steps 05-07, 17, 19-20 |
| ieugwasr | — | LD clumping |
| minfi, wateRmelon, impute | 1.46.0 / 2.6.0 / — | step 01 |
| clusterProfiler, org.Hs.eg.db | 4.6.2 / — | steps 19-20 |
| CMplot, meta, RMediation, ggalluvial, plot1cell, scRNAtoolVis | — | figures |

`sessionInfo()` is captured to `output/qc/25_session_info.txt` by step 25.

---

## Outputs

Everything the pipeline produces is written under `output/`:

```
output/objects/     intermediate and final R objects
output/tables/      every result table as CSV (plus an RDS copy)
output/figures/     all figures as PDF (vector, editable)
output/qc/          quality-control reports and session information
```

---

## Before you publish or reuse this code

`docs/04_known_issues.md` documents two substantive issues that must be
reviewed before the code is made public or before the manuscript's claims
are finalised:

1. **GSMR argument mis-alignment.** The published GSMR runs passed a
   parameter list written for `gsmr` v1 to `gsmr2` 1.1.1, which shifted five
   parameters by one position and effectively disabled HEIDI-outlier
   filtering and LD-FDR pruning. The repository reproduces the published
   numbers by default (`GSMR_PARAMS$legacy_published_arguments = TRUE`) and
   emits a warning; set it to `FALSE` to run the intended analysis.
2. **Single-cell integration** was performed on the receptor-status group
   label (5 batches) rather than the individual donor (47 batches).

---

## Data availability

No data are redistributed here. See
[`docs/03_data_sources.md`](docs/03_data_sources.md) for every accession,
download location and licence.

## Citation

If you use this code, please cite the associated paper (see
`CITATION.cff`) and the underlying resources listed in
`docs/03_data_sources.md`.

## License

MIT — see [LICENSE](LICENSE).

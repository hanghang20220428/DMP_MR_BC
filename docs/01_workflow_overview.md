# Analysis workflow

The pipeline proceeds through six stages. Each stage is a directory of
numbered R scripts; scripts write only to `output/` and read their inputs
from `output/objects/` or from an external data root configured in
`config/paths_local.R`.

```
   ┌─────────────────────────────────────────────────────────────────┐
   │ 0. Data preparation                                             │
   │    outcome and reference GWAS formatted for downstream tools    │
   └─────────────────────────────────────────────────────────────────┘
                              │
   ┌─────────────────────────────────────────────────────────────────┐
   │ 1. DMP discovery  (GSE104942, blood 450K)                       │
   │    cases vs controls -> differentially methylated positions     │
   └─────────────────────────────────────────────────────────────────┘
                              │  19,169 DMPs (FDR < 0.05)
   ┌─────────────────────────────────────────────────────────────────┐
   │ 2. cis-meQTL MR   (GoDMC blood meQTL x 3 breast-cancer GWAS)    │
   │    instruments -> TSMR (FinnGen, BCAC, GCST90011804)            │
   │    -> intersection -> GSMR sensitivity -> forest -> meta        │
   └─────────────────────────────────────────────────────────────────┘
                              │  3 causal CpGs
   ┌─────────────────────────────────────────────────────────────────┐
   │ 3. Colocalisation   coloc.abf in the +/- 1 Mb cis window        │
   └─────────────────────────────────────────────────────────────────┘
   ┌─────────────────────────────────────────────────────────────────┐
   │ 4. Phenome-wide screen   causal CpGs x all FinnGen R11 endpoints│
   └─────────────────────────────────────────────────────────────────┘
   ┌─────────────────────────────────────────────────────────────────┐
   │ 5. Protein mediation   CpG -> UKB-PPP protein -> breast cancer  │
   │    product-of-coefficients mediation                            │
   └─────────────────────────────────────────────────────────────────┘
                              │  15 significant CpG-protein pathways
   ┌─────────────────────────────────────────────────────────────────┐
   │ 6. Single-cell context  (GSE161529, 271,532 cells)              │
   │    QC -> clustering -> annotation -> mediator expression        │
   └─────────────────────────────────────────────────────────────────┘
```

---

## Stage 0 — Data preparation (`data_preparation/`)

| Step | Script | Produces |
|---|---|---|
| 00a | `00a_finngen_r11_bc_to_coloc.R` | coloc-ready FinnGen R11 breast-cancer table |
| 00b | `00b_pancancer17_to_ldsc.R` | LDSC-format pan-cancer GWAS (optional branch) |
| 00c | `00c_pancancer17_to_mrdata.R` | GWAS Catalog summary statistics as TwoSampleMR objects |

## Stage 1 — Differential methylation (`01_dmp/`)

| Step | Script | Key parameters |
|---|---|---|
| 01 | `01_dmp_discovery_GSE104942.R` | k-NN imputation (k = 10); probes with mean beta < 0.005 dropped; `betaqn` quantile normalisation; `dmpFinder` on M-values, FDR < 0.05 |
| 02 | `02_dmp_heatmap.R` | row-scaled beta heat map, samples ordered by group |

Result: **19,169 DMPs** at FDR < 0.05, which become the exposure candidates.

## Stage 2 — cis-meQTL Mendelian randomisation (`02_meqtl_mr/`)

| Step | Script | Notes |
|---|---|---|
| 03 | `03_godmc_meqtl_cleaning.R` | GoDMC meQTL cleaned; `chr:pos` mapped to rsIDs |
| 04 | `04_build_godmc_instruments.R` | one instrument set per CpG: P < 5e-8, LD r² < 0.001, 10,000 kb, F > 10 |
| 05 | `05_tsmr_finngen.R` | IVW / Wald ratio vs FinnGen R11 (20,586 cases) |
| 06 | `06_tsmr_bcac.R` | vs BCAC 2020 (N = 247,173) |
| 07 | `07_tsmr_gcst90011804.R` | vs GCST90011804 (exploratory third outcome) |
| 08 | `08_intersect_union_cpg.R` | CpG set A (intersection), B (union), C (three-way) |
| 09-11 | `09_...`, `10_...`, `11_...` | GSMR + HEIDI-outlier sensitivity — see `docs/04_known_issues.md` §1 |
| 12 | `12_forest_plot.R` | TSMR vs GSMR forest plot |
| 13 | `13_meta_analysis.R` | random-effects meta-analysis across the four method x outcome strata |

Result: **566** CpGs significant against FinnGen, **706** against BCAC
(heterogeneity / pleiotropy / directionality filtered), converging on a
**3-CpG** causal set after meta-analysis.

## Stage 3 — Colocalisation (`03_colocalization/`)

| Step | Script | Notes |
|---|---|---|
| 14 | `14_coloc_finngen.R` | `coloc.abf` in the ± 1 Mb cis window, `type = "quant"` vs `"cc"` |
| 15 | `15_coloc_bcac.R` | as above, BCAC outcome |
| 16 | `16_coloc_dotplot.R` | bubble plot of the lead-SNP posterior per CpG (`top_snp_PP.H4`) |

## Stage 4 — Phenome-wide screen (`04_phewas/`)

| Step | Script | Notes |
|---|---|---|
| 17 | `17_phewas_mr.R` | each causal CpG × every FinnGen R11 endpoint (5,587 tests) |
| 18 | `18_phewas_circos_plot.R` | circos heat map of the top 20 endpoints per CpG |

## Stage 5 — Protein mediation (`05_protein_mediation/`)

| Step | Script | Notes |
|---|---|---|
| 19 | `19_cpg_to_protein_tsmr.R` | CpG → UKB-PPP protein (2,940 aptamers) |
| 20 | `20_protein_to_bc_tsmr.R` | UKB-PPP protein → BCAC breast cancer |
| 21 | `21_mediation_analysis.R` | product of coefficients, 95 % CI by `RMediation::medci` |
| 22 | `22_phenogram_input.R` | PhenoGram input tables |
| 23 | `23_cpg_methylation_boxplot.R` | CpG beta values by group in the discovery cohort |
| 24 | `24_alluvial_plot.R` | SNP → CpG → protein → outcome alluvial diagram |

Result: **15** significant mediated CpG-protein pairs across the three causal
CpGs (13 for `cg15146122`, 1 for `cg24704912`, 1 for `cg01529201`).

## Stage 6 — Single-cell context (`06_single_cell/`)

| Step | Script | Notes |
|---|---|---|
| 25 | `25_scrna_qc_and_clustering.R` | QC, integration, clustering — **all thresholds in `docs/02_scRNA_seq_QC_parameters.md`** |
| 26 | `26_scrna_cell_type_annotation.R` | markers and manual annotation into 7 cell types |
| 27 | `27_scrna_mediator_expression.R` | expression of the mediator proteins, cancer vs normal |

Result: **271,532 cells**, 28,273 genes, 35 clusters, 7 cell types.

---

## Running the pipeline

```bash
git clone <REPOSITORY URL>
cd 4DMP-MR-BC
cp config/paths_local.R.example config/paths_local.R
$EDITOR config/paths_local.R          # point at your local data
```

Then, with the repository root as the working directory:

```r
source("config/config.R")             # sanity check the paths
source("01_dmp/01_dmp_discovery_GSE104942.R")
```

or, non-interactively:

```bash
Rscript 01_dmp/01_dmp_discovery_GSE104942.R
```

Steps must be run in numeric order. Steps 05-07 and 17 are the slow ones
(tens of thousands of MR fits); the rest complete in minutes each.

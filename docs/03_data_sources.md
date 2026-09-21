# Data sources

No third-party data are redistributed with this repository. Every input is
publicly available; download instructions and the corresponding entry in
`config/config.R` are listed below.

---

## Discovery cohort

### GSE104942 — blood DNA methylation, breast cancer cases vs controls
* **Type**  Illumina HumanMethylation450 BeadChip, whole blood
* **Access** GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104942>
* **File**  `GSE104942_series_matrix.txt.gz`
* **Config** `gse104942_matrix`

---

## Exposure: blood cis-meQTL

### GoDMC — Genetics of DNA Methylation Consortium
* **Type**  meta-analysis of cis-meQTLs in whole blood, GRCh37
* **Access** <https://mqtldb.godmc.org.uk/>
* **Files** `assoc_meta_all.csv.gz` (SNP–CpG pairs), `snps.csv.gz`
  (`chr:pos` → rsID lookup)
* **Config** `godmc_assoc`, `godmc_snps`, `godmc_exposure_dir`

---

## Outcome GWAS

### FinnGen release 11 — C3_BREAST_EXALLC
* **Type**  breast cancer, 20,586 cases / 201,494 controls, GRCh38
* **Access** <https://www.finngen.fi/en/access_results>
* **Files** `finngen_R11_C3_BREAST_EXALLC.gz`, `finngen_R11_manifest.csv`
* **Config** `finngen_bc_raw`, `finngen_manifest`
* The full R11 endpoint collection is additionally used for the phenome-wide
  screen (`finngen_all_outcomes`).

### BCAC 2020 — breast cancer susceptibility meta-analysis
* **Type**  iCOGS + OncoArray + UK Biobank, N = 247,173, GRCh37
* **Access** <http://bcac.ccge.medschl.cam.ac.uk/>
* **Files** `icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt`
* **Config** `bcac_sus_raw`

### GCST90011804 — breast cancer (BCAC / BRIDGES / CIMBA)
* **Type**  breast cancer, GRCh37, published without beta / SE / EAF
* **Access** GWAS Catalog: <https://www.ebi.ac.uk/gwas/studies/GCST90011804>
* **Reference** Rashkin SR, et al. Pan-cancer study detects genetic risk
  variants and shared genetic basis in two large cohorts. *Nat Commun*
  2020;11:4423 (PMID 32887889)
* **Config** `gcst90011804_raw`

---

## Mediators: plasma proteomics

### UK Biobank Pharma Proteomics Project (UKB-PPP)
* **Type**  SomaScan plasma proteomics, 2,940 aptamers, cis-pQTL instruments
* **Access** <https://www.synapse.org/#!Synapse:syn51364943>
* **Config** `ukbppp_outcome_dir`, `ukbppp_exposure_dir`

---

## Single-cell atlas

### GSE161529 — single-cell transcriptomic atlas of the human breast
* **Type**  10x Genomics Chromium 3' count matrices
* **Access** GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529>
* **Reference** Pal B, Chen Y, Vaillant F, et al. A single-cell RNA expression
  atlas of normal, preneoplastic and tumorigenic states in the human breast.
  *EMBO J* 2021;40:e107333. PMID 33950524
* **Layout expected by the pipeline**: `<scrna_10x_root>/<GSM>/` containing
  `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`
* **Config** `scrna_10x_root`

The 47 GSM accessions used here are listed in `GSE161529_SAMPLES` in
`config/config.R`.

---

## Reference files

| Resource | Purpose | Config entry |
|---|---|---|
| 1000 Genomes EUR (PLINK 1 binary) | LD clumping, LD correlation matrices | `plink_bfile`, `plink_bin` |
| GPL13534 manifest | 450K probe coordinates (GRCh37) | `anno_450k` |
| GPL21145 manifest | EPIC probe coordinates (GRCh37) | `anno_epic` |
| GRCh37.87 gene locations | gene start positions | `generef_grch37` |
| GRCh38.100 gene locations | gene start positions | `generef_grch38` |
| ALFA European allele frequencies | back-filling missing EAFs | `snp_eaf_eur` |
| EPIC meQTL reference | CpG → chromosome/position lookup | `epic_meqtl` |

The 1000 Genomes reference panel can be obtained from
<https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip> and converted to
PLINK 1 binary format with `plink --make-bed`; `PLINK_BIN` must point at a
PLINK 1.9 executable.

---

## Tools that are not R packages

| Tool | Version used | Purpose | Config entry |
|---|---|---|---|
| PLINK 1.9 | 20230116 | LD clumping and LD matrices | `plink_bin` |
| gsmr2 | R package 1.1.1 | GSMR + HEIDI-outlier | — |
| PhenoGram | web service | chromosome ideograms for step 22 | manual upload |

PhenoGram inputs written by step 22 are tab-separated with the columns
`cpg_pro`, `chr`, `pos`, `phenotype` and are uploaded at
<https://visualization.ritchielab.org/phenograms/plot>.

# Single-cell RNA-seq quality control and processing parameters

This document answers the editor's and reviewers' request for the
"specific threshold choices, filtering criteria, and software parameters"
used for the scRNA-seq analysis. Every value below is set from the named
parameter block `SC_QC_PARAMS` in `config/config.R` and is applied by
`06_single_cell/25_scrna_qc_and_clustering.R`; nothing is left at an implicit
default.

---

## 1. Dataset

| Item | Value |
|---|---|
| Accession | **GSE161529** (GEO) |
| Reference | Pal B, Chen Y, Vaillant F, et al. A single-cell RNA expression atlas of normal, preneoplastic and tumorigenic states in the human breast. *EMBO J* 2021;40:e107333. |
| Platform | 10x Genomics Chromium 3' (Cell Ranger count matrices) |
| Samples retained | 47 (Normal 13, ER+ 19, HER2+ 6, PR+ 1, TNBC 8) |
| Reference genome | GRCh38 |

The 47 GSM accessions retained here are listed explicitly in
`GSE161529_SAMPLES` in `config/config.R`.

---

## 2. Cell calling, filtering and the resulting object

### 2.1 Cell calling — applied per sample by `CreateSeuratObject()`

| Parameter | Value | Meaning |
|---|---|---|
| `min.cells` | **3** | a gene is retained only if detected in at least 3 cells |
| `min.features` | **300** | a cell is retained only if at least 300 genes are detected |

Applied separately to each of the 47 samples **before** merging, so that a
gene or a cell is never rescued by depth contributed by another library.

### 2.2 Cell filtering — applied after merging by `subset()`

| Parameter | Value | Meaning |
|---|---|---|
| `nFeature_RNA` | **> 300** | minimum number of detected genes per cell |
| `percent.mt` | **< 20** | maximum mitochondrial read fraction, computed with `PercentageFeatureSet(pattern = "^MT-")` |

No lower bound was placed on `nCount_RNA` and no upper bound on
`nFeature_RNA`.

> **Note.** The `nFeature_RNA > 300` criterion is largely redundant with the
> `min.features = 300` filter already applied at object creation; it is
> retained because it was part of the pre-specified filtering step.

### 2.3 Resulting object

| Item | Value |
|---|---|
| Cells retained | **271,532** |
| Genes retained | **28,273** |
| Cells — Normal | 58,372 |
| Cells — ER+ | 110,347 |
| Cells — HER2+ | 42,411 |
| Cells — PR+ | 3,176 |
| Cells — TNBC | 57,226 |

Observed distributions **in the retained object** (they confirm that the
stated thresholds are the binding ones — nothing outside them survives):

| Metric | Minimum | Median | 99th pct | Maximum |
|---|---|---|---|---|
| `nFeature_RNA` (genes/cell) | 301 | 1,130 | 5,600 | 10,660 |
| `nCount_RNA` (UMI/cell) | 374 | 2,954 | 38,425 | 232,130 |
| `percent.mt` (%) | 0.00 | 5.46 | 18.91 | 20.00 |

Per-group medians:

| Group | Cells | Median genes/cell | Median UMIs/cell | Median `percent.mt` |
|---|---|---|---|---|
| Normal | 58,372 | 1,549 | 4,557 | 3.44 % |
| ER+ | 110,347 | 850 | 2,000 | 6.15 % |
| HER2+ | 42,411 | 874 | 2,467 | 7.02 % |
| PR+ | 3,176 | 931 | 2,294 | 3.93 % |
| TNBC | 57,226 | 1,419 | 3,992 | 5.74 % |

When the pipeline is re-run from the raw 10X matrices the number of cells
removed by each criterion is written to
`output/qc/25_qc_retention.csv`, and the per-sample metrics to
`output/qc/25_qc_per_sample.csv`.

---

## 3. Normalisation and feature selection

| Step | Function | Parameters |
|---|---|---|
| Normalisation | `NormalizeData` | `normalization.method = "LogNormalize"`, `scale.factor = 10000` |
| Highly variable genes | `FindVariableFeatures` | `selection.method = "vst"`, `nfeatures = 2000` |
| Scaling | `ScaleData` | all genes; no covariates regressed out |

---

## 4. Dimensionality reduction, integration and clustering

| Step | Function | Parameters |
|---|---|---|
| PCA | `RunPCA` | `npcs = 50` |
| Batch integration | `RunHarmony` | `group.by.vars = "group"` (receptor-status label), `theta = 2`, `max_iter = 10` |
| Neighbour graph | `FindNeighbors` | `reduction = "harmony"`, `dims = 1:10`, `k.param = 20` |
| Clustering | `FindClusters` | `resolution = 1` (Louvain) → **35 clusters** |
| Embedding | `RunUMAP` | `reduction = "harmony"`, `dims = 1:10`, `n.neighbors = 30`, `min.dist = 0.3` |
| Random seed | `set.seed(12345)` | `FindClusters` and `RunUMAP` additionally use their own deterministic internal seeds |

> **Important – batch variable.** In the released object the integration
> variable encodes the **receptor-status group** (5 batches: Normal, ER+,
> HER2+, PR+, TNBC). `SC_QC_PARAMS$harmony_group_var` should be changed to
> `"orig.ident"` to integrate on the individual donor (47 batches), which is
> the recommended choice for any re-analysis. See
> `04_known_issues.md` §2.

---

## 5. Cluster marker discovery and annotation

| Step | Function | Parameters |
|---|---|---|
| Marker detection | `FindAllMarkers` | `only.pos = TRUE`, `min.pct = 0.25`, `logfc.threshold = 0.25`, Wilcoxon rank-sum test (Seurat default) |
| Marker significance | — | `|avg_log2FC| > 0.5` and Bonferroni-adjusted `p_val_adj < 0.05` |
| Annotation | manual | canonical lineage markers (`EPCAM`, `CD3E`, `CD4`, `CD8A`, `MS4A1`, `CD79A`, `LYZ`, `CD68`, `C1QA`, `FCER1A`, `CD1C`, `CLEC9A`, `PECAM1`, `VWF`, `COL1A1`, `LUM`, `PDGFRA`, `ACTA2`, `MYH11`) |

Seven cell types were identified:

| Cell type | Cells | % of total |
|---|---|---|
| Epithelial cells | 143,131 | 52.7 % |
| T cells | 46,615 | 17.2 % |
| Fibroblasts | 41,913 | 15.4 % |
| Macrophages | 21,529 | 7.9 % |
| Endothelial cells | 11,147 | 4.1 % |
| B cells | 6,411 | 2.4 % |
| Dendritic cells | 786 | 0.3 % |

The full cluster-to-cell-type mapping is in
`06_single_cell/26_scrna_cell_type_annotation.R` (`cluster_map`).

---

## 6. Software versions

The analysis was performed in **R 4.2.3 (2023-03-15, ucrt)** on Windows.

| Package | Version |
|---|---|
| Seurat | 4.4.0 |
| harmony | 1.2.1 |
| plot1cell | 0.1.0 |
| scRNAtoolVis | 0.0.7 |
| ggplot2 | 3.5.2 |
| dplyr | (see `sessionInfo`) |

`sessionInfo()` is written to `output/qc/25_session_info.txt` every time step
25 is executed, so the exact environment is captured alongside the results.

---

## 7. Ready-to-paste Methods paragraph (English)

> **Single-cell RNA-seq processing and quality control.**
> Publicly available 10x Genomics Chromium 3' count matrices for 47 samples
> of normal breast tissue and breast tumours were obtained from GEO accession
> GSE161529. Count matrices were read with `Read10X()` and Seurat objects
> were created per sample with `CreateSeuratObject(min.cells = 3,
> min.features = 300)`, so that genes detected in fewer than three cells and
> cells expressing fewer than 300 genes were excluded before merging (Seurat
> 4.4.0). Samples were then merged, and the mitochondrial read fraction was
> computed for each cell with `PercentageFeatureSet(pattern = "^MT-")`.
> Cells with fewer than 300 detected genes or with a mitochondrial read
> fraction of 20 % or more were removed (`nFeature_RNA > 300 & percent.mt <
> 20`); 271,532 cells and 28,273 genes were retained (median 1,130 genes and
> 2,954 UMIs per cell; median mitochondrial fraction 5.5 %). Counts were
> normalised with `NormalizeData()` (LogNormalize, scale factor 10,000) and
> the 2,000 most variable genes were selected with
> `FindVariableFeatures(selection.method = "vst")`. After `ScaleData()` and
> `RunPCA(npcs = 50)`, batch effects were corrected with Harmony
> (`RunHarmony(group.by.vars = "group", theta = 2, max_iter = 10)`), where
> "group" denotes the receptor-status annotation (Normal, ER+, HER2+, PR+,
> TNBC). The first 10 Harmony components were used for graph construction
> (`FindNeighbors(k.param = 20)`), clustering (`FindClusters(resolution = 1)`)
> and UMAP embedding (`RunUMAP(n.neighbors = 30, min.dist = 0.3)`), yielding
> 35 clusters. Cluster markers were identified with `FindAllMarkers(only.pos =
> TRUE, min.pct = 0.25, logfc.threshold = 0.25)` and retained when
> |avg_log2FC| > 0.5 and Bonferroni-adjusted P < 0.05. Clusters were
> annotated manually against canonical lineage markers, yielding seven cell
> types (epithelial cells, T cells, fibroblasts, macrophages, endothelial
> cells, B cells and dendritic cells). All analyses used `set.seed(12345)`.
> The complete analysis code, including the exact parameter values listed
> above, is available at <REPOSITORY URL>.

---

## 8. Response-to-reviewer snippet

> **Reviewer comment.** *Please provide any details (parameters) of scRNA-seq
> QC or code for this analysis.*
>
> **Response.** Thank you for asking us to make this explicit. The complete
> single-cell analysis code, with every quality-control threshold and software
> parameter set from a single named parameter block, is now available at
> <REPOSITORY URL> (`06_single_cell/25_scrna_qc_and_clustering.R` and
> `config/config.R`). In brief, genes detected in fewer than 3 cells and cells
> expressing fewer than 300 genes were excluded at object creation
> (`CreateSeuratObject(min.cells = 3, min.features = 300)`); after merging, a
> further filter removed cells with fewer than 300 detected genes or with a
> mitochondrial read fraction ≥ 20 % (`nFeature_RNA > 300 & percent.mt < 20`,
> mitochondrial genes identified by `^MT-`). A total of 271,532 cells and
> 28,273 genes were retained. Normalisation used LogNormalize with a scale
> factor of 10,000, 2,000 highly variable genes were selected by the
> variance-stabilising transformation, and Harmony integration was performed
> with `theta = 2` and `max_iter = 10`. Clustering used the first 10 Harmony
> components with Louvain resolution 1, and the UMAP embedding used
> `n.neighbors = 30` and `min.dist = 0.3`. These values have been added to
> the Methods section (page …, lines …).

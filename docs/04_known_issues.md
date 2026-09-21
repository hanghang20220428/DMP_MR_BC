# Known issues and differences from the originally published run

This file records every place where the released code either
(a) intentionally reproduces a quirk of the originally run scripts, or
(b) corrects a defect. Nothing here is hidden: each item states the impact
and the one-line change that resolves it.

Legend: **[review before publishing]** — a reviewer could reasonably raise
this; the repository currently reproduces the published behaviour so that the
paper remains reproducible.

---

## 1. GSMR: positional argument mis-alignment — **[review before publishing]**

**Where** `helpers/gsmr_functions.R`, used by steps 09, 10 and 11.

**What happened.** The original scripts called `gsmr2::gsmr()` with
positional arguments written against the argument list of **gsmr v1**:

```r
gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp, n_ref,
     heidi_outlier_flag, gwas_thresh,
     multi_snps_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh,
     gsmr2_beta)
```

`gsmr2` 1.1.1 inserts **`single_snp_heidi_thresh`** between `gwas_thresh`
and `multi_snps_heidi_thresh`. Both packages are installed on this machine
and both export a function named `gsmr`; the scripts attach `gsmr2`, so the
`gsmr2` formals are the ones that apply. The first eleven arguments still
land correctly, but from position 12 onwards every value shifted by one and
`gsmr2_beta` fell back to its default.

**Effective values actually used**

| Argument (gsmr2 1.1.1) | Intended | Actually delivered | Consequence |
|---|---|---|---|
| `n_ref` | 503 | 503 | — |
| `heidi_outlier_flag` | TRUE | TRUE | — |
| `gwas_thresh` | 5e-8 | 5e-8 | — |
| `single_snp_heidi_thresh` | 0.01 | 0.01 | — |
| `multi_snps_heidi_thresh` | 0.01 | **5** | **HEIDI-outlier filtering never fires** (a P value can never exceed 5) |
| `nsnps_thresh` | 5 | **0.05** | minimum instrument count effectively removed; GSMR runs even for a single shared SNP |
| `ld_r2_thresh` | 0.05 | 0.05 | — (coincidentally equal) |
| `ld_fdr_thresh` | 0.05 | **1** | LD-FDR pruning of coincidentally correlated SNPs disabled |
| `gsmr2_beta` | 1 | **0** (default) | legacy HEIDI-outlier implementation instead of the gsmr2 method |

**Impact.** Any statement in the manuscript to the effect that GSMR was run
*with HEIDI-outlier filtering* overstates what the code did. The GSMR
estimates themselves are still valid single-instrument / multi-instrument
causal estimates, but the pleiotropy-robustness claim needs re-checking.

**How the repository handles it.** `GSMR_PARAMS$legacy_published_arguments`
defaults to `TRUE`, which reproduces the published numbers exactly and emits
a runtime warning. Set it to `FALSE` in `config/config.R` to run steps 09-11
with the intended parameter values. Doing so **will change the GSMR results**
and requires the corresponding tables and figures to be regenerated.

**Recommended action.** Re-run steps 09-11 with
`legacy_published_arguments = FALSE`, compare the significant CpG sets, and
either (i) update the manuscript, or (ii) delete the HEIDI-outlier claim.

---

## 2. Single-cell integration was performed on the group label, not the donor

**Where** `config/config.R` → `SC_QC_PARAMS$harmony_group_var`.

**What happened.** In the object used for the published analysis, the Seurat
`orig.ident` field held the **receptor-status** label (5 levels: Normal, ER+,
HER2+, PR+, TNBC), and `RunHarmony(group.by.vars = "orig.ident")` was
therefore given 5 batches rather than 47 donors. Correcting on the group
label removes between-group differences but leaves donor-level variation
within each group in the embedding.

**How the repository handles it.** `orig.ident` now holds the per-donor GSM
accession and the receptor-status label lives in a separate `group` column.
`SC_QC_PARAMS$harmony_group_var` selects which of the two Harmony integrates
on:

* `"group"` (default) — reproduces the published integration granularity;
* `"orig.ident"` — per-donor integration, **recommended for any re-analysis**.

**Recommended action.** If the manuscript describes the integration as
removing donor-level batch effects, either re-run step 25 with
`harmony_group_var = "orig.ident"` or reword the Methods to say that
integration was performed across receptor-status groups.

---

## 3. Defensive fixes that do not alter any published result

These changes were made for code quality. For every input that completed
successfully in the original runs they are behaviour-preserving.

| Original construct | Problem | Fix |
|---|---|---|
| `repeat { dat <- try(harmonise_data(...)); if (!inherits(dat, "try-error")) break }` | if `harmonise_data()` fails persistently the loop never terminates | single `try()` plus an explicit failure check |
| `if (is.null(exp)) next` on a data frame | a zero-row data frame is not `NULL`, so the guard never triggers | `if (nrow(exp) == 0) next` |
| `if (is.null(out)) next` on a data frame | same as above | `if (nrow(out) == 0) next` |
| `exp_data %>% filter(cpg == i)` inside a loop | relies on tidy-evaluation of the loop variable | base-subset on an explicit column |
| `setwd()` calls interleaved with analysis | working directory changes mid-script; not portable | all paths resolved from `config/config.R` |
| Absolute paths (`K:/`, `M:/`, `H:/`, `D:/`) | not portable, and the manuscript's code-availability statement would be unverifiable | `config/paths_local.R` |

Two `coloc.abf()` call-site details also belong here, because both fail
**silently** inside the `try()` that wraps the call — the step produces an
empty result table rather than an error:

* **`silent` is not an argument of `coloc.abf()`.** coloc ≥ 5 declares only
  `dataset1, dataset2, MAF, p1, p2, p12` and has no `...`, so passing `silent`
  raises `unused argument`. The released helper no longer passes it, and
  progress output is not suppressed.
* **`dataset2$s` must be a single number.** The coloc-ready outcome tables
  store `s = ncase / samplesize` per row, so the vector used to be handed
  straight to `coloc.abf()`. `check_dataset()` tests
  `!is.numeric(d$s) || d$s <= 0 || d$s >= 1`, and a length > 1 vector only
  survived on R 4.2.3, where `||` silently used the first element; R ≥ 4.3
  raises `'length = 4833' in coercion to 'logical(1)'`. The helper now passes
  `unique(gw$s)` as a scalar. Because `s` is constant within an outcome study
  (verified: 1 unique value across all 20,437,548 FinnGen rows), the two forms
  give **bit-identical** posteriors — max |difference| = 0 over the entire
  per-SNP vector for all three causal CpGs.

---

## 4. Discovery cohort 1 (`01_dmp/01_dmp_discovery_GSE104942.R`)

The original script referred to a `group` object that it never created
(`dmpFinder(..., pheno = group$Group)`, `densityBeanPlot(..., sampGroups =
group$Group)`). In that form the script could not run from a clean session.
The released version uses the `clinical` data frame that the script does
build, which is the evident intent. **Verify this against the original
analysis output before publishing**, since it cannot be excluded that the
published DMP list was produced by a slightly different script revision.

---

## 5. Dataset naming: GSE161529 vs GSE161259

The single-cell scripts and saved objects use both accessions
interchangeably (`GSE161529Normal13.rdata`, `GSE161259_umap.rdata`,
`0BC_GSE161259_Celltype.rdata`). The GSM accessions hard-coded in the scripts
(`GSM4909253`–`GSM4909320`) belong to **GSE161529** (Pal et al., *EMBO J*
2021). GSE161259 is a different study.

The released code uses **GSE161529** throughout. Please confirm the accession
quoted in the manuscript's Data availability statement and Methods.

---

## 6. `circos.par` assignment style

`04_phewas/18_phewas_circos_plot.R` uses the documented
`circos.par(gap.degree = ..., start.degree = ...)` form. The
`circos.par$gap.degree <- 60` form used in the original scripts was verified
to behave identically under circlize 0.4.16, so this is cosmetic only.

---

## 7. Repository infrastructure notes (no effect on analysis results)

These concern the scaffolding that makes the pipeline portable. They are
recorded because each one is an easy trap to re-introduce.

**Reading Illumina probe manifests.** The GEO platform files `GPL13534`
(450K, 485,577 probes) and `GPL21145` (EPIC, 868,564 probes) both open with a
block of comment lines written as `#COLUMN = ...` (e.g. `#ID = IlmnID`)
before the real tab-delimited header, and the two releases do not use exactly
the same column set. `read.delim()` mis-handles these files;
`data.table::fread()` skips the comment block and lands on the true header.
`helpers/mr_functions.R::read_probe_manifest()` wraps `fread()` and validates
that `ID`, `CHR` and `MAPINFO` are present, so a malformed manifest produces a
readable error instead of a silent `NULL` propagating into the coordinates.
Both files were confirmed to parse under `fread()` defaults:

```
GPL13534 : 485,577 x 37   ID[1] = cg00035864   CHR[1] = Y   MAPINFO[1] = 8553009
GPL21145 : 868,564 x 14   ID[1] = cg07881041   CHR[1] = 19  MAPINFO[1] = 5236016
```

**Loading `config/paths_local.R`.** Two Windows-specific pitfalls are handled
in `config/config.R`:

* The environment the file is loaded into must have **`baseenv()` as its
  parent, not `emptyenv()`**. `x <- "..."` is an ordinary call to `` `<-` ``,
  which R resolves like any other function; with an `emptyenv()` parent every
  assignment fails with `could not find function "<-"` and the local paths are
  silently never loaded.
* `source(..., encoding = "UTF-8")` is **not** reliable on builds whose native
  locale is not UTF-8 (`LC_CTYPE=C` or `CP936`). In that configuration R
  attempts an iconv translation into the native charset, fails on the first
  non-ASCII byte and raises `unexpected INCOMPLETE_STRING`. Conversely, plain
  `source()` reads the bytes but produces a string that is *valid UTF-8 and
  yet wrong* — a path such as `K:/4DMP_MR_BC/泛癌17NC2020` stops matching its
  own contents. `.load_local_paths()` therefore tries the UTF-8 `source()`
  first and, on failure, falls back to
  `eval(parse(text = readLines(path, encoding = "UTF-8"), encoding = "UTF-8"))`,
  which preserves the bytes. The fallback path was verified to return a
  string that matches `泛癌` literally, whereas plain `source()` did not.

  The practical recommendation for users is unchanged and simpler: keep
  `paths_local.R` **ASCII-only**, and if that is impossible save it as UTF-8
  without a BOM.

**Working directory.** Every script is written to be run from the repository
root (or from the `4DMP-MR-BC.Rproj` project, whose default working directory
is the root): the helper files are pulled in with
`source("helpers/mr_functions.R")` and the entry point with
`source("config/config.R")`. `config.R` itself locates the repository by
walking up from `getwd()` looking for `config/config.R`, so the output
directories are created in the right place regardless.

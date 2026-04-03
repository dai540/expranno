# expranno <img src="man/figures/expranno-logo.svg" align="right" height="88" alt="expranno logo" />

[![pkgdown](https://img.shields.io/badge/docs-pkgdown-315c86)](https://dai540.github.io/expranno/)
[![R-CMD-check](https://github.com/dai540/expranno/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dai540/expranno/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

`expranno` is an R package for RNA-seq workflows that start from an
expression matrix and sample metadata, then move through gene annotation,
metadata integration, Deconvolution, and Signature analysis.

It is intentionally a downstream orchestration layer rather than a full
RNA-seq pipeline. Alignment, counting, QC, and differential expression
stay outside the package. The goal here is to take an existing Ensembl-ID
expression matrix and make it annotation-aware, metadata-aligned, and
ready for deconvolution or signature scoring.

The package is designed for human and mouse expression matrices where:

- the first column of `expr` is `gene_id`
- the remaining columns of `expr` are samples
- the first column of `meta` is `sample`

`expranno` is built around four jobs:

- gene annotation from Ensembl IDs
- expression and metadata integration
- Deconvolution with `immunedeconv`
- Signature analysis with GSVA or ssGSEA

Rather than forcing a single annotation backend, `expranno` uses a
coverage-first strategy that can combine:

- `biomaRt` with a fixed Ensembl release (`v102` by default)
- `org.Hs.eg.db` or `org.Mm.eg.db`
- optional `EnsDb` packages

This makes it possible to recover more symbols, names, and identifiers
than a single-source annotation pass, while still recording which backend
filled each field.

![expranno workflow](man/figures/workflow-overview.svg)

`expranno` does this in a simple analysis sequence:

1. validate `expr` and `meta`
2. annotate genes and write `expr_anno.csv`
3. save `annotation_report.csv`, `annotation_ambiguity.csv`, and `annotation_provenance.csv`
4. merge expression and metadata into `expr_meta_merged.csv`
5. run Deconvolution and save one table per method
6. run Signature analysis and save one score table per method

## Installation

Install from GitHub:

```r
install.packages("pak")
pak::pak("dai540/expranno")
```

Or:

```r
install.packages("remotes")
remotes::install_github("dai540/expranno")
```

Or install from a source tarball:

```r
install.packages("path/to/expranno_2.1.0.tar.gz", repos = NULL, type = "source")
```

Optional backends:

```r
BiocManager::install(c(
  "biomaRt",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "ensembldb",
  "EnsDb.Hsapiens.v86",
  "EnsDb.Mmusculus.v79",
  "immunedeconv",
  "GSVA"
))
```

Then load the package:

```r
library(expranno)
```

## Citation

If you use `expranno`, cite the package as:

> Dai (2026). *expranno: Expression Annotation, Metadata Integration,
> Deconvolution, and Signature Analysis*. R package.
> <https://dai540.github.io/expranno/>

You can also retrieve the citation from R:

```r
citation("expranno")
```

## Author

`expranno` is developed and maintained by Dai Dai.

## Workflow Overview

The main wrapper is `run_expranno()`. It orchestrates the whole pipeline
from raw inputs to saved outputs.

```r
result <- expranno::run_expranno(
  expr = expr,
  meta = meta,
  species = "human",
  annotation_engine = "hybrid",
  biomart_version = 102,
  expr_scale = "abundance",
  duplicate_strategy = "mean",
  output_dir = "results",
  run_deconvolution = TRUE,
  deconv_args = list(indications = rep("SKCM", ncol(expr) - 1L)),
  run_signature = TRUE,
  geneset_file = "hallmark.gmt",
  signature_method = "both",
  signature_kcdf = "Gaussian"
)
```

`deconv_args` is available for method-specific `immunedeconv` arguments.
This is especially important for `timer` and `consensus_tme`, which
require an `indications` vector.

`annotate_expr()` and `run_expranno()` also accept
`SummarizedExperiment`-like inputs. If your assay and metadata already
live in a Bioconductor container, use `assay_name` and `gene_id_col` to
coerce them into the package contract without manual CSV reshaping.

`expr_scale` and `duplicate_strategy` control how duplicated symbols are
collapsed before downstream analyses. The default `"auto"` strategy uses
`"sum"` for count-like data and `"mean"` for abundance or log-scale data.

## Input Contract

```r
demo <- expranno::example_expranno_data()
```

`expr` must be a gene-by-sample table.

```r
head(demo$expr)
```

| gene_id | sample_a | sample_b |
|---|---:|---:|
| ENSG00000141510.17 | 120 | 140 |
| ENSG00000146648.18 | 80 | 77 |

`meta` must be a sample metadata table.

```r
head(demo$meta)
```

| sample | group | batch |
|---|---|---|
| sample_a | case | b1 |
| sample_b | control | b1 |

## Core Functions

`expranno` is intentionally split into one end-to-end wrapper and four
stepwise building blocks.

- `run_expranno()`
- `annotate_expr()`
- `benchmark_annotation_engines()`
- `merge_expr_meta()`
- `run_cell_deconvolution()`
- `run_signature_analysis()`
- `as_expranno_input()`

## Minimal Example

```r
demo <- expranno::example_expranno_data()

res <- expranno::run_expranno(
  expr = demo$expr,
  meta = demo$meta,
  species = "human",
  annotation_engine = "none",
  output_dir = tempdir(),
  run_deconvolution = FALSE,
  run_signature = FALSE
)
```

## Main Outputs

`expranno` is designed to write analysis-ready CSV files with stable names.

- `expr_anno.csv`
- `annotation_report.csv`
- `annotation_ambiguity.csv`
- `annotation_provenance.csv`
- `expr_meta_merged.csv`
- `cell_deconv_<method>.csv`
- `signature_gsva.csv`
- `signature_ssgsea.csv`
- `session_info.txt`

The high-level wrapper is `run_expranno()`, but the package also exposes
smaller steps:

- `annotate_expr()`
- `benchmark_annotation_engines()`
- `merge_expr_meta()`
- `run_cell_deconvolution()`
- `run_signature_analysis()`

## What To Inspect First

After a real run with `annotation_engine = "hybrid"`, the first files and
tables to inspect are:

- `expr_anno.csv`
- `annotation_report.csv`
- `annotation_ambiguity.csv`
- `annotation_provenance.csv`
- `expr_meta_merged.csv`

The practical checklist is:

- confirm `symbol` coverage is high enough for downstream symbol-based tools
- confirm `gene_name` coverage tracks with `symbol`
- inspect any multi-hit cases in `annotation_ambiguity.csv`
- record the Ensembl release, OrgDb package, and run date from `annotation_provenance.csv`
- confirm merged sample labels line up with the expected `meta` rows
- use the merged table as the common input for Deconvolution and Signature analyses

For signature scoring, also confirm that `signature_kcdf` matches the
expression scale you are using. Count-like inputs usually call for
`"Poisson"`, while continuous log-scale inputs are usually better matched
by `"Gaussian"`.

## Documentation

The package website includes:

- a getting started article
- a benchmarking and reproducibility article
- human and mouse step-by-step case studies
- a theory and design article
- a function reference

Website: <https://dai540.github.io/expranno/>

## Backend Smoke Test

For a local optional-backend smoke test, run:

```r
source(system.file("scripts", "run_smoke_optional_backends.R", package = "expranno"))
```

Or run the installed script directly with an output directory:

```r
script <- system.file("scripts", "run_smoke_optional_backends.R", package = "expranno")
system2(file.path(R.home("bin"), "Rscript"), c(script, "smoke-outputs"))
```

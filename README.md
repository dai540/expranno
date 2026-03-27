# expranno <img src="man/figures/expranno-logo.svg" align="right" height="88" alt="expranno logo" />

[![pkgdown](https://img.shields.io/badge/docs-pkgdown-315c86)](https://dai540.github.io/expranno/)
[![R-CMD-check](https://github.com/dai540/expranno/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dai540/expranno/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

`expranno` is an R package for RNA-seq workflows that start from an
expression matrix and sample metadata, then move through gene annotation,
metadata integration, Deconvolution, and Signature analysis.

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

- `biomaRt`
- `org.Hs.eg.db` or `org.Mm.eg.db`
- optional `EnsDb` packages

This makes it possible to recover more symbols, names, and identifiers
than a single-source annotation pass.

![expranno workflow](man/figures/workflow-overview.svg)

`expranno` does this in a simple analysis sequence:

1. validate `expr` and `meta`
2. annotate genes and write `expr_anno.csv`
3. merge expression and metadata into `expr_meta_merged.csv`
4. run Deconvolution and save one table per method
5. run Signature analysis and save one score table per method

## Installation

Install from GitHub:

```r
# install.packages("pak")
pak::pak("dai540/expranno")
```

Local install:

```r
# install.packages("pak")
# pak::pak("path/to/expranno")
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

## Workflow Overview

The main wrapper is `run_expranno()`. It orchestrates the whole pipeline
from raw inputs to saved outputs.

```r
result <- expranno::run_expranno(
  expr = expr,
  meta = meta,
  species = "human",
  annotation_engine = "hybrid",
  output_dir = "results",
  run_deconvolution = TRUE,
  run_signature = TRUE,
  geneset_file = "hallmark.gmt",
  signature_method = "both"
)
```

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
- `merge_expr_meta()`
- `run_cell_deconvolution()`
- `run_signature_analysis()`

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
- `expr_meta_merged.csv`
- `cell_deconv_<method>.csv`
- `signature_gsva.csv`
- `signature_ssgsea.csv`

The high-level wrapper is `run_expranno()`, but the package also exposes
smaller steps:

- `annotate_expr()`
- `merge_expr_meta()`
- `run_cell_deconvolution()`
- `run_signature_analysis()`

## What To Inspect First

After a real run with `annotation_engine = "hybrid"`, the first files and
tables to inspect are:

- `expr_anno.csv`
- the annotation coverage report returned by `annotate_expr()`
- `expr_meta_merged.csv`

The practical checklist is:

- confirm `symbol` coverage is high enough for downstream symbol-based tools
- confirm `gene_name` coverage tracks with `symbol`
- confirm merged sample labels line up with the expected `meta` rows
- use the merged table as the common input for Deconvolution and Signature analyses

## Documentation

The package website includes:

- a getting started article
- human and mouse step-by-step case studies
- a theory and design article
- a function reference

Website: <https://dai540.github.io/expranno/>

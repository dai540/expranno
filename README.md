# expranno <img src="man/figures/logo.svg" align="right" height="88" alt="expranno logo" />

[![pkgdown](https://img.shields.io/badge/docs-pkgdown-315c86)](https://dai540.github.io/expranno/)
[![R-CMD-check](https://github.com/dai540/expranno/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dai540/expranno/actions/workflows/R-CMD-check.yaml)
[![GitHub release](https://img.shields.io/github/v/release/dai540/expranno)](https://github.com/dai540/expranno/releases)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE.txt)

`expranno` is a downstream RNA-seq orchestration package for human and
mouse expression matrices that already exist as an Ensembl-ID table plus
sample metadata. It focuses on four jobs:

- gene annotation
- metadata integration
- Deconvolution with `immunedeconv`
- Signature analysis with GSVA or ssGSEA

The input contract is intentionally strict:

- `expr`: first column is `gene_id`, remaining columns are samples
- `meta`: first column is `sample`

For repeatable annotation, `expranno` ships fixed presets such as
`human_tpm_v102` and `mouse_tpm_v102`, plus bundled truth tables through
`example_annotation_truth()`.

![expranno workflow](man/figures/workflow-overview.svg)

The package validates inputs, annotates genes, records coverage and
provenance, merges expression with metadata, and can optionally run
Deconvolution and Signature scoring.

## Installation

Install from GitHub with `pak`:

```r
install.packages("pak")
pak::pak("dai540/expranno")
```

or `remotes`:

```r
install.packages("remotes")
remotes::install_github("dai540/expranno")
```

Or install from a source tarball:

```r
install.packages("path/to/expranno_<version>.tar.gz", repos = NULL, type = "source")
```

Optional backends can be installed with Bioconductor:

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

Or use the bundled installer script after the package is installed:

```r
source(system.file("scripts", "install_optional_backends.R", package = "expranno"))
expranno_install_optional_backends()
```

For containerized and reproducible installs, `expranno` also ships a
starter Dockerfile at `system.file("docker", "Dockerfile", package =
"expranno")`.

Then load the package:

```r
library(expranno)
```

## Minimal Example

`run_expranno()` is the main wrapper:

```r
demo <- expranno::example_expranno_data()

result <- expranno::run_expranno(
  expr = demo$expr,
  meta = demo$meta,
  annotation_preset = "human_tpm_v102",
  expr_scale = "abundance",
  duplicate_strategy = "mean",
  output_dir = tempdir(),
  run_deconvolution = FALSE,
  run_signature = FALSE
)
```

## Core Functions

- `run_expranno()`
- `annotate_expr()`
- `benchmark_annotation_engines()`
- `validate_annotation_engines()`
- `merge_expr_meta()`
- `run_cell_deconvolution()`
- `run_signature_analysis()`
- `as_expranno_input()`
- `as_expranno_se()`
- `list_annotation_presets()`
- `example_annotation_truth()`

`annotate_expr()` and `run_expranno()` also accept
`SummarizedExperiment`-like inputs through `as_expranno_input()`. To
move annotated results back into Bioconductor containers, use
`as_expranno_se()`.

## Main Outputs

Stable outputs include:

- `expr_anno.csv`
- `annotation_report.csv`
- `annotation_ambiguity.csv`
- `annotation_provenance.csv`
- `annotation_validation_summary.csv`
- `annotation_validation_detail.csv`
- `expr_meta_merged.csv`
- `cell_deconv_<method>.csv`
- `signature_gsva.csv`
- `signature_ssgsea.csv`
- `session_info.txt`

## Presets and Truth Tables

```r
expranno::list_annotation_presets()
expranno::example_annotation_truth("human")
expranno::example_annotation_truth("mouse")
```

## Documentation

Website: <https://dai540.github.io/expranno/>

## Citation

```r
citation("expranno")
```

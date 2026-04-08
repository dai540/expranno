# expranno

[![pkgdown](https://img.shields.io/badge/docs-pkgdown-315c86)](https://dai540.github.io/expranno/)
[![R-CMD-check](https://github.com/dai540/expranno/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dai540/expranno/actions/workflows/R-CMD-check.yaml)
[![GitHub
release](https://img.shields.io/github/v/release/dai540/expranno)](https://github.com/dai540/expranno/releases)
[![License:
MIT](https://img.shields.io/badge/license-MIT-green.svg)](https://dai540.github.io/expranno/LICENSE.txt)

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
[`example_annotation_truth()`](https://dai540.github.io/expranno/reference/example_annotation_truth.md).

![expranno workflow](reference/figures/workflow-overview.svg)

expranno workflow

The package validates inputs, annotates genes, records coverage and
provenance, merges expression with metadata, and can optionally run
Deconvolution and Signature scoring.

## Installation

Install from GitHub with `pak`:

``` r
install.packages("pak")
pak::pak("dai540/expranno")
```

or `remotes`:

``` r
install.packages("remotes")
remotes::install_github("dai540/expranno")
```

Or install from a source tarball:

``` r
install.packages("path/to/expranno_<version>.tar.gz", repos = NULL, type = "source")
```

Optional backends can be installed with Bioconductor:

``` r
BiocManager::install(c(
  "biomaRt",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "ensembldb",
  "EnsDb.Hsapiens.v86",
  "EnsDb.Mmusculus.v79",
  "GSVA"
))
```

Install `immunedeconv` separately from GitHub if you need deconvolution:

``` r
remotes::install_github("omnideconv/immunedeconv")
```

Then load the package:

``` r
library(expranno)
```

## Minimal Example

[`run_expranno()`](https://dai540.github.io/expranno/reference/run_expranno.md)
is the main wrapper:

``` r
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

- [`run_expranno()`](https://dai540.github.io/expranno/reference/run_expranno.md)
- [`annotate_expr()`](https://dai540.github.io/expranno/reference/annotate_expr.md)
- [`benchmark_annotation_engines()`](https://dai540.github.io/expranno/reference/benchmark_annotation_engines.md)
- [`validate_annotation_engines()`](https://dai540.github.io/expranno/reference/validate_annotation_engines.md)
- [`merge_expr_meta()`](https://dai540.github.io/expranno/reference/merge_expr_meta.md)
- [`run_cell_deconvolution()`](https://dai540.github.io/expranno/reference/run_cell_deconvolution.md)
- [`run_signature_analysis()`](https://dai540.github.io/expranno/reference/run_signature_analysis.md)
- [`as_expranno_input()`](https://dai540.github.io/expranno/reference/as_expranno_input.md)
- [`as_expranno_se()`](https://dai540.github.io/expranno/reference/as_expranno_se.md)
- [`list_annotation_presets()`](https://dai540.github.io/expranno/reference/list_annotation_presets.md)
- [`example_annotation_truth()`](https://dai540.github.io/expranno/reference/example_annotation_truth.md)

[`annotate_expr()`](https://dai540.github.io/expranno/reference/annotate_expr.md)
and
[`run_expranno()`](https://dai540.github.io/expranno/reference/run_expranno.md)
also accept `SummarizedExperiment`-like inputs through
[`as_expranno_input()`](https://dai540.github.io/expranno/reference/as_expranno_input.md).
To move annotated results back into Bioconductor containers, use
[`as_expranno_se()`](https://dai540.github.io/expranno/reference/as_expranno_se.md).

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

``` r
expranno::list_annotation_presets()
expranno::example_annotation_truth("human")
expranno::example_annotation_truth("mouse")
```

## Documentation

Website: <https://dai540.github.io/expranno/>

## Citation

``` r
citation("expranno")
```

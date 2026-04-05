# Convert expranno output to SummarizedExperiment

Creates a `SummarizedExperiment` with annotated expression values in the
main assay, annotation fields in `rowData`, and sample metadata in
`colData`. When available, signature and deconvolution outputs are also
appended to `colData` as sample-level columns.

## Usage

``` r
as_expranno_se(
  x,
  assay_name = "expression",
  include_signatures = TRUE,
  include_deconvolution = TRUE
)
```

## Arguments

- x:

  An `expranno_result` or `expranno_annotation`.

- assay_name:

  Assay name used for the expression matrix.

- include_signatures:

  Whether to append signature scores to `colData` when `x` is an
  `expranno_result`.

- include_deconvolution:

  Whether to append deconvolution scores to `colData` when `x` is an
  `expranno_result`.

## Value

A `SummarizedExperiment`.

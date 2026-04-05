# Coerce Bioconductor containers to expranno inputs

Converts a `SummarizedExperiment`-like object into the strict `expr` and
`meta` data-frame contract used by `expranno`.

## Usage

``` r
as_expranno_input(
  x,
  assay_name = NULL,
  gene_id_col = NULL,
  sample_col = "sample"
)
```

## Arguments

- x:

  A `SummarizedExperiment` or `SingleCellExperiment`.

- assay_name:

  Optional assay name or index. Defaults to the first assay.

- gene_id_col:

  Optional row-data column containing Ensembl gene IDs.

- sample_col:

  Column name to use for sample identifiers in the output metadata
  table.

## Value

A named list with `expr` and `meta`.

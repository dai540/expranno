# Run immune deconvolution

Uses `immunedeconv` on an annotated expression matrix. By default, all
available methods except CIBERSORT-family methods are executed.

## Usage

``` r
run_cell_deconvolution(
  expr_anno,
  species = c("auto", "human", "mouse"),
  methods = "all_except_cibersort",
  gene_column = "symbol",
  output_dir = NULL,
  prefix = "cell_deconv_",
  verbose = TRUE
)
```

## Arguments

- expr_anno:

  Annotated expression table.

- species:

  Either `"auto"`, `"human"`, or `"mouse"`.

- methods:

  Either `"all_except_cibersort"` or a character vector of explicit
  method names.

- gene_column:

  Gene symbol column to use for deconvolution.

- output_dir:

  Optional directory to write one CSV per method.

- prefix:

  File name prefix for output CSV files.

- verbose:

  Whether to emit progress messages.

## Value

A named list of deconvolution result tables.

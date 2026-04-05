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
  expr_scale = c("auto", "count", "abundance", "log"),
  duplicate_strategy = c("auto", "sum", "mean", "max", "first"),
  output_dir = NULL,
  prefix = "cell_deconv_",
  verbose = TRUE,
  ...
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

- expr_scale:

  Expression scale. This affects duplicate symbol handling and triggers
  a warning for `expr_scale = "log"`.

- duplicate_strategy:

  Strategy used when multiple rows map to the same symbol. `"auto"` uses
  `"sum"` for counts and `"mean"` otherwise.

- output_dir:

  Optional directory to write one CSV per method.

- prefix:

  File name prefix for output CSV files.

- verbose:

  Whether to emit progress messages.

- ...:

  Additional arguments forwarded to `immunedeconv`, including
  method-specific values such as `indications`.

## Value

A named list of deconvolution result tables.

# Run GSVA or ssGSEA signature scoring

Uses a symbol-based expression matrix derived from `expr_anno` and a
user-supplied gene set database.

## Usage

``` r
run_signature_analysis(
  expr_anno,
  geneset_file = NULL,
  gene_sets = NULL,
  method = c("gsva", "ssgsea", "both"),
  gene_column = "symbol",
  output_dir = NULL,
  prefix = "signature_"
)
```

## Arguments

- expr_anno:

  Annotated expression table.

- geneset_file:

  Optional GMT path.

- gene_sets:

  Optional named list of gene sets.

- method:

  One of `"gsva"`, `"ssgsea"`, or `"both"`.

- gene_column:

  Symbol column to use.

- output_dir:

  Optional output directory.

- prefix:

  File prefix.

## Value

A named list with score matrices as data frames.

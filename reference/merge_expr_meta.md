# Merge annotated expression and metadata

Converts the annotated expression matrix to long format, joins it to
`meta` by `sample`, and optionally writes the merged table as CSV.

## Usage

``` r
merge_expr_meta(expr_anno, meta, output_file = NULL, value_name = "expression")
```

## Arguments

- expr_anno:

  Annotated expression table returned by
  [`annotate_expr()`](https://example.com/expranno/reference/annotate_expr.md).

- meta:

  Checked or user-supplied metadata table with `sample` first.

- output_file:

  Optional CSV path.

- value_name:

  Name of the expression value column in the long table.

## Value

A merged `data.frame`.

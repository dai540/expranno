# Getting Started with expranno

``` r
library(expranno)
demo <- example_expranno_data()
```

`expranno` expects:

- an expression matrix-like data frame with `gene_id` in the first
  column
- sample columns after `gene_id`
- a metadata table with `sample` in the first column

``` r
annotated <- annotate_expr(
  expr = demo$expr,
  meta = demo$meta,
  species = "human",
  annotation_engine = "none"
)

merged <- merge_expr_meta(
  expr_anno = annotated$expr_anno,
  meta = annotated$meta_checked
)
```

The all-in-one wrapper is
[`run_expranno()`](https://example.com/expranno/reference/run_expranno.md).

``` r
result <- run_expranno(
  expr = demo$expr,
  meta = demo$meta,
  species = "human",
  annotation_engine = "none",
  output_dir = tempdir(),
  run_deconvolution = FALSE,
  run_signature = FALSE
)

result
#> <expranno_result>
#>   annotated genes: 3
#>   merged rows: 9
#>   deconvolution runs: 0
#>   signature runs: 0
```

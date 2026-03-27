# Validate an expression table

`expr` must be a `data.frame` with `gene_id` as the first column and
sample columns after it.

## Usage

``` r
validate_expr(expr)
```

## Arguments

- expr:

  A gene-by-sample expression table.

## Value

`TRUE` invisibly on success.

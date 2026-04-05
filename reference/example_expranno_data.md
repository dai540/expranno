# Example expression and metadata inputs

Returns a small human or mouse RNA-seq example that follows the package
input contract. The toy matrix is suitable for examples, tests, and
documentation.

## Usage

``` r
example_expranno_data(species = c("human", "mouse"))
```

## Arguments

- species:

  Either `"human"` or `"mouse"`.

## Value

A named list with `expr` and `meta`.

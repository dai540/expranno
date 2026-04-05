# Example annotation truth tables

Loads small bundled truth tables for human or mouse annotation
validation. These tables are intended for examples, tests, and
reproducible benchmark demonstrations with
[`validate_annotation_engines()`](https://dai540.github.io/expranno/reference/validate_annotation_engines.md).

## Usage

``` r
example_annotation_truth(species = c("human", "mouse"))
```

## Arguments

- species:

  Either `"human"` or `"mouse"`.

## Value

A `data.frame` with `gene_id`, `symbol`, `gene_name`, and `biotype`.

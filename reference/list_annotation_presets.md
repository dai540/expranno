# List built-in annotation presets

Returns the fixed annotation presets that ship with `expranno`. These
presets are intended to make repeated human and mouse annotation runs
easier to standardize across projects.

## Usage

``` r
list_annotation_presets()
```

## Value

A `data.frame` describing the available presets, including
expression-scale intent, backend order, and bundled truth-table
resource.

# Annotate an expression matrix

Validates `expr` and `meta`, normalizes Ensembl IDs, infers or accepts
the species, runs a coverage-first annotation workflow, binds annotation
columns to the expression table, and optionally writes `expr_anno.csv`.

## Usage

``` r
annotate_expr(
  expr,
  meta,
  species = c("auto", "human", "mouse"),
  annotation_engine = c("hybrid", "biomart", "orgdb", "ensdb", "none"),
  fields = NULL,
  strip_version = TRUE,
  output_file = NULL,
  verbose = TRUE
)
```

## Arguments

- expr:

  Expression table with `gene_id` in the first column.

- meta:

  Metadata table with `sample` in the first column.

- species:

  Either `"auto"`, `"human"`, or `"mouse"`.

- annotation_engine:

  Annotation backend strategy. Use `"hybrid"` to try `biomaRt`, `orgdb`,
  and `EnsDb` in sequence.

- fields:

  Optional annotation fields to keep.

- strip_version:

  Whether to remove Ensembl version suffixes.

- output_file:

  Optional CSV path for the annotated matrix.

- verbose:

  Whether to emit progress messages.

## Value

An `expranno_annotation` object.

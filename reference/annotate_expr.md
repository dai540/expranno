# Annotate an expression matrix

Validates input, normalizes Ensembl IDs, infers or accepts the species,
runs a coverage-first annotation workflow, binds annotation columns to
the expression table, and optionally writes `expr_anno.csv` plus
provenance and ambiguity reports.

## Usage

``` r
annotate_expr(
  expr,
  meta = NULL,
  species = c("auto", "human", "mouse"),
  annotation_preset = NULL,
  annotation_engine = c("hybrid", "biomart", "orgdb", "ensdb", "none"),
  fields = NULL,
  strip_version = TRUE,
  biomart_version = 102,
  biomart_host = NULL,
  biomart_mirror = NULL,
  assay_name = NULL,
  gene_id_col = NULL,
  sample_col = "sample",
  output_file = NULL,
  report_file = NULL,
  ambiguity_file = NULL,
  provenance_file = NULL,
  verbose = TRUE
)
```

## Arguments

- expr:

  Expression table with `gene_id` in the first column, or a
  `SummarizedExperiment`-like object.

- meta:

  Metadata table with `sample` in the first column. Leave `NULL` when
  `expr` is a `SummarizedExperiment` and metadata should be taken from
  `colData`.

- species:

  Either `"auto"`, `"human"`, or `"mouse"`.

- annotation_preset:

  Optional preset that fixes a reproducible annotation configuration.
  Supported values are `"human_v102"`, `"mouse_v102"`,
  `"human_tpm_v102"`, `"mouse_tpm_v102"`, `"human_count_v102"`, and
  `"mouse_count_v102"`.

- annotation_engine:

  Annotation backend strategy. Use `"hybrid"` to try `biomaRt`, `orgdb`,
  and `EnsDb` in sequence.

- fields:

  Optional annotation fields to keep.

- strip_version:

  Whether to remove Ensembl version suffixes.

- biomart_version:

  Fixed Ensembl release used by the `biomaRt` backend.

- biomart_host:

  Optional explicit Ensembl host.

- biomart_mirror:

  Optional Ensembl mirror name.

- assay_name:

  Optional assay name when `expr` is a `SummarizedExperiment`.

- gene_id_col:

  Optional row-data column containing Ensembl IDs when `expr` is a
  `SummarizedExperiment`.

- sample_col:

  Metadata column to use as `sample` when `expr` is a
  `SummarizedExperiment`.

- output_file:

  Optional CSV path for the annotated matrix.

- report_file:

  Optional CSV path for the coverage report.

- ambiguity_file:

  Optional CSV path for the ambiguity report.

- provenance_file:

  Optional CSV path for backend provenance.

- verbose:

  Whether to emit progress messages.

## Value

An `expranno_annotation` object.

# Validate annotation engines against a truth table

Runs one or more annotation engines, then compares the chosen annotation
fields against a user-supplied truth table keyed by Ensembl gene ID.

## Usage

``` r
validate_annotation_engines(
  expr,
  meta = NULL,
  truth,
  truth_gene_col = "gene_id",
  species = c("auto", "human", "mouse"),
  annotation_preset = NULL,
  engines = c("biomart", "orgdb", "ensdb", "hybrid"),
  fields = NULL,
  strip_version = TRUE,
  biomart_version = 102,
  biomart_host = NULL,
  biomart_mirror = NULL,
  assay_name = NULL,
  gene_id_col = NULL,
  sample_col = "sample",
  output_file = NULL,
  detail_file = NULL,
  verbose = TRUE
)
```

## Arguments

- expr:

  Expression table with `gene_id` first, or a
  `SummarizedExperiment`-like object.

- meta:

  Metadata table with `sample` first. Leave `NULL` when `expr` is a
  `SummarizedExperiment`.

- truth:

  A data frame containing a `gene_id` column and one or more truth
  annotation fields such as `symbol` or `gene_name`.

- truth_gene_col:

  Column in `truth` containing the Ensembl gene ID.

- species:

  Either `"auto"`, `"human"`, or `"mouse"`.

- annotation_preset:

  Optional preset that fixes a reproducible annotation configuration.
  Supported values are `"human_v102"`, `"mouse_v102"`,
  `"human_tpm_v102"`, `"mouse_tpm_v102"`, `"human_count_v102"`, and
  `"mouse_count_v102"`.

- engines:

  Character vector of annotation engines to compare.

- fields:

  Optional annotation fields to validate. Defaults to the intersection
  between `truth` and the package's standard annotation fields.

- strip_version:

  Whether to remove Ensembl version suffixes in both prediction and
  truth tables.

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

  Optional CSV path for the validation summary.

- detail_file:

  Optional CSV path for per-gene validation details.

- verbose:

  Whether to emit progress messages.

## Value

An `expranno_validation` object.

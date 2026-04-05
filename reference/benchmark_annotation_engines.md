# Benchmark annotation engines on the same input

Runs one or more annotation engines on the same expression and metadata
inputs, then summarizes coverage, ambiguity, and provenance for
comparison.

## Usage

``` r
benchmark_annotation_engines(
  expr,
  meta = NULL,
  species = c("auto", "human", "mouse"),
  annotation_preset = NULL,
  engines = c("none", "biomart", "orgdb", "ensdb", "hybrid"),
  fields = NULL,
  strip_version = TRUE,
  biomart_version = 102,
  biomart_host = NULL,
  biomart_mirror = NULL,
  assay_name = NULL,
  gene_id_col = NULL,
  sample_col = "sample",
  output_file = NULL,
  coverage_file = NULL,
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

  Optional annotation fields to benchmark.

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

  Optional CSV path for the benchmark summary.

- coverage_file:

  Optional CSV path for long-format coverage results.

- verbose:

  Whether to emit progress messages.

## Value

An `expranno_benchmark` object.

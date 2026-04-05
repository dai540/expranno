# Run the full expranno workflow

This wrapper annotates genes, writes `expr_anno.csv`, writes provenance
and ambiguity reports, merges expression with metadata, optionally runs
deconvolution and signature scoring, optionally benchmarks annotation
engines, and returns all outputs in one structured object.

## Usage

``` r
run_expranno(
  expr,
  meta = NULL,
  species = c("auto", "human", "mouse"),
  annotation_preset = NULL,
  annotation_engine = c("hybrid", "biomart", "orgdb", "ensdb", "none"),
  output_dir = ".",
  biomart_version = 102,
  biomart_host = NULL,
  biomart_mirror = NULL,
  assay_name = NULL,
  gene_id_col = NULL,
  sample_col = "sample",
  run_deconvolution = TRUE,
  deconv_methods = "all_except_cibersort",
  expr_scale = c("auto", "count", "abundance", "log"),
  duplicate_strategy = c("auto", "sum", "mean", "max", "first"),
  deconv_args = list(),
  run_signature = TRUE,
  geneset_file = NULL,
  gene_sets = NULL,
  signature_method = c("gsva", "ssgsea", "both"),
  signature_kcdf = c("auto", "Gaussian", "Poisson", "none"),
  signature_min_size = 1,
  signature_max_size = Inf,
  gsva_args = list(),
  ssgsea_args = list(),
  run_benchmark = FALSE,
  benchmark_engines = c("none", "biomart", "orgdb", "ensdb", "hybrid"),
  run_validation = FALSE,
  validation_truth = NULL,
  validation_truth_gene_col = "gene_id",
  validation_fields = NULL,
  save_session_info = TRUE,
  verbose = TRUE
)
```

## Arguments

- expr:

  Expression table with `gene_id` in the first column, or a
  `SummarizedExperiment`-like object.

- meta:

  Metadata table with `sample` in the first column. Leave `NULL` when
  `expr` is a `SummarizedExperiment`.

- species:

  Either `"auto"`, `"human"`, or `"mouse"`.

- annotation_preset:

  Optional preset that fixes a reproducible annotation configuration.
  Supported values are `"human_v102"`, `"mouse_v102"`,
  `"human_tpm_v102"`, `"mouse_tpm_v102"`, `"human_count_v102"`, and
  `"mouse_count_v102"`.

- annotation_engine:

  Annotation backend strategy.

- output_dir:

  Output directory.

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

- run_deconvolution:

  Whether to run `immunedeconv`.

- deconv_methods:

  Either `"all_except_cibersort"` or an explicit vector.

- expr_scale:

  Expression scale used to choose duplicate-symbol handling.

- duplicate_strategy:

  Strategy used when multiple rows map to the same symbol. `"auto"` uses
  `"sum"` for counts and `"mean"` otherwise.

- deconv_args:

  Optional named list of additional arguments passed to
  [`run_cell_deconvolution()`](https://dai540.github.io/expranno/reference/run_cell_deconvolution.md),
  including method-specific values such as `indications`.

- run_signature:

  Whether to run GSVA or ssGSEA.

- geneset_file:

  Optional GMT path.

- gene_sets:

  Optional named list of gene sets.

- signature_method:

  One of `"gsva"`, `"ssgsea"`, or `"both"`.

- signature_kcdf:

  Kernel choice passed to
  [`run_signature_analysis()`](https://dai540.github.io/expranno/reference/run_signature_analysis.md).

- signature_min_size:

  Minimum gene-set size after ID matching.

- signature_max_size:

  Maximum gene-set size after ID matching.

- gsva_args:

  Optional named list of extra arguments passed to `GSVA::gsvaParam()`.

- ssgsea_args:

  Optional named list of extra arguments passed to
  `GSVA::ssgseaParam()`.

- run_benchmark:

  Whether to run
  [`benchmark_annotation_engines()`](https://dai540.github.io/expranno/reference/benchmark_annotation_engines.md)
  on the same inputs.

- benchmark_engines:

  Annotation engines to benchmark when `run_benchmark = TRUE`.

- run_validation:

  Whether to run
  [`validate_annotation_engines()`](https://dai540.github.io/expranno/reference/validate_annotation_engines.md)
  against a supplied truth table.

- validation_truth:

  Optional truth table keyed by Ensembl gene ID.

- validation_truth_gene_col:

  Column in `validation_truth` containing the Ensembl gene ID.

- validation_fields:

  Optional truth fields to validate, such as `symbol`.

- save_session_info:

  Whether to write `session_info.txt`.

- verbose:

  Whether to emit progress messages.

## Value

An `expranno_result` object.

# Run the full expranno workflow

This wrapper annotates genes, writes `expr_anno.csv`, merges expression
with metadata, optionally runs immune deconvolution and signature
scoring, and returns all outputs in one structured object.

## Usage

``` r
run_expranno(
  expr,
  meta,
  species = c("auto", "human", "mouse"),
  annotation_engine = c("hybrid", "biomart", "orgdb", "ensdb", "none"),
  output_dir = ".",
  run_deconvolution = TRUE,
  deconv_methods = "all_except_cibersort",
  run_signature = TRUE,
  geneset_file = NULL,
  gene_sets = NULL,
  signature_method = c("gsva", "ssgsea", "both"),
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

  Annotation backend strategy.

- output_dir:

  Output directory.

- run_deconvolution:

  Whether to run `immunedeconv`.

- deconv_methods:

  Either `"all_except_cibersort"` or an explicit vector.

- run_signature:

  Whether to run GSVA or ssGSEA.

- geneset_file:

  Optional GMT path.

- gene_sets:

  Optional named list of gene sets.

- signature_method:

  One of `"gsva"`, `"ssgsea"`, or `"both"`.

- verbose:

  Whether to emit progress messages.

## Value

An `expranno_result` object.

# Run GSVA or ssGSEA signature scoring

Uses a symbol-based expression matrix derived from `expr_anno` and a
user-supplied gene set database.

## Usage

``` r
run_signature_analysis(
  expr_anno,
  geneset_file = NULL,
  gene_sets = NULL,
  method = c("gsva", "ssgsea", "both"),
  gene_column = "symbol",
  expr_scale = c("auto", "count", "abundance", "log"),
  duplicate_strategy = c("auto", "sum", "mean", "max", "first"),
  kcdf = c("auto", "Gaussian", "Poisson", "none"),
  min_size = 1,
  max_size = Inf,
  gsva_args = list(),
  ssgsea_args = list(),
  output_dir = NULL,
  prefix = "signature_"
)
```

## Arguments

- expr_anno:

  Annotated expression table.

- geneset_file:

  Optional GMT path.

- gene_sets:

  Optional named list of gene sets.

- method:

  One of `"gsva"`, `"ssgsea"`, or `"both"`.

- gene_column:

  Symbol column to use.

- expr_scale:

  Expression scale. This affects duplicate symbol handling.

- duplicate_strategy:

  Strategy used when multiple rows map to the same symbol. `"auto"` uses
  `"sum"` for counts and `"mean"` otherwise.

- kcdf:

  Kernel choice passed to `GSVA::gsvaParam()`.

- min_size:

  Minimum gene-set size after ID matching.

- max_size:

  Maximum gene-set size after ID matching.

- gsva_args:

  Optional named list of extra arguments passed to `GSVA::gsvaParam()`.

- ssgsea_args:

  Optional named list of extra arguments passed to
  `GSVA::ssgseaParam()`.

- output_dir:

  Optional output directory.

- prefix:

  File prefix.

## Value

A named list with score matrices as data frames.

# Benchmarking and Reproducibility

`expranno` aims to maximize annotation coverage, but public-facing
analysis also needs provenance and benchmarkable behavior.

This article focuses on three things:

1.  comparing annotation engines on the same input
2.  validating exact matches against a truth table
3.  recording backend provenance for reproducibility
4.  keeping a session-level record of the software environment

``` r
library(expranno)
demo <- example_expranno_data("human")
truth <- example_annotation_truth("human")
```

## Benchmark annotation engines

Use
[`benchmark_annotation_engines()`](https://dai540.github.io/expranno/reference/benchmark_annotation_engines.md)
when you want to compare `biomart`, `orgdb`, `ensdb`, and `hybrid` on
exactly the same data.

``` r
benchmark <- benchmark_annotation_engines(
  expr = demo$expr,
  meta = demo$meta,
  species = "human",
  engines = c("none"),
  verbose = FALSE
)

benchmark$summary
#>      engine status species message annotated_genes total_genes
#> none   none     ok   human    <NA>               0           3
#>      ambiguous_gene_fields symbol_coverage gene_name_coverage
#> none                     0               0                  0
```

In a real run with optional backends installed, a more informative call
looks like this:

``` r
benchmark <- benchmark_annotation_engines(
  expr = demo$expr,
  meta = demo$meta,
  species = "human",
  engines = c("biomart", "orgdb", "ensdb", "hybrid"),
  biomart_version = 102
)
```

The benchmark output helps answer practical questions:

- does `hybrid` improve `symbol` coverage over `biomart` alone?
- which engine contributes most to `gene_name` coverage?
- how many ambiguous gene-fields appear under each engine?

## Validate against truth

Coverage is useful, but a public package also needs a way to ask whether
recovered labels are actually correct for a known subset.

``` r
validation <- validate_annotation_engines(
  expr = demo$expr,
  meta = demo$meta,
  truth = truth,
  species = "human",
  engines = "none",
  fields = "symbol",
  verbose = FALSE
)

validation$summary
#>              engine  field truth_rows matched_rows missing_prediction_rows
#> none::symbol   none symbol          3            0                       3
#>              mismatch_rows ambiguous_rows match_rate missing_prediction_rate
#> none::symbol             0              0          0                       1
#>              mismatch_rate
#> none::symbol             0
```

In a production run, the same function can be used to compare `biomart`,
`orgdb`, `ensdb`, and `hybrid` against a curated truth set from the same
species.

The package-level helpers now make that easy to reproduce in examples:

- `example_annotation_truth("human")`
- `example_annotation_truth("mouse")`

## Provenance from annotation

[`annotate_expr()`](https://dai540.github.io/expranno/reference/annotate_expr.md)
now records backend provenance and ambiguity directly.

``` r
annotation <- annotate_expr(
  expr = demo$expr,
  meta = demo$meta,
  species = "human",
  annotation_engine = "none",
  verbose = FALSE
)

annotation$provenance
#>   source species dataset package_name backend_release host mirror genes_queried
#> 1   none   human    none         none            <NA> <NA>   <NA>             3
#>   rows_returned annotation_date  status
#> 1             0      2026-04-05 skipped
#>                                      message
#> 1 Annotation lookup skipped by user request.
annotation$ambiguity_report
#> [1] gene_id         field           chosen_value    candidate_count
#> [5] candidates      source         
#> <0 rows> (or 0-length row.names)
```

For production runs with `annotation_engine = "hybrid"`, the provenance
table becomes much more informative. It records which backend ran, what
species and dataset were used, which Ensembl release was targeted, and
the annotation date.

## Stable files from the full wrapper

[`run_expranno()`](https://dai540.github.io/expranno/reference/run_expranno.md)
writes reproducibility-oriented side files by default.

``` r
outdir <- tempfile("expranno-benchmark-")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

result <- run_expranno(
  expr = demo$expr,
  meta = demo$meta,
  species = "human",
  annotation_engine = "none",
  output_dir = outdir,
  run_deconvolution = FALSE,
  run_signature = FALSE,
  run_benchmark = TRUE,
  benchmark_engines = c("none"),
  verbose = FALSE
)

basename(result$files)
#> [1] "annotation_ambiguity.csv"          "annotation_benchmark_coverage.csv"
#> [3] "annotation_benchmark_summary.csv"  "annotation_provenance.csv"        
#> [5] "annotation_report.csv"             "expr_anno.csv"                    
#> [7] "expr_meta_merged.csv"              "session_info.txt"
```

The wrapper now keeps:

- `annotation_report.csv`
- `annotation_ambiguity.csv`
- `annotation_provenance.csv`
- `annotation_benchmark_summary.csv`
- `annotation_benchmark_coverage.csv`
- `annotation_validation_summary.csv`
- `annotation_validation_detail.csv`
- `session_info.txt`

These files make it easier to re-run the same workflow later and explain
which annotation sources were active in a given result directory.

For a more reproducible machine-level setup, `expranno` also ships:

- `inst/scripts/install_optional_backends.R`
- `inst/docker/Dockerfile`

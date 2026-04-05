# Getting Started with expranno

``` r
library(expranno)
demo <- example_expranno_data()
```

`expranno` expects:

- an expression matrix-like data frame with `gene_id` in the first
  column
- sample columns after `gene_id`
- a metadata table with `sample` in the first column

``` mermaid
flowchart TB
  subgraph inputs["Inputs"]
    expr["expr
gene x sample table"]
    meta["meta
sample x metadata table"]
  end

  subgraph core["Core workflow"]
    validate["Validate
check required columns"]
    annotate["Annotation
normalize Ensembl IDs
hybrid human or mouse mapping"]
    report["Coverage report
annotation_rate by field"]
    merge["Merge
combine expr_anno and meta"]
  end

  subgraph outputs["Saved outputs"]
    expranno["expr_anno.csv"]
    merged["expr_meta_merged.csv"]
  end

  subgraph downstream["Downstream analyses"]
    deconv["Deconvolution
immunedeconv
one CSV per method"]
    signature["Signature
GSVA or ssGSEA
one CSV per method"]
  end

  expr --> validate
  meta --> validate
  validate --> annotate
  annotate --> expranno
  annotate --> report
  expranno --> merge
  meta --> merge
  merge --> merged
  merged --> deconv
  merged --> signature

  classDef input fill:#eef5fb,stroke:#315c86,color:#17324d,stroke-width:1.5px;
  classDef process fill:#ffffff,stroke:#315c86,color:#17324d,stroke-width:1.5px;
  classDef output fill:#eaf6ef,stroke:#3c7a57,color:#17324d,stroke-width:1.5px;

  class expr,meta input;
  class validate,annotate,report,merge,deconv,signature process;
  class expranno,merged output;
```

## The workflow at a glance

The workflow is easiest to read in five checkpoints:

1.  validate the input contract
2.  annotate genes and inspect coverage
3.  inspect provenance and ambiguity
4.  merge expression with metadata
5.  run Deconvolution or Signature analyses from the merged table

If you want fixed human or mouse defaults instead of manually setting
species, version stripping, and Ensembl release each time, use one of
the built-in presets.

``` r
list_annotation_presets()
#>   annotation_preset species  recommended_input expr_scale annotation_engine
#> 1        human_v102   human                any       auto            hybrid
#> 2        mouse_v102   mouse                any       auto            hybrid
#> 3    human_tpm_v102   human TPM-like abundance  abundance            hybrid
#> 4    mouse_tpm_v102   mouse TPM-like abundance  abundance            hybrid
#> 5  human_count_v102   human         raw counts      count            hybrid
#> 6  mouse_count_v102   mouse         raw counts      count            hybrid
#>   strip_version biomart_version                   symbol_priority
#> 1          TRUE             102 hgnc_symbol -> external_gene_name
#> 2          TRUE             102  mgi_symbol -> external_gene_name
#> 3          TRUE             102 hgnc_symbol -> external_gene_name
#> 4          TRUE             102  mgi_symbol -> external_gene_name
#> 5          TRUE             102 hgnc_symbol -> external_gene_name
#> 6          TRUE             102  mgi_symbol -> external_gene_name
#>                                   fallback_order
#> 1  biomaRt -> org.Hs.eg.db -> EnsDb.Hsapiens.v86
#> 2 biomaRt -> org.Mm.eg.db -> EnsDb.Mmusculus.v79
#> 3  biomaRt -> org.Hs.eg.db -> EnsDb.Hsapiens.v86
#> 4 biomaRt -> org.Mm.eg.db -> EnsDb.Mmusculus.v79
#> 5  biomaRt -> org.Hs.eg.db -> EnsDb.Hsapiens.v86
#> 6 biomaRt -> org.Mm.eg.db -> EnsDb.Mmusculus.v79
#>                       bundled_truth
#> 1 example_annotation_truth('human')
#> 2 example_annotation_truth('mouse')
#> 3 example_annotation_truth('human')
#> 4 example_annotation_truth('mouse')
#> 5 example_annotation_truth('human')
#> 6 example_annotation_truth('mouse')
```

If you want the preset choices in a more explicit comparison table,
including recommended input scale, backend cascade, and the matching
bundled truth resource, see the preset reference article.

The bundled truth resources are useful when you want a reproducible
validation example without preparing your own table first.

``` r
example_annotation_truth("human")
#>              gene_id symbol                        gene_name        biotype
#> 1 ENSG00000141510.17   TP53                tumor protein p53 protein_coding
#> 2 ENSG00000146648.18   EGFR epidermal growth factor receptor protein_coding
#> 3 ENSG00000012048.23  BRCA1      BRCA1 DNA repair associated protein_coding
```

## Pick a species

Use `species = "human"` for `ENSG...` input IDs and `species = "mouse"`
for `ENSMUSG...` input IDs. If IDs are clean Ensembl IDs,
`species = "auto"` can infer this for you.

``` r
human_demo <- example_expranno_data("human")
mouse_demo <- example_expranno_data("mouse")

head(human_demo$expr$gene_id)
#> [1] "ENSG00000141510.17" "ENSG00000146648.18" "ENSG00000012048.23"
head(mouse_demo$expr$gene_id)
#> [1] "ENSMUSG00000059552.8"  "ENSMUSG00000020122.15" "ENSMUSG00000017167.16"
```

## Pick an annotation engine

`expranno` exposes five annotation modes:

- `"hybrid"`: recommended production mode
- `"biomart"`: use only `biomaRt`
- `"orgdb"`: use only `org.Hs.eg.db` or `org.Mm.eg.db`
- `"ensdb"`: use only `EnsDb`
- `"none"`: skip annotation lookup and keep the normalized IDs only

`"hybrid"` is the default because it is the most coverage-oriented mode.
It uses a cascade:

1.  `biomaRt` with a fixed Ensembl release (`v102` by default)
2.  species-specific `orgdb`
3.  species-specific `EnsDb`

That mirrors a practical manual workflow: query Ensembl first, fill
unmapped IDs from `org.Hs.eg.db` or `org.Mm.eg.db`, then use `EnsDb` to
enrich structural fields such as biotype and coordinates.

For lightweight docs and tests, the examples below use `"none"`.

For real runs, a reproducible preset is often easier to share than a
long list of arguments. Examples:

- `annotation_preset = "human_tpm_v102"`
- `annotation_preset = "mouse_tpm_v102"`
- `annotation_preset = "human_count_v102"`
- `annotation_preset = "mouse_count_v102"`

The preset table above is the quickest way to compare recommended input
scale, backend order, and bundled truth resource.

``` r
annotated <- annotate_expr(
  expr = demo$expr,
  meta = demo$meta,
  species = "human",
  annotation_engine = "none"
)

merged <- merge_expr_meta(
  expr_anno = annotated$expr_anno,
  meta = annotated$meta_checked
)
```

The all-in-one wrapper is
[`run_expranno()`](https://dai540.github.io/expranno/reference/run_expranno.md).

``` r
result <- run_expranno(
  expr = demo$expr,
  meta = demo$meta,
  species = "human",
  annotation_engine = "none",
  output_dir = tempdir(),
  run_deconvolution = FALSE,
  run_signature = FALSE
)

result
#> <expranno_result>
#>   annotated genes: 3
#>   merged rows: 9
#>   deconvolution runs: 0
#>   signature runs: 0
#>   benchmark runs: 0
#>   validation runs: 0
```

For a reproducible production-style run, switch to a preset explicitly:

``` r
run_expranno(
  expr = demo$expr,
  meta = demo$meta,
  annotation_preset = "human_count_v102",
  output_dir = "results"
)
```

After a run, you can also return to a Bioconductor-native container.

``` r
annotated <- annotate_expr(
  expr = demo$expr,
  meta = demo$meta,
  species = "human",
  annotation_engine = "none",
  verbose = FALSE
)

if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  se <- as_expranno_se(annotated)
  se
} else {
  "Install SummarizedExperiment to use as_expranno_se()."
}
#> class: SummarizedExperiment 
#> dim: 3 3 
#> metadata(1): expranno
#> assays(1): expression
#> rownames(3): ENSG00000141510 ENSG00000146648 ENSG00000012048
#> rowData names(48): gene_id_raw gene_id ... annotation_backend_mirror
#>   annotation_date
#> colnames(3): sample_a sample_b sample_c
#> colData names(4): sample group batch species
```

If your data already live in a `SummarizedExperiment`, you can start
there directly.

``` r
result <- run_expranno(
  expr = se,
  species = "human",
  assay_name = "counts",
  gene_id_col = "gene_id",
  annotation_engine = "hybrid"
)
```

## Control duplicate symbols and score kernels

Downstream methods usually operate on gene symbols, so duplicated
symbols have to be collapsed. `expranno` now makes that rule explicit.

- `expr_scale = "count"` usually pairs with `duplicate_strategy = "sum"`
- `expr_scale = "abundance"` or `"log"` usually pairs with
  `duplicate_strategy = "mean"`

For signature scoring, `signature_kcdf` should also match the expression
scale:

- use `"Poisson"` for count-like input
- use `"Gaussian"` for continuous log-scale input
- use `"auto"` if you want GSVA to choose

## Method-specific deconvolution arguments

Some `immunedeconv` methods require extra arguments. In particular,
`timer` and `consensus_tme` need an `indications` vector with one tumor
type per sample.

``` r
run_expranno(
  expr = demo$expr,
  meta = demo$meta,
  species = "human",
  annotation_engine = "hybrid",
  output_dir = "results",
  run_deconvolution = TRUE,
  deconv_args = list(indications = c("SKCM", "SKCM")),
  run_signature = TRUE,
  geneset_file = "hallmark.gmt",
  signature_method = "both",
  signature_kcdf = "Gaussian"
)
```

## What to inspect after a run

For a production analysis with `"hybrid"`, the first things to check
are:

- `expr_anno.csv`
- `annotation_report.csv`
- `annotation_ambiguity.csv`
- `annotation_provenance.csv`
- `expr_meta_merged.csv`

If `symbol` coverage is weak, downstream Deconvolution and Signature
methods will usually be the first places where the problem appears.

If you want to compare engines on the same input before standardizing a
lab workflow, benchmark them explicitly.

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

If you also have a truth table, you can validate exact matches instead
of only coverage.

``` r
truth <- data.frame(
  gene_id = demo$expr$gene_id,
  symbol = c("TP53", "EGFR", "BRCA1"),
  stringsAsFactors = FALSE
)

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

For reproducible examples or CI runs, you do not need to write that
truth table by hand. `expranno` ships with small human and mouse truth
resources that line up with the demo Ensembl IDs.

``` r
example_annotation_truth("human")
#>              gene_id symbol                        gene_name        biotype
#> 1 ENSG00000141510.17   TP53                tumor protein p53 protein_coding
#> 2 ENSG00000146648.18   EGFR epidermal growth factor receptor protein_coding
#> 3 ENSG00000012048.23  BRCA1      BRCA1 DNA repair associated protein_coding
```

If your workflow starts and ends in Bioconductor containers, convert the
result back into a `SummarizedExperiment`.

``` r
se <- as_expranno_se(result)
```

That keeps the annotated expression matrix in the assay, writes
annotation fields to `rowData`, and appends sample-level Deconvolution
or Signature outputs to `colData`.

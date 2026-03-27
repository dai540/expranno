# expranno <img src="man/figures/expranno-logo.svg" align="right" height="88" alt="expranno logo" />

`expranno` is an R package for RNA-seq workflows that start from an
expression matrix and sample metadata, then move through gene annotation,
metadata integration, Deconvolution, and Signature analysis.

The package is designed for human and mouse expression matrices where:

- the first column of `expr` is `gene_id`
- the remaining columns of `expr` are samples
- the first column of `meta` is `sample`

`expranno` is built around four jobs:

- gene annotation from Ensembl IDs
- expression and metadata integration
- Deconvolution with `immunedeconv`
- Signature analysis with GSVA or ssGSEA

Rather than forcing a single annotation backend, `expranno` uses a
coverage-first strategy that can combine:

- `biomaRt`
- `org.Hs.eg.db` or `org.Mm.eg.db`
- optional `EnsDb` packages

This makes it possible to recover more symbols, names, and identifiers
than a single-source annotation pass.

<style>
.mermaid {
  margin: 1.5rem auto 2rem auto;
  max-width: 980px;
  padding: 0.75rem;
  background: #f7f9fc;
  border: 1px solid #d9e4ef;
  border-radius: 16px;
}
</style>
<pre class="mermaid">
flowchart TB
  subgraph inputs["Inputs"]
    expr["expr<br/>gene x sample table"]
    meta["meta<br/>sample x metadata table"]
  end

  subgraph core["Core workflow"]
    validate["Validate<br/>check required columns"]
    annotate["Annotation<br/>normalize Ensembl IDs<br/>hybrid human or mouse mapping"]
    report["Coverage report<br/>annotation_rate by field"]
    merge["Merge<br/>combine expr_anno and meta"]
  end

  subgraph outputs["Saved outputs"]
    expranno["expr_anno.csv"]
    merged["expr_meta_merged.csv"]
  end

  subgraph downstream["Downstream analyses"]
    deconv["Deconvolution<br/>immunedeconv<br/>one CSV per method"]
    signature["Signature<br/>GSVA or ssGSEA<br/>one CSV per method"]
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
</pre>
<script src="docs/deps/mermaid/mermaid.min.js"></script>
<script>
if (window.mermaid) {
  window.mermaid.initialize({
    startOnLoad: true,
    securityLevel: "loose",
    theme: "base",
    themeVariables: {
      primaryColor: "#ffffff",
      primaryTextColor: "#17324d",
      primaryBorderColor: "#315c86",
      lineColor: "#315c86",
      secondaryColor: "#eef5fb",
      tertiaryColor: "#eaf6ef",
      clusterBkg: "#f7f9fc",
      clusterBorder: "#cfddeb",
      fontFamily: "Lato, Arial, sans-serif"
    },
    flowchart: {
      curve: "basis",
      htmlLabels: true
    }
  });
  window.mermaid.run({ querySelector: ".mermaid" });
}
</script>

`expranno` does this in a simple analysis sequence:

1. validate `expr` and `meta`
2. annotate genes and write `expr_anno.csv`
3. merge expression and metadata into `expr_meta_merged.csv`
4. run Deconvolution and save one table per method
5. run Signature analysis and save one score table per method

## Installation

Local install:

```r
# install.packages("pak")
# pak::pak("path/to/expranno")
```

Optional backends:

```r
BiocManager::install(c(
  "biomaRt",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "ensembldb",
  "EnsDb.Hsapiens.v86",
  "EnsDb.Mmusculus.v79",
  "immunedeconv",
  "GSVA"
))
```

## Workflow Overview

The main wrapper is `run_expranno()`. It orchestrates the whole pipeline
from raw inputs to saved outputs.

```r
result <- expranno::run_expranno(
  expr = expr,
  meta = meta,
  species = "human",
  annotation_engine = "hybrid",
  output_dir = "results",
  run_deconvolution = TRUE,
  run_signature = TRUE,
  geneset_file = "hallmark.gmt",
  signature_method = "both"
)
```

## Input Contract

```r
demo <- expranno::example_expranno_data()
```

`expr` must be a gene-by-sample table.

```r
head(demo$expr)
```

| gene_id | sample_a | sample_b |
|---|---:|---:|
| ENSG00000141510.17 | 120 | 140 |
| ENSG00000146648.18 | 80 | 77 |

`meta` must be a sample metadata table.

```r
head(demo$meta)
```

| sample | group | batch |
|---|---|---|
| sample_a | case | b1 |
| sample_b | control | b1 |

## Core Functions

`expranno` is intentionally split into one end-to-end wrapper and four
stepwise building blocks.

- `run_expranno()`
- `annotate_expr()`
- `merge_expr_meta()`
- `run_cell_deconvolution()`
- `run_signature_analysis()`

## Minimal Example

```r
demo <- expranno::example_expranno_data()

res <- expranno::run_expranno(
  expr = demo$expr,
  meta = demo$meta,
  species = "human",
  annotation_engine = "none",
  output_dir = tempdir(),
  run_deconvolution = FALSE,
  run_signature = FALSE
)
```

## Main Outputs

`expranno` is designed to write analysis-ready CSV files with stable names.

- `expr_anno.csv`
- `expr_meta_merged.csv`
- `cell_deconv_<method>.csv`
- `signature_gsva.csv`
- `signature_ssgsea.csv`

The high-level wrapper is `run_expranno()`, but the package also exposes
smaller steps:

- `annotate_expr()`
- `merge_expr_meta()`
- `run_cell_deconvolution()`
- `run_signature_analysis()`

## What To Inspect First

After a real run with `annotation_engine = "hybrid"`, the first files and
tables to inspect are:

- `expr_anno.csv`
- the annotation coverage report returned by `annotate_expr()`
- `expr_meta_merged.csv`

The practical checklist is:

- confirm `symbol` coverage is high enough for downstream symbol-based tools
- confirm `gene_name` coverage tracks with `symbol`
- confirm merged sample labels line up with the expected `meta` rows
- use the merged table as the common input for Deconvolution and Signature analyses

## Documentation

The package website includes:

- a getting started article
- a step-by-step case study tutorial
- a theory and design article
- a function reference

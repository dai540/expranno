# expranno

`expranno` is a minimal R package for downstream RNA-seq workflows built around
two strict input tables:

- `expr`: the first column is `gene_id`, the remaining columns are samples
- `meta`: the first column is `sample`, the remaining columns are metadata

The package is intentionally narrow. It does not perform alignment, read-level
quality control, count generation, normalization theory, or differential
expression analysis. It starts after an expression matrix already exists and
focuses on making that matrix easier to annotate, merge, and pass into common
downstream tools.

The current source tree is deliberately small. Large example data, generated
site output, temporary build products, smoke logs, tarballs, and `.Rcheck`
directories are not kept in the repository.

Documentation site: <https://dai540.github.io/expranno/>

## Installation

Install the package from GitHub:

```r
install.packages("pak")
pak::pak("dai540/expranno")
```

or:

```r
install.packages("remotes")
remotes::install_github("dai540/expranno")
```

Optional backends are kept out of `Imports` so the package stays light.
Install only the backends you need.

Annotation backends:

```r
BiocManager::install(c(
  "biomaRt",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "org.Mm.eg.db"
))
```

Signature scoring:

```r
BiocManager::install("GSVA")
```

Deconvolution:

```r
remotes::install_github("omnideconv/immunedeconv")
```

## Package scope

`expranno` standardizes four downstream steps:

1. gene annotation from Ensembl identifiers
2. expression and metadata integration
3. deconvolution through `immunedeconv`
4. GSVA or ssGSEA signature scoring through `GSVA`

The package targets human and mouse expression matrices. Species can be set
explicitly or inferred from the Ensembl ID prefix:

- `ENSG...` -> human
- `ENSMUSG...` -> mouse

To keep the repository and installed package small, the package ships only tiny
code-generated example tables. It does not bundle large demo data or download
large datasets as part of examples or documentation.

## Core functions

- `annotate_expr()` annotates an expression table and optionally writes
  `expr_anno.csv` and `annotation_report.csv`
- `merge_expr_meta()` converts the annotated matrix into a long table and joins
  it to `meta`
- `run_cell_deconvolution()` runs selected `immunedeconv` methods on the
  annotated matrix
- `run_signature_analysis()` runs GSVA or ssGSEA from a gene-set list or GMT
  file
- `run_expranno()` orchestrates the full workflow from `expr` and `meta`

Support functions:

- `validate_expr()`
- `validate_meta()`
- `read_genesets()`
- `example_expranno_data()`

## Input contract

The package is strict about input shape.

For `expr`:

- first column must be named `gene_id`
- all remaining columns must be numeric
- sample column names must be unique

For `meta`:

- first column must be named `sample`
- each sample ID must be unique

This is intentional. The package uses a narrow contract so downstream behavior
stays predictable and CSV outputs remain stable.

## Annotation design

`annotate_expr()` supports four engines:

- `"hybrid"`: `biomaRt` first, then `OrgDb` fallback
- `"biomart"`: `biomaRt` only
- `"orgdb"`: `OrgDb` only
- `"none"`: no external annotation, useful for tests and lightweight examples

By default, version suffixes such as `.17` are stripped from Ensembl gene IDs.
The output always keeps both the raw ID and the cleaned ID:

- `gene_id_raw`
- `gene_id`
- `symbol`
- `gene_name`
- `annotation_source`
- `annotation_status`

The annotation report summarizes:

- inferred or requested species
- total rows
- unique Ensembl IDs
- annotated rows
- unannotated rows
- annotation rate

## Output files

When `output_dir` is supplied, the package writes stable file names:

- `expr_anno.csv`
- `annotation_report.csv`
- `expr_meta_merged.csv`
- `cell_deconv_<method>.csv`
- `signature_gsva.csv`
- `signature_ssgsea.csv`

## Minimal example

```r
library(expranno)

demo <- example_expranno_data("human")

result <- run_expranno(
  expr = demo$expr,
  meta = demo$meta,
  annotation_engine = "none",
  output_dir = tempdir(),
  run_deconvolution = FALSE,
  run_signature = FALSE
)

result$annotation$report
head(result$expr_meta_merged)
```

This example uses `annotation_engine = "none"` so it runs without optional
annotation packages. For real annotation, use `"hybrid"`, `"biomart"`, or
`"orgdb"` after installing the corresponding backends.

## Documentation map

The pkgdown site is organized into four top-level sections:

- `Getting Started`: package scope, installation, and a first run
- `Guides`: focused guidance for annotation and analysis choices
- `Tutorials`: step-by-step human and mouse walkthroughs using small built-in
  examples
- `Reference`: function-level documentation generated from the source

## Design decisions

The package is intentionally conservative in a few places.

1. Optional heavy backends are not installed automatically.
   This keeps the package small and avoids shipping unnecessary dependencies.

2. The repository does not keep generated site files.
   `pkgdown` is built in GitHub Actions and published from `gh-pages`.

3. The package does not bundle large datasets.
   Examples use tiny code-generated tables so the source tree stays compact.

4. The codebase is split into only a few source files.
   This keeps navigation simple and avoids excessive file fragmentation.

## Repository layout

- `R/package.R`: package documentation and imports
- `R/utils.R`: validation and internal helpers
- `R/annotate.R`: example data and annotation workflow
- `R/analysis.R`: merge, deconvolution, signature scoring, and main wrapper
- `tests/`: minimal unit tests
- `vignettes/`: pkgdown articles for Getting Started, Guides, and Tutorials


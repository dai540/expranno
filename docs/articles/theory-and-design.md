# Theory and Design

`expranno` is designed around a CSV-first RNA-seq analysis workflow.

## Design Goals

The package aims to do four things well:

1.  enforce a simple and strict input contract
2.  maximize annotation coverage for human and mouse
3.  keep expression and metadata aligned
4.  make downstream immune and signature analyses reproducible

## Why the package starts with annotation

Many downstream tools expect gene symbols rather than Ensembl IDs. A
reliable annotation layer makes it possible to:

- keep the original `gene_id`
- recover a display-ready `symbol`
- attach gene-level metadata such as name and biotype
- create symbol-based matrices for deconvolution and scoring

## Annotation strategy

The default design is a coverage-first cascade:

1.  `biomaRt`
2.  `org.Hs.eg.db` or `org.Mm.eg.db`
3.  optional `EnsDb` packages

Each stage only fills missing fields left by the previous stage. This is
why `expranno` can often recover more usable annotation than a
single-source lookup.

## Data shape choices

`expranno` keeps two complementary representations:

- a wide annotated matrix for file export
- a long merged table for downstream analysis

The wide representation is best for exchange and reproducibility:

- `expr_anno.csv`

The long representation is best for plotting, grouped summaries, and
sample-level joins:

- `expr_meta_merged.csv`

## Downstream analyses

After annotation and merging, the package supports:

- immune deconvolution with `immunedeconv`
- pathway or signature scoring with GSVA or ssGSEA

Both are designed to write explicit CSV outputs so that results are easy
to inspect outside R.

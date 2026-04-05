# Preset Reference

`expranno` ships with fixed annotation presets so repeated human and
mouse workflows can share the same release, species, and version
stripping rules.

``` r
library(expranno)

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

## How to read the preset table

- `recommended_input` tells you which expression scale the preset is
  designed around.
- `biomart_version` fixes the targeted Ensembl release.
- `fallback_order` shows the human or mouse backend cascade used by the
  hybrid engine.
- `bundled_truth` tells you which built-in truth table pairs naturally
  with the preset for validation examples.

## Typical choices

For human TPM workflows:

``` r
result <- run_expranno(
  expr = expr,
  meta = meta,
  annotation_preset = "human_tpm_v102",
  expr_scale = "abundance",
  duplicate_strategy = "mean",
  run_deconvolution = TRUE,
  run_signature = TRUE,
  geneset_file = "hallmark.gmt",
  signature_kcdf = "Gaussian"
)
```

For mouse count workflows:

``` r
result <- run_expranno(
  expr = expr,
  meta = meta,
  annotation_preset = "mouse_count_v102",
  expr_scale = "count",
  duplicate_strategy = "sum",
  run_deconvolution = TRUE,
  run_signature = TRUE,
  geneset_file = "hallmark.gmt",
  signature_kcdf = "Poisson"
)
```

## Built-in truth resources

The package now bundles small human and mouse truth tables that match
the example Ensembl IDs and are useful for demos, tests, and
reproducible validation tutorials.

``` r
example_annotation_truth("human")
#>              gene_id symbol                        gene_name        biotype
#> 1 ENSG00000141510.17   TP53                tumor protein p53 protein_coding
#> 2 ENSG00000146648.18   EGFR epidermal growth factor receptor protein_coding
#> 3 ENSG00000012048.23  BRCA1      BRCA1 DNA repair associated protein_coding
example_annotation_truth("mouse")
#>                 gene_id symbol                         gene_name        biotype
#> 1  ENSMUSG00000059552.8  Trp53 transformation related protein 53 protein_coding
#> 2 ENSMUSG00000020122.15   Egfr  epidermal growth factor receptor protein_coding
#> 3 ENSMUSG00000017167.16  Brca1                   breast cancer 1 protein_coding
```

These are intentionally small. They are not meant to replace a curated
project-specific truth set, but they make the validation workflow easy
to reproduce across machines and in documentation.

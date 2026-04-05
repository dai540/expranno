# Case Study: Human RNA-seq Workflow

## Goal

This case study shows a simple human RNA-seq workflow in `expranno`,
starting from `expr` and `meta`, then moving through annotation,
integration, and output generation.

## Step 1: Input tables

``` r
demo$expr
#>              gene_id sample_a sample_b sample_c
#> 1 ENSG00000141510.17      120      140      118
#> 2 ENSG00000146648.18       80       77       91
#> 3 ENSG00000012048.23       25       30       21
demo$meta
#>     sample   group batch species
#> 1 sample_a    case    b1   human
#> 2 sample_b control    b1   human
#> 3 sample_c    case    b2   human
```

## Step 2: Validate the input contract

``` r
validate_expr(demo$expr)
validate_meta(demo$meta)
```

## Step 3: Choose an annotation engine

For human data, the recommended production choice is:

- `annotation_preset = "human_tpm_v102"` for TPM-like input
- `annotation_preset = "human_count_v102"` for count-like input

That setting tries:

1.  `biomaRt`
2.  `org.Hs.eg.db`
3.  optional `EnsDb.Hsapiens.v86`

For the website build, we use `"none"` so the article stays lightweight.
To still show what a successful annotation pass looks like, the article
also includes a small mock output table and coverage report.

``` r
anno <- annotate_expr(
  expr = demo$expr,
  meta = demo$meta,
  species = "human",
  annotation_engine = "none"
)

anno$expr_anno
#>          gene_id_raw         gene_id symbol symbol_candidates
#> 1 ENSG00000141510.17 ENSG00000141510   <NA>              <NA>
#> 2 ENSG00000146648.18 ENSG00000146648   <NA>              <NA>
#> 3 ENSG00000012048.23 ENSG00000012048   <NA>              <NA>
#>   symbol_candidate_count symbol_is_ambiguous symbol_source gene_name
#> 1                      0               FALSE          <NA>      <NA>
#> 2                      0               FALSE          <NA>      <NA>
#> 3                      0               FALSE          <NA>      <NA>
#>   gene_name_candidates gene_name_candidate_count gene_name_is_ambiguous
#> 1                 <NA>                         0                  FALSE
#> 2                 <NA>                         0                  FALSE
#> 3                 <NA>                         0                  FALSE
#>   gene_name_source entrez_id entrez_id_candidates entrez_id_candidate_count
#> 1             <NA>      <NA>                 <NA>                         0
#> 2             <NA>      <NA>                 <NA>                         0
#> 3             <NA>      <NA>                 <NA>                         0
#>   entrez_id_is_ambiguous entrez_id_source biotype biotype_candidates
#> 1                  FALSE             <NA>    <NA>               <NA>
#> 2                  FALSE             <NA>    <NA>               <NA>
#> 3                  FALSE             <NA>    <NA>               <NA>
#>   biotype_candidate_count biotype_is_ambiguous biotype_source chromosome
#> 1                       0                FALSE           <NA>       <NA>
#> 2                       0                FALSE           <NA>       <NA>
#> 3                       0                FALSE           <NA>       <NA>
#>   chromosome_candidates chromosome_candidate_count chromosome_is_ambiguous
#> 1                  <NA>                          0                   FALSE
#> 2                  <NA>                          0                   FALSE
#> 3                  <NA>                          0                   FALSE
#>   chromosome_source start start_candidates start_candidate_count
#> 1              <NA>  <NA>             <NA>                     0
#> 2              <NA>  <NA>             <NA>                     0
#> 3              <NA>  <NA>             <NA>                     0
#>   start_is_ambiguous start_source  end end_candidates end_candidate_count
#> 1              FALSE         <NA> <NA>           <NA>                   0
#> 2              FALSE         <NA> <NA>           <NA>                   0
#> 3              FALSE         <NA> <NA>           <NA>                   0
#>   end_is_ambiguous end_source strand strand_candidates strand_candidate_count
#> 1            FALSE       <NA>   <NA>              <NA>                      0
#> 2            FALSE       <NA>   <NA>              <NA>                      0
#> 3            FALSE       <NA>   <NA>              <NA>                      0
#>   strand_is_ambiguous strand_source annotation_source annotation_status
#> 1               FALSE          <NA>              <NA>       unannotated
#> 2               FALSE          <NA>              <NA>       unannotated
#> 3               FALSE          <NA>              <NA>       unannotated
#>   annotation_backend_release annotation_backend_host annotation_backend_mirror
#> 1                       <NA>                    <NA>                      <NA>
#> 2                       <NA>                    <NA>                      <NA>
#> 3                       <NA>                    <NA>                      <NA>
#>   annotation_date sample_a sample_b sample_c
#> 1      2026-04-05      120      140      118
#> 2      2026-04-05       80       77       91
#> 3      2026-04-05       25       30       21
```

## Step 4: Inspect `expr_anno.csv`

In a real run, this table is what gets written to `expr_anno.csv`. Below
is a compact mock excerpt that shows the kind of annotated output you
should expect from a successful human run with `"hybrid"`.

``` r
mock_expr_anno
#>          gene_id_raw         gene_id symbol
#> 1 ENSG00000141510.17 ENSG00000141510   TP53
#> 2 ENSG00000146648.18 ENSG00000146648   EGFR
#> 3 ENSG00000157764.15 ENSG00000157764   BRAF
#>                                       gene_name annotation_status
#> 1                             tumor protein p53         annotated
#> 2              epidermal growth factor receptor         annotated
#> 3 B-Raf proto-oncogene, serine/threonine kinase         annotated
```

The main columns to check are:

- `gene_id_raw`
- `gene_id`
- `symbol`
- `gene_name`
- `annotation_status`

## Step 5: Inspect the annotation report

`expranno` returns a compact report table that summarizes how many genes
received each annotation field.

``` r
mock_report
#>       field annotation_rate
#> 1    symbol            1.00
#> 2 gene_name            1.00
#> 3 entrez_id            0.67
#> 4   biotype            1.00
```

## Step 6: Inspect provenance

In a real hybrid run, `annotation_provenance.csv` records which backend
ran and which release was used.

``` r
mock_provenance
#>    source           backend_release annotation_date status
#> 1 biomart              Ensembl v102      2026-04-03     ok
#> 2   orgdb       org.Hs.eg.db 3.22.0      2026-04-03     ok
#> 3   ensdb EnsDb.Hsapiens.v86 2.99.0      2026-04-03     ok
```

## Step 7: How to read annotation coverage

The `annotation_rate` column is the proportion of genes with a
non-missing value for each field in the annotated gene table.

The main practical checks are:

1.  `symbol` coverage should be high enough for downstream symbol-based
    tools
2.  `gene_name` coverage should usually track with symbol coverage
3.  `entrez_id` can be a little lower, but large gaps still deserve
    attention
4.  low rates suggest rerunning with `"hybrid"` and checking installed
    backends

In this mock example, `symbol` and `gene_name` are fully covered, while
`entrez_id` coverage is partial. That pattern is often acceptable if the
next steps mainly depend on symbols, but it would still be worth noting
in the analysis record.

## Step 8: Merge annotated expression with metadata

``` r
merged <- merge_expr_meta(
  expr_anno = anno$expr_anno,
  meta = anno$meta_checked
)

merged
#>     sample        gene_id_raw         gene_id symbol symbol_candidates
#> 1 sample_a ENSG00000141510.17 ENSG00000141510   <NA>              <NA>
#> 2 sample_a ENSG00000146648.18 ENSG00000146648   <NA>              <NA>
#> 3 sample_a ENSG00000012048.23 ENSG00000012048   <NA>              <NA>
#> 4 sample_b ENSG00000141510.17 ENSG00000141510   <NA>              <NA>
#> 5 sample_b ENSG00000146648.18 ENSG00000146648   <NA>              <NA>
#> 6 sample_b ENSG00000012048.23 ENSG00000012048   <NA>              <NA>
#> 7 sample_c ENSG00000141510.17 ENSG00000141510   <NA>              <NA>
#> 8 sample_c ENSG00000146648.18 ENSG00000146648   <NA>              <NA>
#> 9 sample_c ENSG00000012048.23 ENSG00000012048   <NA>              <NA>
#>   symbol_candidate_count symbol_is_ambiguous symbol_source gene_name
#> 1                      0               FALSE          <NA>      <NA>
#> 2                      0               FALSE          <NA>      <NA>
#> 3                      0               FALSE          <NA>      <NA>
#> 4                      0               FALSE          <NA>      <NA>
#> 5                      0               FALSE          <NA>      <NA>
#> 6                      0               FALSE          <NA>      <NA>
#> 7                      0               FALSE          <NA>      <NA>
#> 8                      0               FALSE          <NA>      <NA>
#> 9                      0               FALSE          <NA>      <NA>
#>   gene_name_candidates gene_name_candidate_count gene_name_is_ambiguous
#> 1                 <NA>                         0                  FALSE
#> 2                 <NA>                         0                  FALSE
#> 3                 <NA>                         0                  FALSE
#> 4                 <NA>                         0                  FALSE
#> 5                 <NA>                         0                  FALSE
#> 6                 <NA>                         0                  FALSE
#> 7                 <NA>                         0                  FALSE
#> 8                 <NA>                         0                  FALSE
#> 9                 <NA>                         0                  FALSE
#>   gene_name_source entrez_id entrez_id_candidates entrez_id_candidate_count
#> 1             <NA>      <NA>                 <NA>                         0
#> 2             <NA>      <NA>                 <NA>                         0
#> 3             <NA>      <NA>                 <NA>                         0
#> 4             <NA>      <NA>                 <NA>                         0
#> 5             <NA>      <NA>                 <NA>                         0
#> 6             <NA>      <NA>                 <NA>                         0
#> 7             <NA>      <NA>                 <NA>                         0
#> 8             <NA>      <NA>                 <NA>                         0
#> 9             <NA>      <NA>                 <NA>                         0
#>   entrez_id_is_ambiguous entrez_id_source biotype biotype_candidates
#> 1                  FALSE             <NA>    <NA>               <NA>
#> 2                  FALSE             <NA>    <NA>               <NA>
#> 3                  FALSE             <NA>    <NA>               <NA>
#> 4                  FALSE             <NA>    <NA>               <NA>
#> 5                  FALSE             <NA>    <NA>               <NA>
#> 6                  FALSE             <NA>    <NA>               <NA>
#> 7                  FALSE             <NA>    <NA>               <NA>
#> 8                  FALSE             <NA>    <NA>               <NA>
#> 9                  FALSE             <NA>    <NA>               <NA>
#>   biotype_candidate_count biotype_is_ambiguous biotype_source chromosome
#> 1                       0                FALSE           <NA>       <NA>
#> 2                       0                FALSE           <NA>       <NA>
#> 3                       0                FALSE           <NA>       <NA>
#> 4                       0                FALSE           <NA>       <NA>
#> 5                       0                FALSE           <NA>       <NA>
#> 6                       0                FALSE           <NA>       <NA>
#> 7                       0                FALSE           <NA>       <NA>
#> 8                       0                FALSE           <NA>       <NA>
#> 9                       0                FALSE           <NA>       <NA>
#>   chromosome_candidates chromosome_candidate_count chromosome_is_ambiguous
#> 1                  <NA>                          0                   FALSE
#> 2                  <NA>                          0                   FALSE
#> 3                  <NA>                          0                   FALSE
#> 4                  <NA>                          0                   FALSE
#> 5                  <NA>                          0                   FALSE
#> 6                  <NA>                          0                   FALSE
#> 7                  <NA>                          0                   FALSE
#> 8                  <NA>                          0                   FALSE
#> 9                  <NA>                          0                   FALSE
#>   chromosome_source start start_candidates start_candidate_count
#> 1              <NA>  <NA>             <NA>                     0
#> 2              <NA>  <NA>             <NA>                     0
#> 3              <NA>  <NA>             <NA>                     0
#> 4              <NA>  <NA>             <NA>                     0
#> 5              <NA>  <NA>             <NA>                     0
#> 6              <NA>  <NA>             <NA>                     0
#> 7              <NA>  <NA>             <NA>                     0
#> 8              <NA>  <NA>             <NA>                     0
#> 9              <NA>  <NA>             <NA>                     0
#>   start_is_ambiguous start_source  end end_candidates end_candidate_count
#> 1              FALSE         <NA> <NA>           <NA>                   0
#> 2              FALSE         <NA> <NA>           <NA>                   0
#> 3              FALSE         <NA> <NA>           <NA>                   0
#> 4              FALSE         <NA> <NA>           <NA>                   0
#> 5              FALSE         <NA> <NA>           <NA>                   0
#> 6              FALSE         <NA> <NA>           <NA>                   0
#> 7              FALSE         <NA> <NA>           <NA>                   0
#> 8              FALSE         <NA> <NA>           <NA>                   0
#> 9              FALSE         <NA> <NA>           <NA>                   0
#>   end_is_ambiguous end_source strand strand_candidates strand_candidate_count
#> 1            FALSE       <NA>   <NA>              <NA>                      0
#> 2            FALSE       <NA>   <NA>              <NA>                      0
#> 3            FALSE       <NA>   <NA>              <NA>                      0
#> 4            FALSE       <NA>   <NA>              <NA>                      0
#> 5            FALSE       <NA>   <NA>              <NA>                      0
#> 6            FALSE       <NA>   <NA>              <NA>                      0
#> 7            FALSE       <NA>   <NA>              <NA>                      0
#> 8            FALSE       <NA>   <NA>              <NA>                      0
#> 9            FALSE       <NA>   <NA>              <NA>                      0
#>   strand_is_ambiguous strand_source annotation_source annotation_status
#> 1               FALSE          <NA>              <NA>       unannotated
#> 2               FALSE          <NA>              <NA>       unannotated
#> 3               FALSE          <NA>              <NA>       unannotated
#> 4               FALSE          <NA>              <NA>       unannotated
#> 5               FALSE          <NA>              <NA>       unannotated
#> 6               FALSE          <NA>              <NA>       unannotated
#> 7               FALSE          <NA>              <NA>       unannotated
#> 8               FALSE          <NA>              <NA>       unannotated
#> 9               FALSE          <NA>              <NA>       unannotated
#>   annotation_backend_release annotation_backend_host annotation_backend_mirror
#> 1                       <NA>                    <NA>                      <NA>
#> 2                       <NA>                    <NA>                      <NA>
#> 3                       <NA>                    <NA>                      <NA>
#> 4                       <NA>                    <NA>                      <NA>
#> 5                       <NA>                    <NA>                      <NA>
#> 6                       <NA>                    <NA>                      <NA>
#> 7                       <NA>                    <NA>                      <NA>
#> 8                       <NA>                    <NA>                      <NA>
#> 9                       <NA>                    <NA>                      <NA>
#>   annotation_date expression   group batch species
#> 1      2026-04-05        120    case    b1   human
#> 2      2026-04-05         80    case    b1   human
#> 3      2026-04-05         25    case    b1   human
#> 4      2026-04-05        140 control    b1   human
#> 5      2026-04-05         77 control    b1   human
#> 6      2026-04-05         30 control    b1   human
#> 7      2026-04-05        118    case    b2   human
#> 8      2026-04-05         91    case    b2   human
#> 9      2026-04-05         21    case    b2   human
```

## Step 9: Inspect a simple sample summary

``` r
aggregate(expression ~ sample + group, data = merged, FUN = mean)
#>     sample   group expression
#> 1 sample_a    case   75.00000
#> 2 sample_c    case   76.66667
#> 3 sample_b control   82.33333
```

## Step 10: Run the standard wrapper

``` r
out <- run_expranno(
  expr = demo$expr,
  meta = demo$meta,
  species = "human",
  annotation_engine = "none",
  output_dir = tempdir(),
  run_deconvolution = FALSE,
  run_signature = FALSE
)

out
#> <expranno_result>
#>   annotated genes: 3
#>   merged rows: 9
#>   deconvolution runs: 0
#>   signature runs: 0
#>   benchmark runs: 0
#>   validation runs: 0
```

## Human workflow takeaway

For a real human analysis, the usual sequence is:

1.  annotate with a fixed preset such as `"human_tpm_v102"`
2.  inspect `expr_anno.csv`
3.  inspect `annotation_report.csv`
4.  inspect `annotation_provenance.csv`
5.  merge with `meta`
6.  run deconvolution
7.  run GSVA or ssGSEA

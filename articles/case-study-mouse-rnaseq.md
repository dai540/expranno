# Case Study: Mouse RNA-seq Workflow

## Goal

This case study shows the same workflow for mouse RNA-seq data. The key
difference is that the species and mouse annotation backends need to be
set correctly.

## Step 1: Input tables

``` r
demo$expr
#>                 gene_id sample_a sample_b sample_c
#> 1  ENSMUSG00000059552.8      220      205      214
#> 2 ENSMUSG00000020122.15      150      163      158
#> 3 ENSMUSG00000017167.16       44       47       40
demo$meta
#>     sample   group batch species
#> 1 sample_a    case    b1   mouse
#> 2 sample_b control    b1   mouse
#> 3 sample_c    case    b2   mouse
```

## Step 2: Validate the input contract

``` r
validate_expr(demo$expr)
validate_meta(demo$meta)
```

## Step 3: Choose an annotation engine

For mouse data, the recommended production choice is:

- `annotation_preset = "mouse_tpm_v102"` for TPM-like input
- `annotation_preset = "mouse_count_v102"` for count-like input

That setting tries:

1.  `biomaRt`
2.  `org.Mm.eg.db`
3.  optional `EnsDb.Mmusculus.v79`

This matters because mouse Ensembl IDs begin with `ENSMUSG`, and the
symbol mapping must come from mouse-specific resources. As in the human
article, a small mock output is shown below so the docs can display a
nonzero coverage example without heavy annotation backends.

``` r
anno <- annotate_expr(
  expr = demo$expr,
  meta = demo$meta,
  species = "mouse",
  annotation_engine = "none"
)

anno$expr_anno
#>             gene_id_raw            gene_id symbol symbol_candidates
#> 1  ENSMUSG00000059552.8 ENSMUSG00000059552   <NA>              <NA>
#> 2 ENSMUSG00000020122.15 ENSMUSG00000020122   <NA>              <NA>
#> 3 ENSMUSG00000017167.16 ENSMUSG00000017167   <NA>              <NA>
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
#> 1      2026-04-05      220      205      214
#> 2      2026-04-05      150      163      158
#> 3      2026-04-05       44       47       40
```

## Step 4: Inspect `expr_anno.csv`

``` r
mock_expr_anno
#>             gene_id_raw            gene_id symbol
#> 1  ENSMUSG00000059552.8 ENSMUSG00000059552  Trp53
#> 2  ENSMUSG00000064341.1 ENSMUSG00000064341 mt-Nd1
#> 3 ENSMUSG00000017167.15 ENSMUSG00000017167   Egfr
#>                                      gene_name annotation_status
#> 1            transformation related protein 53         annotated
#> 2 mitochondrially encoded NADH dehydrogenase 1         annotated
#> 3             epidermal growth factor receptor         annotated
```

For mouse data, this is where you check that the normalized IDs still
follow the expected `ENSMUSG...` pattern and that downstream annotation
fields were filled when a real backend is used.

## Step 5: Inspect the annotation report

``` r
mock_report
#>       field annotation_rate
#> 1    symbol            1.00
#> 2 gene_name            1.00
#> 3 entrez_id            0.67
#> 4   biotype            1.00
```

## Step 6: Inspect provenance

``` r
mock_provenance
#>    source            backend_release annotation_date status
#> 1 biomart               Ensembl v102      2026-04-03     ok
#> 2   orgdb        org.Mm.eg.db 3.22.0      2026-04-03     ok
#> 3   ensdb EnsDb.Mmusculus.v79 2.99.0      2026-04-03     ok
```

## Step 7: How to read annotation coverage

For mouse annotation, the report helps answer three questions quickly:

1.  was `symbol` coverage high enough for downstream tools?
2.  were `gene_name` and `entrez_id` filled as expected?
3.  do the fields look consistent with mouse-specific resources?

In this mock example, the key symbol and gene name fields are fully
covered, while `entrez_id` is only partially covered. That is a useful
signal: the annotation is good enough for symbol-based Deconvolution and
Signature workflows, but not every identifier family is equally
complete.

If the report is unexpectedly weak, check:

- `species = "mouse"`
- mouse annotation backends are installed
- the input IDs are real mouse Ensembl gene IDs

## Step 8: Merge annotated expression with metadata

``` r
merged <- merge_expr_meta(
  expr_anno = anno$expr_anno,
  meta = anno$meta_checked
)

merged
#>     sample           gene_id_raw            gene_id symbol symbol_candidates
#> 1 sample_a  ENSMUSG00000059552.8 ENSMUSG00000059552   <NA>              <NA>
#> 2 sample_a ENSMUSG00000020122.15 ENSMUSG00000020122   <NA>              <NA>
#> 3 sample_a ENSMUSG00000017167.16 ENSMUSG00000017167   <NA>              <NA>
#> 4 sample_b  ENSMUSG00000059552.8 ENSMUSG00000059552   <NA>              <NA>
#> 5 sample_b ENSMUSG00000020122.15 ENSMUSG00000020122   <NA>              <NA>
#> 6 sample_b ENSMUSG00000017167.16 ENSMUSG00000017167   <NA>              <NA>
#> 7 sample_c  ENSMUSG00000059552.8 ENSMUSG00000059552   <NA>              <NA>
#> 8 sample_c ENSMUSG00000020122.15 ENSMUSG00000020122   <NA>              <NA>
#> 9 sample_c ENSMUSG00000017167.16 ENSMUSG00000017167   <NA>              <NA>
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
#> 1      2026-04-05        220    case    b1   mouse
#> 2      2026-04-05        150    case    b1   mouse
#> 3      2026-04-05         44    case    b1   mouse
#> 4      2026-04-05        205 control    b1   mouse
#> 5      2026-04-05        163 control    b1   mouse
#> 6      2026-04-05         47 control    b1   mouse
#> 7      2026-04-05        214    case    b2   mouse
#> 8      2026-04-05        158    case    b2   mouse
#> 9      2026-04-05         40    case    b2   mouse
```

## Step 9: Inspect a simple sample summary

``` r
aggregate(expression ~ sample + group, data = merged, FUN = mean)
#>     sample   group expression
#> 1 sample_a    case   138.0000
#> 2 sample_c    case   137.3333
#> 3 sample_b control   138.3333
```

## Step 10: Prepare signature inputs

``` r
geneset_path <- system.file("extdata", "hallmark_demo.gmt", package = "expranno")
read_genesets(geneset_path)
#> $HALLMARK_P53_PATHWAY
#> [1] "TP53" "MDM2" "BAX" 
#> 
#> $HALLMARK_EGFR_SIGNALING
#> [1] "EGFR"  "ERBB2" "GRB2"
```

## Mouse workflow takeaway

The structure is the same as the human case study, but the species must
be fixed to mouse and the annotation backend needs mouse-compatible
resources. Using `"mouse_tpm_v102"` is the clearest way to preserve that
intent. The annotation report and provenance table are the fastest way
to confirm that the species-specific mapping worked and which release
was used.

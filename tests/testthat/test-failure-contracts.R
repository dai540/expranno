test_that("SummarizedExperiment output conversion rejects unsupported inputs", {
  expect_error(
    as_expranno_se(list()),
    "`x` must be an `expranno_result` or `expranno_annotation`"
  )
})

test_that("example_annotation_truth validates species choices", {
  expect_error(
    example_annotation_truth("rat"),
    "should be one of"
  )
})

test_that("run_signature_analysis fails clearly when the chosen gene column is empty", {
  expr_anno <- data.frame(
    gene_id_raw = c("ENSG1.1", "ENSG2.1"),
    gene_id = c("ENSG1", "ENSG2"),
    symbol = c(NA_character_, ""),
    sample_a = c(10, 20),
    sample_b = c(5, 6),
    stringsAsFactors = FALSE
  )

  expect_error(
    run_signature_analysis(
      expr_anno = expr_anno,
      gene_sets = list(pathway_a = c("TP53", "EGFR")),
      method = "gsva",
      expr_scale = "count"
    ),
    "No non-missing values were found in `expr_anno\\$symbol`"
  )
})

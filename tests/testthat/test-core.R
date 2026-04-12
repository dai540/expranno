test_that("example data follows the package contract", {
  demo <- example_expranno_data("human")
  expect_true(validate_expr(demo$expr))
  expect_true(validate_meta(demo$meta))
})

test_that("annotation works in none mode without optional packages", {
  demo <- example_expranno_data("mouse")
  out <- annotate_expr(demo$expr, annotation_engine = "none")

  expect_s3_class(out, "expranno_annotation")
  expect_true(all(c("gene_id_raw", "gene_id", "annotation_status") %in% names(out$expr_anno)))
  expect_equal(out$species, "mouse")
})

test_that("merge_expr_meta creates a long table", {
  demo <- example_expranno_data("human")
  anno <- annotate_expr(demo$expr, annotation_engine = "none")
  merged <- merge_expr_meta(anno$expr_anno, demo$meta)

  expect_true("expression" %in% names(merged))
  expect_equal(length(unique(merged$sample)), 3)
})

test_that("run_expranno handles the minimal workflow", {
  demo <- example_expranno_data("human")
  result <- run_expranno(
    expr = demo$expr,
    meta = demo$meta,
    annotation_engine = "none",
    run_deconvolution = FALSE,
    run_signature = FALSE
  )

  expect_true(is.list(result))
  expect_true("annotation" %in% names(result))
  expect_true("expr_meta_merged" %in% names(result))
})


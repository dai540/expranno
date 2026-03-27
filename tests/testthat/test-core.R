test_that("example data follow package contract", {
  demo <- example_expranno_data()
  mouse_demo <- example_expranno_data("mouse")

  expect_true(is.data.frame(demo$expr))
  expect_true(is.data.frame(demo$meta))
  expect_identical(names(demo$expr)[1], "gene_id")
  expect_identical(names(demo$meta)[1], "sample")
  expect_match(mouse_demo$expr$gene_id[1], "^ENSMUSG")
})

test_that("validation accepts well-formed inputs", {
  demo <- example_expranno_data()

  expect_true(validate_expr(demo$expr))
  expect_true(validate_meta(demo$meta))
})

test_that("annotate_expr with no backend still returns aligned objects", {
  demo <- example_expranno_data()
  result <- annotate_expr(
    expr = demo$expr,
    meta = demo$meta,
    species = "human",
    annotation_engine = "none"
  )

  expect_s3_class(result, "expranno_annotation")
  expect_true("expr_anno" %in% names(result))
  expect_true("meta_checked" %in% names(result))
  expect_equal(nrow(result$expr_anno), nrow(demo$expr))
})

test_that("merge_expr_meta pivots sample columns into long format", {
  demo <- example_expranno_data()
  annotated <- annotate_expr(
    expr = demo$expr,
    meta = demo$meta,
    species = "human",
    annotation_engine = "none"
  )

  merged <- merge_expr_meta(annotated$expr_anno, annotated$meta_checked)

  expect_true(is.data.frame(merged))
  expect_true(all(c("sample", "expression", "group") %in% names(merged)))
  expect_equal(nrow(merged), nrow(demo$expr) * (ncol(demo$expr) - 1L))
})

test_that("read_genesets parses GMT files", {
  path <- tempfile(fileext = ".gmt")
  writeLines(
    c(
      "pathway_a\tdemo\tTP53\tEGFR",
      "pathway_b\tdemo\tBRCA1\tBRCA2"
    ),
    con = path
  )

  sets <- read_genesets(path)
  expect_true(is.list(sets))
  expect_identical(names(sets), c("pathway_a", "pathway_b"))
  expect_identical(sets[[1]], c("TP53", "EGFR"))
})

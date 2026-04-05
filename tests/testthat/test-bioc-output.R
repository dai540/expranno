test_that("annotation-only outputs can be converted to SummarizedExperiment", {
  testthat::skip_if_not_installed("SummarizedExperiment")

  demo <- example_expranno_data()
  annotated <- annotate_expr(
    expr = demo$expr,
    meta = demo$meta,
    species = "human",
    annotation_engine = "none",
    verbose = FALSE
  )

  se <- as_expranno_se(annotated)

  expect_s4_class(se, "SummarizedExperiment")
  expect_identical(rownames(SummarizedExperiment::assay(se)), sub("\\..*$", "", demo$expr$gene_id))
  expect_identical(colnames(SummarizedExperiment::assay(se)), demo$meta$sample)
  expect_true("gene_id_raw" %in% names(as.data.frame(SummarizedExperiment::rowData(se))))
  expect_true("sample" %in% names(as.data.frame(SummarizedExperiment::colData(se))))
  expect_true("expranno" %in% names(S4Vectors::metadata(se)))
})

test_that("Bioconductor coercion fails clearly when row gene identifiers are missing", {
  testthat::skip_if_not_installed("SummarizedExperiment")

  demo <- example_expranno_data()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = as.matrix(demo$expr[, -1, drop = FALSE])),
    rowData = S4Vectors::DataFrame(feature_id = demo$expr$gene_id),
    colData = S4Vectors::DataFrame(sample = demo$meta$sample)
  )

  expect_error(
    as_expranno_input(se, assay_name = "counts", gene_id_col = "gene_id"),
    "Unable to find gene identifiers"
  )
})

test_that("full expranno results carry sample-level outputs into colData", {
  testthat::skip_if_not_installed("SummarizedExperiment")

  demo <- example_expranno_data()
  annotated <- annotate_expr(
    expr = demo$expr,
    meta = demo$meta,
    species = "human",
    annotation_engine = "none",
    verbose = FALSE
  )
  merged <- merge_expr_meta(annotated$expr_anno, annotated$meta_checked)
  result <- new_expranno_result(
    annotation = annotated,
    expr_meta_merged = merged,
    deconvolution = list(
      quantiseq = data.frame(
        cell_type = c("B cell", "T cell"),
        sample_a = c(0.1, 0.2),
        sample_b = c(0.3, 0.4),
        sample_c = c(0.5, 0.6),
        stringsAsFactors = FALSE
      )
    ),
    signatures = list(
      gsva = data.frame(
        signature = c("pathway_a", "pathway_b"),
        sample_a = c(1, 2),
        sample_b = c(3, 4),
        sample_c = c(5, 6),
        stringsAsFactors = FALSE
      )
    ),
    files = character(0),
    benchmark = NULL,
    validation = NULL,
    session_info = NULL
  )

  se <- as_expranno_se(result)
  col_data_names <- names(as.data.frame(SummarizedExperiment::colData(se)))

  expect_true(any(grepl("^signature__", col_data_names)))
  expect_true(any(grepl("^deconv__", col_data_names)))
  expect_identical(colnames(SummarizedExperiment::assay(se)), demo$meta$sample)
})

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
  expect_true("ambiguity_report" %in% names(result))
  expect_true("provenance" %in% names(result))
  expect_equal(nrow(result$expr_anno), nrow(demo$expr))
})

test_that("benchmark_annotation_engines compares engines on the same inputs", {
  demo <- example_expranno_data()

  benchmark <- benchmark_annotation_engines(
    expr = demo$expr,
    meta = demo$meta,
    species = "human",
    engines = c("none"),
    verbose = FALSE
  )

  expect_s3_class(benchmark, "expranno_benchmark")
  expect_true(all(c("engine", "status", "annotated_genes") %in% names(benchmark$summary)))
  expect_identical(benchmark$summary$engine, "none")
})

test_that("annotation presets are listed and can override defaults", {
  presets <- list_annotation_presets()
  mouse_demo <- example_expranno_data("mouse")

  expect_true(all(c("annotation_preset", "species", "biomart_version") %in% names(presets)))
  expect_true("mouse_tpm_v102" %in% presets$annotation_preset)

  annotated <- annotate_expr(
    expr = mouse_demo$expr,
    meta = mouse_demo$meta,
    species = "human",
    annotation_preset = "mouse_tpm_v102",
    annotation_engine = "none",
    verbose = FALSE
  )

  expect_identical(annotated$params$species, "mouse")
  expect_identical(annotated$params$annotation_engine, "hybrid")
  expect_identical(annotated$params$annotation_preset, "mouse_tpm_v102")
  expect_identical(annotated$params$biomart_version, 102)
})

test_that("validate_annotation_engines returns truth-based summary and details", {
  demo <- example_expranno_data()
  truth <- data.frame(
    gene_id = demo$expr$gene_id,
    symbol = c("TP53", "EGFR", "BRCA1"),
    stringsAsFactors = FALSE
  )

  validation <- validate_annotation_engines(
    expr = demo$expr,
    meta = demo$meta,
    truth = truth,
    species = "human",
    engines = c("none"),
    fields = "symbol",
    verbose = FALSE
  )

  expect_s3_class(validation, "expranno_validation")
  expect_true(all(c("engine", "field", "match_rate") %in% names(validation$summary)))
  expect_true(all(c("gene_id", "truth_value", "predicted_value") %in% names(validation$details)))
  expect_identical(validation$summary$engine, "none")
  expect_equal(validation$summary$missing_prediction_rows, 3)
})

test_that("SummarizedExperiment inputs can be coerced to expranno contract", {
  testthat::skip_if_not_installed("SummarizedExperiment")

  demo <- example_expranno_data()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = as.matrix(demo$expr[, -1, drop = FALSE])),
    rowData = S4Vectors::DataFrame(gene_id = demo$expr$gene_id),
    colData = S4Vectors::DataFrame(
      sample = demo$meta$sample,
      group = demo$meta$group,
      batch = demo$meta$batch
    )
  )
  colnames(se) <- demo$meta$sample

  coerced <- as_expranno_input(se, assay_name = "counts")
  annotated <- annotate_expr(
    expr = se,
    species = "human",
    annotation_engine = "none",
    assay_name = "counts",
    verbose = FALSE
  )

  expect_identical(names(coerced$expr)[1], "gene_id")
  expect_identical(names(coerced$meta)[1], "sample")
  expect_s3_class(annotated, "expranno_annotation")
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

test_that("duplicate symbol handling depends on expression scale", {
  expr_anno <- data.frame(
    gene_id_raw = c("ENSG1.1", "ENSG2.1", "ENSG3.1"),
    gene_id = c("ENSG1", "ENSG2", "ENSG3"),
    symbol = c("TP53", "TP53", "EGFR"),
    sample_a = c(10, 20, 5),
    sample_b = c(2, 4, 7),
    stringsAsFactors = FALSE
  )

  count_matrix <- expranno:::collapse_symbol_matrix(
    expr_anno,
    expr_scale = "count",
    duplicate_strategy = "auto"
  )
  abundance_matrix <- expranno:::collapse_symbol_matrix(
    expr_anno,
    expr_scale = "abundance",
    duplicate_strategy = "auto"
  )
  first_matrix <- expranno:::collapse_symbol_matrix(
    expr_anno,
    expr_scale = "count",
    duplicate_strategy = "first"
  )

  expect_equal(count_matrix["TP53", "sample_a"], 30)
  expect_equal(abundance_matrix["TP53", "sample_a"], 15)
  expect_equal(first_matrix["TP53", "sample_a"], 10)
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

test_that("indication-specific deconvolution methods are handled explicitly", {
  expect_identical(
    expranno:::sanitize_deconvolution_methods(
      methods = c("timer", "quantiseq"),
      auto_discovered = TRUE,
      verbose = FALSE
    ),
    "quantiseq"
  )

  expect_error(
    expranno:::sanitize_deconvolution_methods(
      methods = c("timer", "quantiseq"),
      auto_discovered = FALSE,
      verbose = FALSE
    ),
    "require `indications`"
  )

  expect_identical(
    expranno:::sanitize_deconvolution_methods(
      methods = c("timer", "quantiseq"),
      indications = c("SKCM", "SKCM"),
      auto_discovered = FALSE,
      verbose = FALSE
    ),
    c("timer", "quantiseq")
  )
})

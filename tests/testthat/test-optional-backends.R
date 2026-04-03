optional_backend_packages <- function() {
  c(
    "AnnotationDbi",
    "biomaRt",
    "ensembldb",
    "AnnotationFilter",
    "EnsDb.Hsapiens.v86",
    "EnsDb.Mmusculus.v79",
    "GSVA",
    "immunedeconv",
    "org.Hs.eg.db",
    "org.Mm.eg.db"
  )
}

skip_if_optional_backends_missing <- function() {
  for (pkg in optional_backend_packages()) {
    testthat::skip_if_not_installed(pkg)
  }
}

real_like_expr_anno <- function(expr_mat) {
  data.frame(
    gene_id_raw = rownames(expr_mat),
    gene_id = rownames(expr_mat),
    symbol = rownames(expr_mat),
    expr_mat,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

test_that("hybrid annotation smoke test works for human and mouse", {
  skip_if_optional_backends_missing()

  human_demo <- example_expranno_data("human")
  mouse_demo <- example_expranno_data("mouse")

  human_ann <- annotate_expr(
    expr = human_demo$expr,
    meta = human_demo$meta,
    species = "human",
    annotation_engine = "hybrid",
    verbose = FALSE
  )
  mouse_ann <- annotate_expr(
    expr = mouse_demo$expr,
    meta = mouse_demo$meta,
    species = "mouse",
    annotation_engine = "hybrid",
    verbose = FALSE
  )

  expect_true(all(human_ann$report$annotation_rate[human_ann$report$field %in% c("symbol", "gene_name")] > 0))
  expect_true(all(mouse_ann$report$annotation_rate[mouse_ann$report$field %in% c("symbol", "gene_name")] > 0))
  expect_true(any(grepl("orgdb|ensdb|biomart", human_ann$annotation$annotation_source)))
  expect_true(any(grepl("orgdb|ensdb|biomart", mouse_ann$annotation$annotation_source)))
  expect_true("biomart" %in% human_ann$provenance$source)
  expect_true("orgdb" %in% mouse_ann$provenance$source)
  expect_true(any(grepl("Ensembl v102", human_ann$provenance$backend_release, fixed = TRUE)))
  expect_true(all(c("annotation_backend_release", "annotation_date") %in% names(human_ann$annotation)))
})

test_that("annotation benchmark smoke test summarizes hybrid versus single backends", {
  skip_if_optional_backends_missing()

  demo <- example_expranno_data("mouse")
  benchmark <- benchmark_annotation_engines(
    expr = demo$expr,
    meta = demo$meta,
    species = "mouse",
    engines = c("biomart", "orgdb", "hybrid"),
    verbose = FALSE
  )

  expect_s3_class(benchmark, "expranno_benchmark")
  expect_true(all(c("biomart", "orgdb", "hybrid") %in% benchmark$summary$engine))
  expect_true(any(benchmark$coverage$field == "symbol"))
})

test_that("real-like deconvolution smoke tests work for human and mouse", {
  skip_if_optional_backends_missing()

  data("dataset_racle", package = "immunedeconv")
  data("dataset_petitprez", package = "immunedeconv")

  human_expr_anno <- real_like_expr_anno(dataset_racle$expr_mat)
  mouse_expr_anno <- real_like_expr_anno(dataset_petitprez$expr_mat)

  human_res <- run_cell_deconvolution(
    human_expr_anno,
    species = "human",
    methods = c("mcp_counter", "quantiseq"),
    expr_scale = "abundance",
    duplicate_strategy = "first",
    verbose = FALSE
  )
  mouse_res <- run_cell_deconvolution(
    mouse_expr_anno,
    species = "mouse",
    methods = c("mmcp_counter", "base"),
    expr_scale = "count",
    duplicate_strategy = "first",
    verbose = FALSE
  )

  for (nm in names(human_res)) {
    expect_false("error" %in% names(human_res[[nm]]), info = nm)
    expect_true(nrow(human_res[[nm]]) > 0, info = nm)
  }
  for (nm in names(mouse_res)) {
    expect_false("error" %in% names(mouse_res[[nm]]), info = nm)
    expect_true(nrow(mouse_res[[nm]]) > 0, info = nm)
  }
})

test_that("run_expranno smoke test writes annotation, merge, and signature outputs", {
  skip_if_optional_backends_missing()

  demo <- example_expranno_data("human")
  outdir <- tempfile("expranno-smoke-")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  result <- run_expranno(
    expr = demo$expr,
    meta = demo$meta,
    species = "human",
    annotation_engine = "hybrid",
    expr_scale = "count",
    duplicate_strategy = "sum",
    output_dir = outdir,
    run_deconvolution = FALSE,
    run_signature = TRUE,
    geneset_file = system.file("extdata", "hallmark_demo.gmt", package = "expranno"),
    signature_method = "both",
    signature_kcdf = "Poisson",
    verbose = FALSE
  )

  expect_s3_class(result, "expranno_result")
  expect_true(file.exists(file.path(outdir, "expr_anno.csv")))
  expect_true(file.exists(file.path(outdir, "expr_meta_merged.csv")))
  expect_true(file.exists(file.path(outdir, "annotation_report.csv")))
  expect_true(file.exists(file.path(outdir, "annotation_ambiguity.csv")))
  expect_true(file.exists(file.path(outdir, "annotation_provenance.csv")))
  expect_true(file.exists(file.path(outdir, "session_info.txt")))
  expect_true(file.exists(file.path(outdir, "signature_gsva.csv")))
  expect_true(file.exists(file.path(outdir, "signature_ssgsea.csv")))
  expect_true(all(c("gsva", "ssgsea") %in% names(result$signatures)))
})

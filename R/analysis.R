#' Merge annotated expression and metadata
#'
#' @param expr_anno Annotated expression table returned by `annotate_expr()`.
#' @param meta Metadata table with `sample` as the first column.
#' @param output_file Optional CSV path.
#' @param value_name Output column name for expression values.
#'
#' @return A merged `data.frame`.
#' @export
merge_expr_meta <- function(expr_anno, meta, output_file = NULL, value_name = "expression") {
  assert_data_frame(expr_anno, "expr_anno")
  validate_meta(meta)

  sample_cols <- intersect(sample_columns(expr_anno), meta$sample)
  if (length(sample_cols) == 0L) {
    stop("No sample columns in `expr_anno` match `meta$sample`.", call. = FALSE)
  }

  info_cols <- setdiff(names(expr_anno), sample_cols)
  pieces <- lapply(sample_cols, function(sample_name) {
    out <- expr_anno[, info_cols, drop = FALSE]
    out$sample <- sample_name
    out[[value_name]] <- expr_anno[[sample_name]]
    out
  })
  merged <- merge(do.call(rbind, pieces), meta, by = "sample", all.x = TRUE, sort = FALSE)
  rownames(merged) <- NULL
  write_if_requested(merged, output_file)
  merged
}

#' Read gene sets
#'
#' @param x A named list of gene sets or a path to a GMT file.
#'
#' @return A named list of character vectors.
#' @export
read_genesets <- function(x) {
  if (is.list(x)) {
    return(x)
  }
  if (!is.character(x) || length(x) != 1L || !file.exists(x)) {
    stop("`genesets` must be a named list or a valid GMT file path.", call. = FALSE)
  }

  lines <- readLines(x, warn = FALSE)
  out <- lapply(lines, function(line) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    parts[-c(1, 2)]
  })
  names(out) <- vapply(strsplit(lines, "\t", fixed = TRUE), `[`, character(1), 1)
  out
}

#' Run deconvolution on an annotated matrix
#'
#' @param expr_anno Annotated expression table.
#' @param species `"human"` or `"mouse"`.
#' @param methods Character vector of deconvolution methods. If `NULL`, a small
#'   default set is used.
#' @param gene_column Annotation column used as feature names.
#' @param duplicate_strategy `"sum"` or `"mean"`.
#' @param output_dir Optional output directory.
#' @param ... Additional arguments passed to `immunedeconv::deconvolute()`.
#'
#' @return A named list of deconvolution tables.
#' @export
run_cell_deconvolution <- function(
    expr_anno,
    species = c("human", "mouse"),
    methods = NULL,
    gene_column = "symbol",
    duplicate_strategy = c("sum", "mean"),
    output_dir = NULL,
    ...) {
  require_namespace("immunedeconv", reason = "deconvolution")
  species <- match.arg(species)
  duplicate_strategy <- match.arg(duplicate_strategy)
  methods <- methods %||% default_deconvolution_methods(species)
  expr_mat <- build_symbol_matrix(expr_anno, gene_column = gene_column, duplicate_strategy = duplicate_strategy)

  results <- setNames(vector("list", length(methods)), methods)
  for (method in methods) {
    results[[method]] <- immunedeconv::deconvolute(expr_mat, method = method, ...)
    if (!is.null(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      write_if_requested(results[[method]], file.path(output_dir, sprintf("cell_deconv_%s.csv", method)))
    }
  }
  results
}

score_gsva <- function(expr_mat, genesets, method, kcdf) {
  require_namespace("GSVA", reason = "signature scoring")

  if (identical(method, "gsva")) {
    if (exists("gsvaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
      param <- GSVA::gsvaParam(expr_mat, genesets, kcdf = kcdf)
      return(as.data.frame(GSVA::gsva(param, verbose = FALSE), check.names = FALSE))
    }
    return(as.data.frame(GSVA::gsva(expr_mat, genesets, method = "gsva", kcdf = kcdf, verbose = FALSE), check.names = FALSE))
  }

  if (exists("ssgseaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
    param <- GSVA::ssgseaParam(expr_mat, genesets)
    return(as.data.frame(GSVA::gsva(param, verbose = FALSE), check.names = FALSE))
  }
  as.data.frame(GSVA::gsva(expr_mat, genesets, method = "ssgsea", kcdf = kcdf, verbose = FALSE), check.names = FALSE)
}

#' Run GSVA or ssGSEA signature scoring
#'
#' @param expr_anno Annotated expression table.
#' @param genesets Named list of gene sets or a GMT file path.
#' @param method `"gsva"`, `"ssgsea"`, or `"both"`.
#' @param gene_column Annotation column used as feature names.
#' @param duplicate_strategy `"sum"` or `"mean"`.
#' @param kcdf Kernel for GSVA.
#' @param output_dir Optional output directory.
#'
#' @return A named list of signature-score tables.
#' @export
run_signature_analysis <- function(
    expr_anno,
    genesets,
    method = c("gsva", "ssgsea", "both"),
    gene_column = "symbol",
    duplicate_strategy = c("sum", "mean"),
    kcdf = c("Gaussian", "Poisson"),
    output_dir = NULL) {
  method <- match.arg(method)
  duplicate_strategy <- match.arg(duplicate_strategy)
  kcdf <- match.arg(kcdf)

  sets <- read_genesets(genesets)
  expr_mat <- build_symbol_matrix(expr_anno, gene_column = gene_column, duplicate_strategy = duplicate_strategy)
  methods <- if (identical(method, "both")) c("gsva", "ssgsea") else method
  results <- setNames(vector("list", length(methods)), methods)

  for (one_method in methods) {
    score_tbl <- score_gsva(expr_mat, sets, method = one_method, kcdf = kcdf)
    score_tbl$signature <- rownames(score_tbl)
    score_tbl <- score_tbl[, c("signature", colnames(expr_mat)), drop = FALSE]
    rownames(score_tbl) <- NULL
    results[[one_method]] <- score_tbl
    if (!is.null(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      write_if_requested(score_tbl, file.path(output_dir, sprintf("signature_%s.csv", one_method)))
    }
  }

  results
}

#' Run the full expranno workflow
#'
#' @param expr Expression table.
#' @param meta Metadata table.
#' @param species `"auto"`, `"human"`, or `"mouse"`.
#' @param annotation_engine `"hybrid"`, `"biomart"`, `"orgdb"`, or `"none"`.
#' @param strip_version Whether to strip Ensembl version suffixes.
#' @param biomart_version Ensembl release used by `biomaRt`.
#' @param output_dir Optional output directory.
#' @param run_deconvolution Whether to run deconvolution.
#' @param deconv_methods Optional vector of deconvolution methods.
#' @param run_signature Whether to run signature scoring.
#' @param genesets Named list of gene sets or a GMT file path.
#' @param signature_method `"gsva"`, `"ssgsea"`, or `"both"`.
#' @param signature_kcdf `"Gaussian"` or `"Poisson"`.
#' @param duplicate_strategy `"sum"` or `"mean"`.
#' @param ... Additional arguments passed to deconvolution.
#'
#' @return A list with annotation, merged data, and optional analysis outputs.
#' @export
run_expranno <- function(
    expr,
    meta,
    species = c("auto", "human", "mouse"),
    annotation_engine = c("hybrid", "biomart", "orgdb", "none"),
    strip_version = TRUE,
    biomart_version = 102,
    output_dir = NULL,
    run_deconvolution = FALSE,
    deconv_methods = NULL,
    run_signature = FALSE,
    genesets = NULL,
    signature_method = c("gsva", "ssgsea", "both"),
    signature_kcdf = c("Gaussian", "Poisson"),
    duplicate_strategy = c("sum", "mean"),
    ...) {
  validate_meta(meta)
  species <- match.arg(species)
  annotation_engine <- match.arg(annotation_engine)
  signature_method <- match.arg(signature_method)
  signature_kcdf <- match.arg(signature_kcdf)
  duplicate_strategy <- match.arg(duplicate_strategy)

  annotation <- annotate_expr(
    expr = expr,
    species = species,
    annotation_engine = annotation_engine,
    strip_version = strip_version,
    biomart_version = biomart_version,
    output_dir = output_dir
  )

  merged <- merge_expr_meta(
    expr_anno = annotation$expr_anno,
    meta = meta,
    output_file = if (is.null(output_dir)) NULL else file.path(output_dir, "expr_meta_merged.csv")
  )

  deconvolution <- NULL
  if (isTRUE(run_deconvolution)) {
    deconvolution <- run_cell_deconvolution(
      expr_anno = annotation$expr_anno,
      species = annotation$species,
      methods = deconv_methods,
      duplicate_strategy = duplicate_strategy,
      output_dir = output_dir,
      ...
    )
  }

  signatures <- NULL
  if (isTRUE(run_signature)) {
    if (is.null(genesets)) {
      stop("`genesets` must be supplied when `run_signature = TRUE`.", call. = FALSE)
    }
    signatures <- run_signature_analysis(
      expr_anno = annotation$expr_anno,
      genesets = genesets,
      method = signature_method,
      duplicate_strategy = duplicate_strategy,
      kcdf = signature_kcdf,
      output_dir = output_dir
    )
  }

  list(
    annotation = annotation,
    expr_meta_merged = merged,
    deconvolution = deconvolution,
    signatures = signatures
  )
}

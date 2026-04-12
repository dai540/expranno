`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

require_namespace <- function(pkg, reason) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for %s.", pkg, reason), call. = FALSE)
  }
}

assert_data_frame <- function(x, name) {
  if (!is.data.frame(x)) {
    stop(sprintf("`%s` must be a data.frame.", name), call. = FALSE)
  }
}

assert_first_column <- function(x, expected, name) {
  if (!identical(names(x)[1], expected)) {
    stop(sprintf("The first column of `%s` must be `%s`.", name, expected), call. = FALSE)
  }
}

write_if_requested <- function(x, path) {
  if (!is.null(path)) {
    utils::write.csv(x, path, row.names = FALSE, na = "")
  }
  invisible(x)
}

detect_species <- function(gene_ids) {
  gene_ids <- as.character(gene_ids)
  keep <- !is.na(gene_ids) & nzchar(gene_ids)
  if (all(grepl("^ENSG", gene_ids[keep]))) {
    return("human")
  }
  if (all(grepl("^ENSMUSG", gene_ids[keep]))) {
    return("mouse")
  }
  stop("Unable to infer species from `gene_id`. Set `species` explicitly.", call. = FALSE)
}

sample_columns <- function(x) {
  names(x)[vapply(x, is.numeric, logical(1))]
}

build_symbol_matrix <- function(expr_anno, gene_column = "symbol", duplicate_strategy = c("sum", "mean")) {
  duplicate_strategy <- match.arg(duplicate_strategy)
  if (!gene_column %in% names(expr_anno)) {
    stop(sprintf("`%s` is not present in `expr_anno`.", gene_column), call. = FALSE)
  }

  sample_cols <- sample_columns(expr_anno)
  genes <- as.character(expr_anno[[gene_column]])
  keep <- !is.na(genes) & nzchar(genes)
  if (!any(keep)) {
    stop(sprintf("`expr_anno$%s` does not contain usable gene symbols.", gene_column), call. = FALSE)
  }

  tbl <- data.frame(
    gene = genes[keep],
    expr_anno[keep, sample_cols, drop = FALSE],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  fun <- if (identical(duplicate_strategy, "sum")) sum else mean
  collapsed <- stats::aggregate(. ~ gene, data = tbl, FUN = fun)
  mat <- as.matrix(collapsed[, -1, drop = FALSE])
  rownames(mat) <- collapsed$gene
  mode(mat) <- "numeric"
  mat
}

#' Validate an expression table
#'
#' @param expr A data frame whose first column is `gene_id` and remaining
#'   columns are numeric sample columns.
#'
#' @return `TRUE` invisibly.
#' @export
validate_expr <- function(expr) {
  assert_data_frame(expr, "expr")
  assert_first_column(expr, "gene_id", "expr")

  if (ncol(expr) < 2L) {
    stop("`expr` must contain at least one sample column after `gene_id`.", call. = FALSE)
  }
  if (anyDuplicated(names(expr)[-1]) > 0L) {
    stop("Sample column names in `expr` must be unique.", call. = FALSE)
  }
  if (!all(vapply(expr[-1], is.numeric, logical(1)))) {
    stop("All sample columns in `expr` must be numeric.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate a metadata table
#'
#' @param meta A data frame whose first column is `sample`.
#'
#' @return `TRUE` invisibly.
#' @export
validate_meta <- function(meta) {
  assert_data_frame(meta, "meta")
  assert_first_column(meta, "sample", "meta")

  if (anyDuplicated(meta[[1]]) > 0L) {
    stop("`meta$sample` must contain unique sample IDs.", call. = FALSE)
  }
  invisible(TRUE)
}

default_deconvolution_methods <- function(species) {
  if (identical(species, "human")) {
    c("quantiseq", "mcp_counter", "xcell")
  } else {
    c("mmcp_counter")
  }
}

#' Small built-in example inputs
#'
#' @param species Either `"human"` or `"mouse"`.
#'
#' @return A list with `expr` and `meta`.
#' @export
example_expranno_data <- function(species = c("human", "mouse")) {
  species <- match.arg(species)

  expr <- switch(
    species,
    human = data.frame(
      gene_id = c("ENSG00000141510.17", "ENSG00000146648.18", "ENSG00000012048.23"),
      sample_a = c(120, 80, 25),
      sample_b = c(140, 77, 30),
      sample_c = c(118, 91, 21),
      stringsAsFactors = FALSE
    ),
    mouse = data.frame(
      gene_id = c("ENSMUSG00000059552.8", "ENSMUSG00000020122.15", "ENSMUSG00000017167.16"),
      sample_a = c(220, 150, 44),
      sample_b = c(205, 163, 47),
      sample_c = c(214, 158, 40),
      stringsAsFactors = FALSE
    )
  )

  meta <- data.frame(
    sample = c("sample_a", "sample_b", "sample_c"),
    group = c("case", "control", "case"),
    batch = c("b1", "b1", "b2"),
    stringsAsFactors = FALSE
  )

  list(expr = expr, meta = meta)
}


#' Validate an expression table
#'
#' `expr` must be a `data.frame` with `gene_id` as the first column and sample
#' columns after it.
#'
#' @param expr A gene-by-sample expression table.
#'
#' @return `TRUE` invisibly on success.
#' @export
validate_expr <- function(expr) {
  assert_data_frame(expr, "expr")
  assert_first_col(expr, "gene_id", "expr")

  if (ncol(expr) < 2L) {
    stop("`expr` must contain at least one sample column after `gene_id`.", call. = FALSE)
  }
  if (anyDuplicated(names(expr)[-1]) > 0L) {
    stop("Sample column names in `expr` must be unique.", call. = FALSE)
  }
  numeric_ok <- vapply(expr[-1], is.numeric, logical(1))
  if (!all(numeric_ok)) {
    stop("All sample columns in `expr` must be numeric.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate a metadata table
#'
#' `meta` must be a `data.frame` with `sample` as the first column and one row
#' per sample.
#'
#' @param meta A sample metadata table.
#'
#' @return `TRUE` invisibly on success.
#' @export
validate_meta <- function(meta) {
  assert_data_frame(meta, "meta")
  assert_first_col(meta, "sample", "meta")

  if (anyDuplicated(meta[[1]]) > 0L) {
    stop("`meta$sample` must contain unique sample IDs.", call. = FALSE)
  }
  invisible(TRUE)
}

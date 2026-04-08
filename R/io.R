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

#' Merge annotated expression and metadata
#'
#' Converts the annotated expression matrix to long format, joins it to
#' `meta` by `sample`, and optionally writes the merged table as CSV.
#'
#' @param expr_anno Annotated expression table returned by [annotate_expr()].
#' @param meta Checked or user-supplied metadata table with `sample` first.
#' @param output_file Optional CSV path.
#' @param value_name Name of the expression value column in the long table.
#'
#' @return A merged `data.frame`.
#' @export
merge_expr_meta <- function(expr_anno, meta, output_file = NULL, value_name = "expression") {
  assert_data_frame(expr_anno, "expr_anno")
  validate_meta(meta)

  sample_cols <- intersect(names(expr_anno), as.character(meta[[1]]))
  if (length(sample_cols) == 0L) {
    stop("No sample columns in `expr_anno` match `meta$sample`.", call. = FALSE)
  }

  annotation_cols <- setdiff(names(expr_anno), sample_cols)
  pieces <- lapply(sample_cols, function(sample_col) {
    out <- expr_anno[, annotation_cols, drop = FALSE]
    out$sample <- sample_col
    out[[value_name]] <- expr_anno[[sample_col]]
    out
  })
  long_expr <- do.call(rbind, pieces)
  rownames(long_expr) <- NULL

  merged <- merge(long_expr, meta, by = "sample", all.x = TRUE, sort = FALSE)
  write_if_requested(merged, output_file = output_file)
  merged
}

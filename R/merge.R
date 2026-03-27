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

#' @export
print.expranno_annotation <- function(x, ...) {
  cat("<expranno_annotation>\n")
  cat(sprintf("  genes: %s\n", nrow(x$expr_anno)))
  cat(sprintf("  samples: %s\n", nrow(x$meta_checked)))
  cat(sprintf("  species: %s\n", x$params$species))
  invisible(x)
}

#' @export
summary.expranno_annotation <- function(object, ...) {
  object$report
}

#' @export
print.expranno_result <- function(x, ...) {
  cat("<expranno_result>\n")
  cat(sprintf("  annotated genes: %s\n", nrow(x$annotation$expr_anno)))
  cat(sprintf("  merged rows: %s\n", nrow(x$expr_meta_merged)))
  cat(sprintf("  deconvolution runs: %s\n", length(x$deconvolution %||% list())))
  cat(sprintf("  signature runs: %s\n", length(x$signatures %||% list())))
  invisible(x)
}

#' @export
print.expranno_annotation <- function(x, ...) {
  cat("<expranno_annotation>\n")
  cat(sprintf("  genes: %s\n", nrow(x$expr_anno)))
  cat(sprintf("  samples: %s\n", nrow(x$meta_checked)))
  cat(sprintf("  species: %s\n", x$params$species))
  cat(sprintf("  annotation engine: %s\n", x$params$annotation_engine))
  if (!is.null(x$params$annotation_preset)) {
    cat(sprintf("  annotation preset: %s\n", x$params$annotation_preset))
  }
  cat(sprintf("  ambiguous gene-fields: %s\n", nrow(x$ambiguity_report)))
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
  cat(sprintf("  benchmark runs: %s\n", if (is.null(x$benchmark)) 0 else nrow(x$benchmark$summary)))
  cat(sprintf("  validation runs: %s\n", if (is.null(x$validation)) 0 else nrow(x$validation$summary)))
  invisible(x)
}

#' @export
print.expranno_benchmark <- function(x, ...) {
  cat("<expranno_benchmark>\n")
  cat(sprintf("  engines: %s\n", paste(x$params$engines, collapse = ", ")))
  cat(sprintf("  fields: %s\n", paste(x$params$fields, collapse = ", ")))
  cat(sprintf("  successful runs: %s\n", sum(x$summary$status == "ok")))
  invisible(x)
}

#' @export
print.expranno_validation <- function(x, ...) {
  cat("<expranno_validation>\n")
  cat(sprintf("  engines: %s\n", paste(x$params$engines, collapse = ", ")))
  cat(sprintf("  fields: %s\n", paste(x$params$fields, collapse = ", ")))
  cat(sprintf("  compared rows: %s\n", nrow(x$details)))
  invisible(x)
}

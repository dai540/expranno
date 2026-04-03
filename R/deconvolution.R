detect_deconvolution_methods <- function() {
  if (!requireNamespace("immunedeconv", quietly = TRUE)) {
    return(character(0))
  }

  candidates <- c(
    "deconvolution_methods",
    "deconvolution_methods_mouse"
  )
  found <- character(0)
  for (fn in candidates) {
    if (!exists(fn, where = asNamespace("immunedeconv"), inherits = FALSE)) {
      next
    }
    values <- tryCatch(
      {
        obj <- getExportedValue("immunedeconv", fn)
        if (is.function(obj)) obj() else obj
      },
      error = function(e) character(0)
    )
    found <- c(found, values)
  }

  found <- unique(found)
  found[!grepl("cibersort", found, ignore.case = TRUE)]
}

indication_required_methods <- function() {
  c("timer", "consensus_tme")
}

sanitize_deconvolution_methods <- function(methods, indications = NULL, auto_discovered = FALSE, verbose = TRUE) {
  required <- indication_required_methods()
  needs_indications <- methods[tolower(methods) %in% required]
  has_indications <- !is.null(indications)

  if (length(needs_indications) == 0L || has_indications) {
    return(methods)
  }

  if (!auto_discovered) {
    stop(
      paste0(
        "The following methods require `indications`: ",
        paste(needs_indications, collapse = ", "),
        ". Supply them through `run_cell_deconvolution(..., indications = ...)` ",
        "or `run_expranno(..., deconv_args = list(indications = ...))`."
      ),
      call. = FALSE
    )
  }

  kept <- methods[!(tolower(methods) %in% required)]
  if (isTRUE(verbose)) {
    message(
      sprintf(
        "Skipping indication-specific methods without `indications`: %s",
        paste(needs_indications, collapse = ", ")
      )
    )
  }
  kept
}

run_one_deconvolution <- function(matrix_input, method, species, ...) {
  args <- c(list(matrix_input), list(method = method), list(...))
  if (species == "mouse" && exists("deconvolute_mouse", where = asNamespace("immunedeconv"), inherits = FALSE)) {
    return(do.call(getExportedValue("immunedeconv", "deconvolute_mouse"), args))
  }
  do.call(getExportedValue("immunedeconv", "deconvolute"), args)
}

#' Run immune deconvolution
#'
#' Uses `immunedeconv` on an annotated expression matrix. By default, all
#' available methods except CIBERSORT-family methods are executed.
#'
#' @param expr_anno Annotated expression table.
#' @param species Either `"auto"`, `"human"`, or `"mouse"`.
#' @param methods Either `"all_except_cibersort"` or a character vector of
#'   explicit method names.
#' @param gene_column Gene symbol column to use for deconvolution.
#' @param expr_scale Expression scale. This affects duplicate symbol handling
#'   and triggers a warning for `expr_scale = "log"`.
#' @param duplicate_strategy Strategy used when multiple rows map to the same
#'   symbol. `"auto"` uses `"sum"` for counts and `"mean"` otherwise.
#' @param output_dir Optional directory to write one CSV per method.
#' @param prefix File name prefix for output CSV files.
#' @param verbose Whether to emit progress messages.
#' @param ... Additional arguments forwarded to `immunedeconv`, including
#'   method-specific values such as `indications`.
#'
#' @return A named list of deconvolution result tables.
#' @export
run_cell_deconvolution <- function(
    expr_anno,
    species = c("auto", "human", "mouse"),
    methods = "all_except_cibersort",
    gene_column = "symbol",
    expr_scale = c("auto", "count", "abundance", "log"),
    duplicate_strategy = c("auto", "sum", "mean", "max", "first"),
    output_dir = NULL,
    prefix = "cell_deconv_",
    verbose = TRUE,
    ...) {
  require_namespace("immunedeconv", reason = "immune deconvolution")

  species <- match.arg(species)
  if (species == "auto") {
    species <- resolve_species("auto", expr_anno$gene_id)
  }
  expr_scale <- match.arg(expr_scale)
  duplicate_strategy <- match.arg(duplicate_strategy)

  if (identical(expr_scale, "log")) {
    warning(
      "Deconvolution methods generally expect non-log expression values. ",
      "Consider providing counts or abundance-like values instead.",
      call. = FALSE
    )
  }

  matrix_input <- collapse_symbol_matrix(
    expr_anno,
    gene_column = gene_column,
    expr_scale = expr_scale,
    duplicate_strategy = duplicate_strategy
  )

  extra_args <- list(...)
  auto_discovered <- is.character(methods) && length(methods) == 1L && methods == "all_except_cibersort"

  if (auto_discovered) {
    methods <- detect_deconvolution_methods()
  }
  methods <- sanitize_deconvolution_methods(
    methods = methods,
    indications = extra_args$indications,
    auto_discovered = auto_discovered,
    verbose = verbose
  )
  if (length(methods) == 0L) {
    stop("No immunedeconv methods were discovered.", call. = FALSE)
  }

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  results <- list()
  for (method in methods) {
    if (verbose) {
      message(sprintf("Running deconvolution method: %s", method))
    }
    out <- tryCatch(
      do.call(
        run_one_deconvolution,
        c(
          list(matrix_input = matrix_input, method = method, species = species),
          extra_args
        )
      ),
      error = function(e) {
        data.frame(
          method = method,
          error = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
    results[[method]] <- as.data.frame(out)
    if (!is.null(output_dir)) {
      utils::write.csv(
        results[[method]],
        file = file.path(output_dir, sprintf("%s%s.csv", prefix, method)),
        row.names = FALSE
      )
    }
  }

  results
}

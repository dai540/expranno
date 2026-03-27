detect_deconvolution_methods <- function() {
  if (!requireNamespace("immunedeconv", quietly = TRUE)) {
    return(character(0))
  }

  candidates <- c(
    "deconvolution_methods",
    "deconvolution_methods_human",
    "deconvolution_methods_mouse"
  )
  found <- character(0)
  for (fn in candidates) {
    if (!exists(fn, where = asNamespace("immunedeconv"), inherits = FALSE)) {
      next
    }
    values <- tryCatch(getExportedValue("immunedeconv", fn)(), error = function(e) character(0))
    found <- c(found, values)
  }

  found <- unique(found)
  found[!grepl("cibersort", found, ignore.case = TRUE)]
}

run_one_deconvolution <- function(matrix_input, method, species) {
  if (species == "mouse" && exists("deconvolute_mouse", where = asNamespace("immunedeconv"), inherits = FALSE)) {
    return(immunedeconv::deconvolute_mouse(matrix_input, method = method))
  }
  immunedeconv::deconvolute(matrix_input, method = method)
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
#' @param output_dir Optional directory to write one CSV per method.
#' @param prefix File name prefix for output CSV files.
#' @param verbose Whether to emit progress messages.
#'
#' @return A named list of deconvolution result tables.
#' @export
run_cell_deconvolution <- function(
    expr_anno,
    species = c("auto", "human", "mouse"),
    methods = "all_except_cibersort",
    gene_column = "symbol",
    output_dir = NULL,
    prefix = "cell_deconv_",
    verbose = TRUE) {
  require_namespace("immunedeconv", reason = "immune deconvolution")

  species <- match.arg(species)
  if (species == "auto") {
    species <- resolve_species("auto", expr_anno$gene_id)
  }

  matrix_input <- collapse_symbol_matrix(expr_anno, gene_column = gene_column)

  if (is.character(methods) && length(methods) == 1L && methods == "all_except_cibersort") {
    methods <- detect_deconvolution_methods()
  }
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
      run_one_deconvolution(matrix_input, method = method, species = species),
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

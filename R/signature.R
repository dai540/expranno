#' Read gene sets from a GMT file
#'
#' @param geneset_file A `.gmt` file path.
#'
#' @return A named list of character vectors.
#' @export
read_genesets <- function(geneset_file) {
  lines <- readLines(geneset_file, warn = FALSE)
  parsed <- strsplit(lines, "\t", fixed = TRUE)
  sets <- lapply(parsed, function(x) x[-c(1, 2)])
  names(sets) <- vapply(parsed, `[`, character(1), 1)
  sets
}

run_gsva_method <- function(expr_matrix, gene_sets, method) {
  require_namespace("GSVA", reason = "signature scoring")
  GSVA::gsva(
    expr = expr_matrix,
    gset.idx.list = gene_sets,
    method = method,
    kcdf = "Gaussian",
    abs.ranking = FALSE,
    verbose = FALSE
  )
}

#' Run GSVA or ssGSEA signature scoring
#'
#' Uses a symbol-based expression matrix derived from `expr_anno` and a
#' user-supplied gene set database.
#'
#' @param expr_anno Annotated expression table.
#' @param geneset_file Optional GMT path.
#' @param gene_sets Optional named list of gene sets.
#' @param method One of `"gsva"`, `"ssgsea"`, or `"both"`.
#' @param gene_column Symbol column to use.
#' @param output_dir Optional output directory.
#' @param prefix File prefix.
#'
#' @return A named list with score matrices as data frames.
#' @export
run_signature_analysis <- function(
    expr_anno,
    geneset_file = NULL,
    gene_sets = NULL,
    method = c("gsva", "ssgsea", "both"),
    gene_column = "symbol",
    output_dir = NULL,
    prefix = "signature_") {
  method <- match.arg(method)
  if (is.null(gene_sets)) {
    if (is.null(geneset_file)) {
      stop("Provide either `geneset_file` or `gene_sets`.", call. = FALSE)
    }
    gene_sets <- read_genesets(geneset_file)
  }
  expr_matrix <- collapse_symbol_matrix(expr_anno, gene_column = gene_column)

  methods <- if (method == "both") c("gsva", "ssgsea") else method
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  results <- list()
  for (one_method in methods) {
    scores <- run_gsva_method(expr_matrix = expr_matrix, gene_sets = gene_sets, method = one_method)
    results[[one_method]] <- as.data.frame(scores)
    results[[one_method]]$signature <- rownames(scores)
    results[[one_method]] <- results[[one_method]][, c("signature", colnames(scores)), drop = FALSE]
    if (!is.null(output_dir)) {
      utils::write.csv(
        results[[one_method]],
        file = file.path(output_dir, sprintf("%s%s.csv", prefix, one_method)),
        row.names = FALSE
      )
    }
  }

  results
}

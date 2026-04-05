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

has_gsva_param_api <- function() {
  all(
    c("gsvaParam", "ssgseaParam", "gsva") %in%
      getNamespaceExports("GSVA")
  )
}

build_gsva_param <- function(
    expr_matrix,
    gene_sets,
    method,
    kcdf,
    min_size,
    max_size,
    gsva_args,
    ssgsea_args) {
  if (identical(method, "gsva")) {
    param_args <- utils::modifyList(
      list(
        exprData = expr_matrix,
        geneSets = gene_sets,
        minSize = min_size,
        maxSize = max_size,
        kcdf = kcdf,
        tau = 1,
        maxDiff = TRUE,
        absRanking = FALSE,
        sparse = TRUE,
        checkNA = "auto",
        use = "everything",
        verbose = FALSE
      ),
      gsva_args
    )
    return(do.call(getExportedValue("GSVA", "gsvaParam"), param_args))
  }

  param_args <- utils::modifyList(
    list(
      exprData = expr_matrix,
      geneSets = gene_sets,
      minSize = min_size,
      maxSize = max_size,
      alpha = 0.25,
      normalize = TRUE,
      checkNA = "auto",
      use = "everything"
    ),
    ssgsea_args
  )
  do.call(getExportedValue("GSVA", "ssgseaParam"), param_args)
}

run_gsva_method <- function(
    expr_matrix,
    gene_sets,
    method,
    kcdf = c("auto", "Gaussian", "Poisson", "none"),
    min_size = 1,
    max_size = Inf,
    gsva_args = list(),
    ssgsea_args = list()) {
  require_namespace("GSVA", reason = "signature scoring")
  kcdf <- match.arg(kcdf)

  if (has_gsva_param_api()) {
    param <- build_gsva_param(
      expr_matrix = expr_matrix,
      gene_sets = gene_sets,
      method = method,
      kcdf = kcdf,
      min_size = min_size,
      max_size = max_size,
      gsva_args = gsva_args,
      ssgsea_args = ssgsea_args
    )
    return(getExportedValue("GSVA", "gsva")(param, verbose = FALSE))
  }

  warning(
    "Falling back to the legacy GSVA API because parameter-object constructors ",
    "were not found in the installed GSVA version.",
    call. = FALSE
  )
  GSVA::gsva(
    expr = expr_matrix,
    gset.idx.list = gene_sets,
    method = method,
    kcdf = kcdf,
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
#' @param expr_scale Expression scale. This affects duplicate symbol handling.
#' @param duplicate_strategy Strategy used when multiple rows map to the same
#'   symbol. `"auto"` uses `"sum"` for counts and `"mean"` otherwise.
#' @param kcdf Kernel choice passed to `GSVA::gsvaParam()`.
#' @param min_size Minimum gene-set size after ID matching.
#' @param max_size Maximum gene-set size after ID matching.
#' @param gsva_args Optional named list of extra arguments passed to
#'   `GSVA::gsvaParam()`.
#' @param ssgsea_args Optional named list of extra arguments passed to
#'   `GSVA::ssgseaParam()`.
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
    expr_scale = c("auto", "count", "abundance", "log"),
    duplicate_strategy = c("auto", "sum", "mean", "max", "first"),
    kcdf = c("auto", "Gaussian", "Poisson", "none"),
    min_size = 1,
    max_size = Inf,
    gsva_args = list(),
    ssgsea_args = list(),
    output_dir = NULL,
    prefix = "signature_") {
  method <- match.arg(method)
  expr_scale <- match.arg(expr_scale)
  duplicate_strategy <- match.arg(duplicate_strategy)
  kcdf <- match.arg(kcdf)
  if (is.null(gene_sets)) {
    if (is.null(geneset_file)) {
      stop("Provide either `geneset_file` or `gene_sets`.", call. = FALSE)
    }
    gene_sets <- read_genesets(geneset_file)
  }
  expr_matrix <- collapse_symbol_matrix(
    expr_anno,
    gene_column = gene_column,
    expr_scale = expr_scale,
    duplicate_strategy = duplicate_strategy
  )
  if (nrow(expr_matrix) == 0L) {
    stop(
      sprintf(
        "No non-missing values were found in `expr_anno$%s` after annotation filtering. Run annotation with a backend that populates `%s`, or choose a populated `gene_column`.",
        gene_column,
        gene_column
      ),
      call. = FALSE
    )
  }

  methods <- if (method == "both") c("gsva", "ssgsea") else method
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  results <- list()
  for (one_method in methods) {
    scores <- run_gsva_method(
      expr_matrix = expr_matrix,
      gene_sets = gene_sets,
      method = one_method,
      kcdf = kcdf,
      min_size = min_size,
      max_size = max_size,
      gsva_args = gsva_args,
      ssgsea_args = ssgsea_args
    )
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

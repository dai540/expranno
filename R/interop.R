#' Coerce Bioconductor containers to expranno inputs
#'
#' Converts a `SummarizedExperiment`-like object into the strict `expr` and
#' `meta` data-frame contract used by `expranno`.
#'
#' @param x A `SummarizedExperiment` or `SingleCellExperiment`.
#' @param assay_name Optional assay name or index. Defaults to the first assay.
#' @param gene_id_col Optional row-data column containing Ensembl gene IDs.
#' @param sample_col Column name to use for sample identifiers in the output
#'   metadata table.
#'
#' @return A named list with `expr` and `meta`.
#' @export
as_expranno_input <- function(
    x,
    assay_name = NULL,
    gene_id_col = NULL,
    sample_col = "sample") {
  require_namespace("SummarizedExperiment", reason = "SummarizedExperiment input support")

  if (!inherits(x, "SummarizedExperiment")) {
    stop("`x` must inherit from SummarizedExperiment.", call. = FALSE)
  }

  assay_value <- if (is.null(assay_name)) {
    SummarizedExperiment::assay(x)
  } else {
    SummarizedExperiment::assay(x, assay_name)
  }

  row_data <- as.data.frame(SummarizedExperiment::rowData(x), stringsAsFactors = FALSE)
  if (is.null(gene_id_col)) {
    candidates <- c("gene_id", "ensembl_gene_id", "ENSEMBL", "Geneid")
    gene_id_col <- candidates[candidates %in% names(row_data)][1] %||% NA_character_
  }

  if (!is.na(gene_id_col) && gene_id_col %in% names(row_data)) {
    gene_ids <- row_data[[gene_id_col]]
  } else {
    gene_ids <- rownames(assay_value)
  }

  if (is.null(gene_ids) || all(is.na(gene_ids) | !nzchar(as.character(gene_ids)))) {
    stop(
      "Unable to find gene identifiers in rowData or row names. Set `gene_id_col` explicitly.",
      call. = FALSE
    )
  }

  expr <- data.frame(
    gene_id = as.character(gene_ids),
    as.data.frame(assay_value, check.names = FALSE, stringsAsFactors = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  meta <- as.data.frame(SummarizedExperiment::colData(x), stringsAsFactors = FALSE)
  if (sample_col %in% names(meta)) {
    meta[[sample_col]] <- as.character(meta[[sample_col]])
  } else {
    meta[[sample_col]] <- colnames(assay_value)
  }

  meta <- meta[, c(sample_col, setdiff(names(meta), sample_col)), drop = FALSE]
  names(meta)[1] <- "sample"

  list(expr = expr, meta = meta)
}

coerce_expranno_inputs <- function(
    expr,
    meta = NULL,
    assay_name = NULL,
    gene_id_col = NULL,
    sample_col = "sample") {
  if (is.data.frame(expr)) {
    if (is.null(meta)) {
      stop("`meta` must be supplied when `expr` is a data.frame.", call. = FALSE)
    }
    return(list(expr = expr, meta = meta, input_type = "data.frame"))
  }

  coerced <- as_expranno_input(
    x = expr,
    assay_name = assay_name,
    gene_id_col = gene_id_col,
    sample_col = sample_col
  )
  coerced$input_type <- class(expr)[1]
  coerced
}

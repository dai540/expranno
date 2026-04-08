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

sanitize_output_feature_names <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  out <- gsub("^_+|_+$", "", out)
  out[out == ""] <- "feature"
  make.unique(out, sep = "_")
}

flatten_sample_level_table <- function(tbl, method, prefix, sample_order) {
  if (is.null(tbl) || !is.data.frame(tbl) || ncol(tbl) < 2L || "error" %in% names(tbl)) {
    return(NULL)
  }

  feature_col <- names(tbl)[1]
  sample_cols <- setdiff(names(tbl), feature_col)
  if (length(sample_cols) == 0L) {
    return(NULL)
  }

  mat <- as.matrix(tbl[, sample_cols, drop = FALSE])
  mode(mat) <- "numeric"
  rownames(mat) <- sanitize_output_feature_names(tbl[[feature_col]])

  out <- as.data.frame(t(mat), check.names = FALSE, stringsAsFactors = FALSE)
  out$sample <- rownames(out)
  rownames(out) <- NULL
  out <- out[match(sample_order, out$sample), , drop = FALSE]
  names(out)[names(out) != "sample"] <- paste(
    prefix,
    sanitize_output_feature_names(method),
    sanitize_output_feature_names(names(out)[names(out) != "sample"]),
    sep = "__"
  )
  out
}

bind_sample_level_outputs <- function(outputs, prefix, sample_order) {
  if (is.null(outputs) || length(outputs) == 0L) {
    return(NULL)
  }

  flattened <- mapply(
    FUN = flatten_sample_level_table,
    tbl = outputs,
    method = names(outputs),
    MoreArgs = list(prefix = prefix, sample_order = sample_order),
    SIMPLIFY = FALSE,
    USE.NAMES = TRUE
  )
  flattened <- Filter(Negate(is.null), flattened)
  if (length(flattened) == 0L) {
    return(NULL)
  }

  merged <- Reduce(function(x, y) merge(x, y, by = "sample", all = TRUE, sort = FALSE), flattened)
  merged[match(sample_order, merged$sample), , drop = FALSE]
}

#' Convert expranno output to SummarizedExperiment
#'
#' Creates a `SummarizedExperiment` with annotated expression values in the main
#' assay, annotation fields in `rowData`, and sample metadata in `colData`.
#' When available, signature and deconvolution outputs are also appended to
#' `colData` as sample-level columns.
#'
#' @param x An `expranno_result` or `expranno_annotation`.
#' @param assay_name Assay name used for the expression matrix.
#' @param include_signatures Whether to append signature scores to `colData`
#'   when `x` is an `expranno_result`.
#' @param include_deconvolution Whether to append deconvolution scores to
#'   `colData` when `x` is an `expranno_result`.
#'
#' @return A `SummarizedExperiment`.
#' @export
as_expranno_se <- function(
    x,
    assay_name = "expression",
    include_signatures = TRUE,
    include_deconvolution = TRUE) {
  if (inherits(x, "expranno_result")) {
    annotation <- x$annotation
    signatures <- if (isTRUE(include_signatures)) x$signatures else NULL
    deconvolution <- if (isTRUE(include_deconvolution)) x$deconvolution else NULL
    metadata_payload <- list(
      files = x$files,
      benchmark = x$benchmark,
      validation = x$validation,
      session_info = x$session_info
    )
  } else if (inherits(x, "expranno_annotation")) {
    annotation <- x
    signatures <- NULL
    deconvolution <- NULL
    metadata_payload <- list()
  } else {
    stop("`x` must be an `expranno_result` or `expranno_annotation`.", call. = FALSE)
  }

  require_namespace("SummarizedExperiment", reason = "Bioconductor output support")
  require_namespace("S4Vectors", reason = "Bioconductor output support")

  expr_anno <- annotation$expr_anno
  meta <- annotation$meta_checked
  sample_order <- as.character(meta$sample)
  sample_cols <- intersect(names(expr_anno), sample_order)
  if (length(sample_cols) == 0L) {
    stop("No sample columns were found in `expr_anno`.", call. = FALSE)
  }

  expr_matrix <- as.matrix(expr_anno[, sample_cols, drop = FALSE])
  mode(expr_matrix) <- "numeric"
  rownames(expr_matrix) <- expr_anno$gene_id
  colnames(expr_matrix) <- sample_cols

  row_data <- expr_anno[, setdiff(names(expr_anno), sample_cols), drop = FALSE]
  col_data <- meta[match(sample_cols, sample_order), , drop = FALSE]

  signature_cols <- bind_sample_level_outputs(signatures, prefix = "signature", sample_order = sample_cols)
  if (!is.null(signature_cols)) {
    col_data <- merge(col_data, signature_cols, by = "sample", all.x = TRUE, sort = FALSE)
    col_data <- col_data[match(sample_cols, col_data$sample), , drop = FALSE]
  }

  deconv_cols <- bind_sample_level_outputs(deconvolution, prefix = "deconv", sample_order = sample_cols)
  if (!is.null(deconv_cols)) {
    col_data <- merge(col_data, deconv_cols, by = "sample", all.x = TRUE, sort = FALSE)
    col_data <- col_data[match(sample_cols, col_data$sample), , drop = FALSE]
  }

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = stats::setNames(list(expr_matrix), assay_name),
    rowData = S4Vectors::DataFrame(row_data, check.names = FALSE),
    colData = S4Vectors::DataFrame(col_data, check.names = FALSE)
  )
  S4Vectors::metadata(se)$expranno <- utils::modifyList(
    list(
      provenance = annotation$provenance,
      report = annotation$report,
      ambiguity_report = annotation$ambiguity_report,
      params = annotation$params
    ),
    metadata_payload
  )
  se
}

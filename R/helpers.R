`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

require_namespace <- function(pkg, reason = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message <- sprintf("Package '%s' is required", pkg)
    if (!is.null(reason) && nzchar(reason)) {
      message <- sprintf("%s for %s", message, reason)
    }
    stop(message, call. = FALSE)
  }
}

assert_data_frame <- function(x, name) {
  if (!is.data.frame(x)) {
    stop(sprintf("`%s` must be a data.frame.", name), call. = FALSE)
  }
}

assert_first_col <- function(x, expected, name) {
  if (ncol(x) < 1L || names(x)[1] != expected) {
    stop(
      sprintf("`%s` must have `%s` as its first column.", name, expected),
      call. = FALSE
    )
  }
}

normalize_gene_ids <- function(gene_ids, strip_version = TRUE) {
  raw_ids <- as.character(gene_ids)
  normalized <- raw_ids
  if (isTRUE(strip_version)) {
    normalized <- sub("\\..*$", "", normalized)
  }

  data.frame(
    gene_id_raw = raw_ids,
    gene_id = normalized,
    stringsAsFactors = FALSE
  )
}

detect_species_from_ids <- function(gene_ids) {
  if (length(gene_ids) == 0L) {
    return("unknown")
  }

  human <- grepl("^ENSG", gene_ids)
  mouse <- grepl("^ENSMUSG", gene_ids)

  if (all(human)) {
    return("human")
  }
  if (all(mouse)) {
    return("mouse")
  }
  if (any(human) && any(mouse)) {
    stop("Mixed human and mouse Ensembl IDs were detected.", call. = FALSE)
  }
  "unknown"
}

resolve_species <- function(species, gene_ids) {
  species <- match.arg(species, c("auto", "human", "mouse"))
  if (species != "auto") {
    return(species)
  }

  detected <- detect_species_from_ids(gene_ids)
  if (detected == "unknown") {
    stop("Unable to infer species from `gene_id`. Please set `species` explicitly.", call. = FALSE)
  }
  detected
}

match_meta_to_expr <- function(expr, meta) {
  sample_cols <- names(expr)[-1]
  meta_samples <- as.character(meta[[1]])

  missing_samples <- setdiff(sample_cols, meta_samples)
  if (length(missing_samples) > 0L) {
    stop(
      sprintf(
        "Sample columns in `expr` missing from `meta$sample`: %s",
        paste(missing_samples, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  matched <- meta[match(sample_cols, meta_samples), , drop = FALSE]
  rownames(matched) <- NULL
  matched
}

write_if_requested <- function(x, output_file = NULL, row.names = FALSE) {
  if (!is.null(output_file)) {
    utils::write.csv(x, file = output_file, row.names = row.names)
  }
  invisible(x)
}

merge_fill <- function(base, incoming, field) {
  if (!field %in% names(incoming)) {
    return(base)
  }
  idx <- match(base$gene_id, incoming$gene_id)
  replacement <- incoming[[field]][idx]
  if (!field %in% names(base)) {
    base[[field]] <- replacement
    return(base)
  }
  take <- is.na(base[[field]]) | base[[field]] == ""
  base[[field]][take] <- replacement[take]
  base
}

collapse_symbol_matrix <- function(expr_anno, gene_column = "symbol") {
  if (!gene_column %in% names(expr_anno)) {
    stop(sprintf("`expr_anno` does not contain `%s`.", gene_column), call. = FALSE)
  }

  known_annotation <- c(
    "gene_id_raw", "gene_id", "symbol", "gene_name", "entrez_id", "biotype",
    "chromosome", "start", "end", "strand", "annotation_source",
    "annotation_status", "symbol_source", "gene_name_source", "entrez_id_source",
    "biotype_source", "chromosome_source", "start_source", "end_source",
    "strand_source"
  )
  sample_cols <- setdiff(names(expr_anno), known_annotation)
  keep <- !is.na(expr_anno[[gene_column]]) & expr_anno[[gene_column]] != ""
  expr_mat <- as.matrix(expr_anno[keep, sample_cols, drop = FALSE])
  mode(expr_mat) <- "numeric"
  aggregated <- stats::aggregate(expr_mat, by = list(gene = expr_anno[[gene_column]][keep]), FUN = sum)
  rownames(aggregated) <- aggregated$gene
  as.matrix(aggregated[, -1, drop = FALSE])
}

default_annotation_fields <- function() {
  c("symbol", "gene_name", "entrez_id", "biotype", "chromosome", "start", "end", "strand")
}

annotation_coverage_report <- function(annotation, fields) {
  rates <- vapply(fields, function(field) {
    if (!field %in% names(annotation)) {
      return(0)
    }
    mean(!is.na(annotation[[field]]) & annotation[[field]] != "")
  }, numeric(1))

  data.frame(
    field = fields,
    annotation_rate = rates,
    stringsAsFactors = FALSE
  )
}

new_expranno_annotation <- function(expr_anno, annotation, meta_checked, report, params) {
  out <- list(
    expr_anno = expr_anno,
    annotation = annotation,
    meta_checked = meta_checked,
    report = report,
    params = params
  )
  class(out) <- c("expranno_annotation", "list")
  out
}

new_expranno_result <- function(annotation, expr_meta_merged, deconvolution, signatures, files) {
  out <- list(
    annotation = annotation,
    expr_meta_merged = expr_meta_merged,
    deconvolution = deconvolution,
    signatures = signatures,
    files = files
  )
  class(out) <- c("expranno_result", "list")
  out
}

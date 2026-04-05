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

clean_annotation_value <- function(x) {
  out <- trimws(as.character(x))
  out[out %in% c("", "NA", "NANA", "NULL")] <- NA_character_
  out
}

non_missing_unique <- function(x) {
  cleaned <- clean_annotation_value(x)
  unique(cleaned[!is.na(cleaned)])
}

candidate_values_from_string <- function(x) {
  if (length(x) == 0L || is.na(x) || !nzchar(x)) {
    return(character(0))
  }
  non_missing_unique(strsplit(x, ";", fixed = TRUE)[[1]])
}

collapse_candidate_strings <- function(values) {
  flattened <- unlist(lapply(values, candidate_values_from_string), use.names = FALSE)
  unique_vals <- non_missing_unique(flattened)
  list(
    value = if (length(unique_vals) > 0L) unique_vals[[1]] else NA_character_,
    candidates = if (length(unique_vals) > 0L) paste(unique_vals, collapse = ";") else NA_character_,
    candidate_count = length(unique_vals),
    is_ambiguous = length(unique_vals) > 1L
  )
}

collapse_annotation_field <- function(values) {
  unique_vals <- non_missing_unique(values)
  list(
    value = if (length(unique_vals) > 0L) unique_vals[[1]] else NA_character_,
    candidates = if (length(unique_vals) > 0L) paste(unique_vals, collapse = ";") else NA_character_,
    candidate_count = length(unique_vals),
    is_ambiguous = length(unique_vals) > 1L
  )
}

write_text_if_requested <- function(text, output_file = NULL) {
  if (!is.null(output_file)) {
    writeLines(text, con = output_file, useBytes = TRUE)
  }
  invisible(text)
}

session_info_text <- function() {
  paste(capture.output(utils::sessionInfo()), collapse = "\n")
}

known_annotation_columns <- function() {
  c(
    "gene_id_raw", "gene_id", "symbol", "gene_name", "entrez_id", "biotype",
    "chromosome", "start", "end", "strand", "annotation_source",
    "annotation_status", "annotation_backend_release", "annotation_backend_host",
    "annotation_backend_mirror", "annotation_date", "symbol_source", "gene_name_source", "entrez_id_source",
    "biotype_source", "chromosome_source", "start_source", "end_source",
    "strand_source", "symbol_candidates", "gene_name_candidates",
    "entrez_id_candidates", "biotype_candidates", "chromosome_candidates",
    "start_candidates", "end_candidates", "strand_candidates",
    "symbol_candidate_count", "gene_name_candidate_count",
    "entrez_id_candidate_count", "biotype_candidate_count",
    "chromosome_candidate_count", "start_candidate_count",
    "end_candidate_count", "strand_candidate_count",
    "symbol_is_ambiguous", "gene_name_is_ambiguous", "entrez_id_is_ambiguous",
    "biotype_is_ambiguous", "chromosome_is_ambiguous", "start_is_ambiguous",
    "end_is_ambiguous", "strand_is_ambiguous"
  )
}

sample_columns_from_expr_anno <- function(expr_anno) {
  setdiff(names(expr_anno), known_annotation_columns())
}

numeric_expression_matrix <- function(expr_anno, sample_cols, keep) {
  expr_mat <- as.matrix(expr_anno[keep, sample_cols, drop = FALSE])
  mode(expr_mat) <- "numeric"
  expr_mat
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

normalize_expr_scale <- function(expr_scale, expr_matrix = NULL) {
  expr_scale <- match.arg(expr_scale, c("auto", "count", "abundance", "log"))
  if (expr_scale != "auto" || is.null(expr_matrix)) {
    return(expr_scale)
  }

  values <- as.numeric(expr_matrix)
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return("abundance")
  }

  if (all(values >= 0) && all(abs(values - round(values)) < 1e-8)) {
    return("count")
  }
  if (any(values < 0)) {
    return("log")
  }
  "abundance"
}

resolve_duplicate_strategy <- function(duplicate_strategy, expr_scale, expr_matrix = NULL) {
  duplicate_strategy <- match.arg(
    duplicate_strategy,
    c("auto", "sum", "mean", "max", "first")
  )
  if (duplicate_strategy != "auto") {
    return(duplicate_strategy)
  }

  normalized_scale <- normalize_expr_scale(expr_scale = expr_scale, expr_matrix = expr_matrix)
  if (identical(normalized_scale, "count")) {
    return("sum")
  }
  "mean"
}

aggregate_duplicate_rows <- function(expr_matrix, gene_ids, duplicate_strategy) {
  if (length(gene_ids) == 0L || nrow(expr_matrix) == 0L) {
    out <- expr_matrix[0, , drop = FALSE]
    rownames(out) <- character(0)
    return(out)
  }

  if (identical(duplicate_strategy, "first")) {
    keep <- !duplicated(gene_ids)
    out <- expr_matrix[keep, , drop = FALSE]
    rownames(out) <- gene_ids[keep]
    return(out)
  }

  aggregate_fun <- switch(
    duplicate_strategy,
    sum = base::sum,
    mean = base::mean,
    max = base::max,
    stop(sprintf("Unsupported duplicate strategy: %s", duplicate_strategy), call. = FALSE)
  )

  aggregated <- stats::aggregate(expr_matrix, by = list(gene = gene_ids), FUN = aggregate_fun)
  rownames(aggregated) <- aggregated$gene
  as.matrix(aggregated[, -1, drop = FALSE])
}

collapse_symbol_matrix <- function(
    expr_anno,
    gene_column = "symbol",
    expr_scale = c("auto", "count", "abundance", "log"),
    duplicate_strategy = c("auto", "sum", "mean", "max", "first")) {
  if (!gene_column %in% names(expr_anno)) {
    stop(sprintf("`expr_anno` does not contain `%s`.", gene_column), call. = FALSE)
  }

  sample_cols <- sample_columns_from_expr_anno(expr_anno)
  if (length(sample_cols) == 0L) {
    stop("`expr_anno` does not contain any sample columns.", call. = FALSE)
  }
  keep <- !is.na(expr_anno[[gene_column]]) & expr_anno[[gene_column]] != ""
  expr_mat <- numeric_expression_matrix(expr_anno, sample_cols = sample_cols, keep = keep)
  strategy <- resolve_duplicate_strategy(
    duplicate_strategy = duplicate_strategy,
    expr_scale = expr_scale,
    expr_matrix = expr_mat
  )
  aggregate_duplicate_rows(
    expr_matrix = expr_mat,
    gene_ids = expr_anno[[gene_column]][keep],
    duplicate_strategy = strategy
  )
}

default_annotation_fields <- function() {
  c("symbol", "gene_name", "entrez_id", "biotype", "chromosome", "start", "end", "strand")
}

available_annotation_presets <- function() {
  c(
    "human_v102",
    "mouse_v102",
    "human_tpm_v102",
    "mouse_tpm_v102",
    "human_count_v102",
    "mouse_count_v102"
  )
}

annotation_preset_defaults <- function(annotation_preset) {
  if (is.null(annotation_preset)) {
    return(NULL)
  }

  annotation_preset <- match.arg(annotation_preset, available_annotation_presets())
  species <- if (grepl("^human", annotation_preset)) "human" else "mouse"
  expr_scale <- if (grepl("_tpm_", annotation_preset)) "abundance" else if (grepl("_count_", annotation_preset)) "count" else "auto"

  list(
    annotation_preset = annotation_preset,
    species = species,
    annotation_engine = "hybrid",
    strip_version = TRUE,
    biomart_version = 102,
    biomart_host = NULL,
    biomart_mirror = NULL,
    expr_scale = expr_scale,
    duplicate_strategy = "auto"
  )
}

truth_fields_from_truth <- function(truth, truth_gene_col) {
  setdiff(names(truth), truth_gene_col)
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

annotation_ambiguity_report <- function(annotation, fields) {
  rows <- list()
  idx <- 1L
  for (field in fields) {
    count_col <- paste0(field, "_candidate_count")
    value_col <- field
    candidate_col <- paste0(field, "_candidates")
    source_col <- paste0(field, "_source")
    if (!all(c(count_col, value_col, candidate_col, source_col) %in% names(annotation))) {
      next
    }
    keep <- !is.na(annotation[[count_col]]) & annotation[[count_col]] > 1L
    if (!any(keep)) {
      next
    }
    rows[[idx]] <- data.frame(
      gene_id = annotation$gene_id[keep],
      field = field,
      chosen_value = annotation[[value_col]][keep],
      candidate_count = annotation[[count_col]][keep],
      candidates = annotation[[candidate_col]][keep],
      source = annotation[[source_col]][keep],
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
  if (length(rows) == 0L) {
    return(data.frame(
      gene_id = character(0),
      field = character(0),
      chosen_value = character(0),
      candidate_count = integer(0),
      candidates = character(0),
      source = character(0),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

annotation_summary_metrics <- function(annotation, report, ambiguity_report) {
  data.frame(
    annotated_genes = sum(annotation$annotation_status == "annotated", na.rm = TRUE),
    total_genes = nrow(annotation),
    ambiguous_gene_fields = nrow(ambiguity_report),
    symbol_coverage = report$annotation_rate[match("symbol", report$field)] %||% NA_real_,
    gene_name_coverage = report$annotation_rate[match("gene_name", report$field)] %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

new_expranno_benchmark <- function(summary, coverage, runs, params) {
  out <- list(
    summary = summary,
    coverage = coverage,
    runs = runs,
    params = params
  )
  class(out) <- c("expranno_benchmark", "list")
  out
}

new_expranno_validation <- function(summary, details, runs, params) {
  out <- list(
    summary = summary,
    details = details,
    runs = runs,
    params = params
  )
  class(out) <- c("expranno_validation", "list")
  out
}

new_expranno_annotation <- function(expr_anno, annotation, meta_checked, report, ambiguity_report, provenance, params) {
  out <- list(
    expr_anno = expr_anno,
    annotation = annotation,
    meta_checked = meta_checked,
    report = report,
    ambiguity_report = ambiguity_report,
    provenance = provenance,
    params = params
  )
  class(out) <- c("expranno_annotation", "list")
  out
}

new_expranno_result <- function(annotation, expr_meta_merged, deconvolution, signatures, files, benchmark = NULL, validation = NULL, session_info = NULL) {
  out <- list(
    annotation = annotation,
    expr_meta_merged = expr_meta_merged,
    deconvolution = deconvolution,
    signatures = signatures,
    files = files,
    benchmark = benchmark,
    validation = validation,
    session_info = session_info
  )
  class(out) <- c("expranno_result", "list")
  out
}

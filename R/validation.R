normalize_truth_table <- function(truth, truth_gene_col = "gene_id", strip_version = TRUE) {
  assert_data_frame(truth, "truth")
  if (!truth_gene_col %in% names(truth)) {
    stop(sprintf("`truth` must contain `%s`.", truth_gene_col), call. = FALSE)
  }

  truth_gene_info <- normalize_gene_ids(truth[[truth_gene_col]], strip_version = strip_version)
  out <- data.frame(
    gene_id = truth_gene_info$gene_id,
    truth[, setdiff(names(truth), truth_gene_col), drop = FALSE],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out <- out[!duplicated(out$gene_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

annotation_validation_detail <- function(annotation, truth, fields, engine) {
  merged <- merge(
    annotation,
    truth,
    by = "gene_id",
    all.y = TRUE,
    sort = FALSE,
    suffixes = c("_predicted", "_truth")
  )

  rows <- list()
  idx <- 1L
  for (field in fields) {
    predicted_col <- if (paste0(field, "_predicted") %in% names(merged)) paste0(field, "_predicted") else field
    truth_col <- if (paste0(field, "_truth") %in% names(merged)) paste0(field, "_truth") else field
    if (!truth_col %in% names(merged) || !predicted_col %in% names(merged)) {
      next
    }
    predicted <- clean_annotation_value(merged[[predicted_col]])
    expected <- clean_annotation_value(merged[[truth_col]])
    has_truth <- !is.na(expected)
    is_match <- has_truth & !is.na(predicted) & predicted == expected
    is_missing_prediction <- has_truth & is.na(predicted)
    is_mismatch <- has_truth & !is.na(predicted) & predicted != expected
    candidate_col <- paste0(field, "_candidates")
    source_col <- paste0(field, "_source")
    ambiguous_col <- paste0(field, "_is_ambiguous")

    rows[[idx]] <- data.frame(
      engine = engine,
      gene_id = merged$gene_id,
      field = field,
      truth_value = expected,
      predicted_value = predicted,
      is_match = is_match,
      is_missing_prediction = is_missing_prediction,
      is_mismatch = is_mismatch,
      is_ambiguous = if (ambiguous_col %in% names(merged)) merged[[ambiguous_col]] else FALSE,
      candidates = if (candidate_col %in% names(merged)) merged[[candidate_col]] else NA_character_,
      source = if (source_col %in% names(merged)) merged[[source_col]] else NA_character_,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }

  if (length(rows) == 0L) {
    return(data.frame(
      engine = character(0),
      gene_id = character(0),
      field = character(0),
      truth_value = character(0),
      predicted_value = character(0),
      is_match = logical(0),
      is_missing_prediction = logical(0),
      is_mismatch = logical(0),
      is_ambiguous = logical(0),
      candidates = character(0),
      source = character(0),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

annotation_validation_summary <- function(details) {
  split_details <- split(details, paste(details$engine, details$field, sep = "::"))
  rows <- lapply(split_details, function(one) {
    truth_rows <- sum(!is.na(one$truth_value) & one$truth_value != "")
    matched_rows <- sum(one$is_match, na.rm = TRUE)
    missing_prediction_rows <- sum(one$is_missing_prediction, na.rm = TRUE)
    mismatch_rows <- sum(one$is_mismatch, na.rm = TRUE)
    ambiguous_rows <- sum(one$is_ambiguous & !is.na(one$truth_value) & one$truth_value != "", na.rm = TRUE)

    data.frame(
      engine = one$engine[[1]],
      field = one$field[[1]],
      truth_rows = truth_rows,
      matched_rows = matched_rows,
      missing_prediction_rows = missing_prediction_rows,
      mismatch_rows = mismatch_rows,
      ambiguous_rows = ambiguous_rows,
      match_rate = if (truth_rows > 0L) matched_rows / truth_rows else NA_real_,
      missing_prediction_rate = if (truth_rows > 0L) missing_prediction_rows / truth_rows else NA_real_,
      mismatch_rate = if (truth_rows > 0L) mismatch_rows / truth_rows else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Validate annotation engines against a truth table
#'
#' Runs one or more annotation engines, then compares the chosen annotation
#' fields against a user-supplied truth table keyed by Ensembl gene ID.
#'
#' @param expr Expression table with `gene_id` first, or a
#'   `SummarizedExperiment`-like object.
#' @param meta Metadata table with `sample` first. Leave `NULL` when `expr` is a
#'   `SummarizedExperiment`.
#' @param truth A data frame containing a `gene_id` column and one or more truth
#'   annotation fields such as `symbol` or `gene_name`.
#' @param truth_gene_col Column in `truth` containing the Ensembl gene ID.
#' @param species Either `"auto"`, `"human"`, or `"mouse"`.
#' @param annotation_preset Optional preset that fixes a reproducible
#'   annotation configuration. Supported values are `"human_v102"`,
#'   `"mouse_v102"`, `"human_tpm_v102"`, `"mouse_tpm_v102"`,
#'   `"human_count_v102"`, and `"mouse_count_v102"`.
#' @param engines Character vector of annotation engines to compare.
#' @param fields Optional annotation fields to validate. Defaults to the
#'   intersection between `truth` and the package's standard annotation fields.
#' @param strip_version Whether to remove Ensembl version suffixes in both
#'   prediction and truth tables.
#' @param biomart_version Fixed Ensembl release used by the `biomaRt` backend.
#' @param biomart_host Optional explicit Ensembl host.
#' @param biomart_mirror Optional Ensembl mirror name.
#' @param assay_name Optional assay name when `expr` is a
#'   `SummarizedExperiment`.
#' @param gene_id_col Optional row-data column containing Ensembl IDs when
#'   `expr` is a `SummarizedExperiment`.
#' @param sample_col Metadata column to use as `sample` when `expr` is a
#'   `SummarizedExperiment`.
#' @param output_file Optional CSV path for the validation summary.
#' @param detail_file Optional CSV path for per-gene validation details.
#' @param verbose Whether to emit progress messages.
#'
#' @return An `expranno_validation` object.
#' @export
validate_annotation_engines <- function(
    expr,
    meta = NULL,
    truth,
    truth_gene_col = "gene_id",
    species = c("auto", "human", "mouse"),
    annotation_preset = NULL,
    engines = c("biomart", "orgdb", "ensdb", "hybrid"),
    fields = NULL,
    strip_version = TRUE,
    biomart_version = 102,
    biomart_host = NULL,
    biomart_mirror = NULL,
    assay_name = NULL,
    gene_id_col = NULL,
    sample_col = "sample",
    output_file = NULL,
    detail_file = NULL,
    verbose = TRUE) {
  preset <- annotation_preset_defaults(annotation_preset)
  species <- if (is.null(preset)) match.arg(species) else preset$species
  if (!is.null(preset)) {
    strip_version <- preset$strip_version
    biomart_version <- preset$biomart_version
    biomart_host <- preset$biomart_host
    biomart_mirror <- preset$biomart_mirror
  }

  truth_norm <- normalize_truth_table(
    truth = truth,
    truth_gene_col = truth_gene_col,
    strip_version = strip_version
  )
  fields <- fields %||% intersect(default_annotation_fields(), truth_fields_from_truth(truth_norm, "gene_id"))
  if (length(fields) == 0L) {
    stop("No overlapping truth fields were found to validate.", call. = FALSE)
  }

  runs <- list()
  detail_rows <- list()
  idx <- 1L

  for (engine in engines) {
    if (isTRUE(verbose)) {
      message(sprintf("Validating annotation engine: %s", engine))
    }

    run <- annotate_expr(
      expr = expr,
      meta = meta,
      species = species,
      annotation_preset = annotation_preset,
      annotation_engine = engine,
      fields = unique(c(fields, default_annotation_fields())),
      strip_version = strip_version,
      biomart_version = biomart_version,
      biomart_host = biomart_host,
      biomart_mirror = biomart_mirror,
      assay_name = assay_name,
      gene_id_col = gene_id_col,
      sample_col = sample_col,
      verbose = verbose
    )
    runs[[engine]] <- run
    detail_rows[[idx]] <- annotation_validation_detail(
      annotation = run$annotation,
      truth = truth_norm,
      fields = fields,
      engine = engine
    )
    idx <- idx + 1L
  }

  details <- do.call(rbind, detail_rows)
  rownames(details) <- NULL
  summary <- annotation_validation_summary(details)

  write_if_requested(summary, output_file = output_file)
  write_if_requested(details, output_file = detail_file)

  new_expranno_validation(
    summary = summary,
    details = details,
    runs = runs,
    params = list(
      species = species,
      annotation_preset = annotation_preset,
      engines = engines,
      fields = fields,
      strip_version = strip_version,
      biomart_version = biomart_version,
      biomart_host = biomart_host,
      biomart_mirror = biomart_mirror
    )
  )
}

#' Benchmark annotation engines on the same input
#'
#' Runs one or more annotation engines on the same expression and metadata
#' inputs, then summarizes coverage, ambiguity, and provenance for comparison.
#'
#' @param expr Expression table with `gene_id` first, or a
#'   `SummarizedExperiment`-like object.
#' @param meta Metadata table with `sample` first. Leave `NULL` when `expr` is a
#'   `SummarizedExperiment`.
#' @param species Either `"auto"`, `"human"`, or `"mouse"`.
#' @param annotation_preset Optional preset that fixes a reproducible
#'   annotation configuration. Supported values are `"human_v102"`,
#'   `"mouse_v102"`, `"human_tpm_v102"`, `"mouse_tpm_v102"`,
#'   `"human_count_v102"`, and `"mouse_count_v102"`.
#' @param engines Character vector of annotation engines to compare.
#' @param fields Optional annotation fields to benchmark.
#' @param strip_version Whether to remove Ensembl version suffixes.
#' @param biomart_version Fixed Ensembl release used by the `biomaRt` backend.
#' @param biomart_host Optional explicit Ensembl host.
#' @param biomart_mirror Optional Ensembl mirror name.
#' @param assay_name Optional assay name when `expr` is a
#'   `SummarizedExperiment`.
#' @param gene_id_col Optional row-data column containing Ensembl IDs when
#'   `expr` is a `SummarizedExperiment`.
#' @param sample_col Metadata column to use as `sample` when `expr` is a
#'   `SummarizedExperiment`.
#' @param output_file Optional CSV path for the benchmark summary.
#' @param coverage_file Optional CSV path for long-format coverage results.
#' @param verbose Whether to emit progress messages.
#'
#' @return An `expranno_benchmark` object.
#' @export
benchmark_annotation_engines <- function(
    expr,
    meta = NULL,
    species = c("auto", "human", "mouse"),
    annotation_preset = NULL,
    engines = c("none", "biomart", "orgdb", "ensdb", "hybrid"),
    fields = NULL,
    strip_version = TRUE,
    biomart_version = 102,
    biomart_host = NULL,
    biomart_mirror = NULL,
    assay_name = NULL,
    gene_id_col = NULL,
    sample_col = "sample",
    output_file = NULL,
    coverage_file = NULL,
    verbose = TRUE) {
  preset <- annotation_preset_defaults(annotation_preset)
  fields <- fields %||% default_annotation_fields()
  species <- if (is.null(preset)) match.arg(species) else preset$species
  if (!is.null(preset)) {
    strip_version <- preset$strip_version
    biomart_version <- preset$biomart_version
    biomart_host <- preset$biomart_host
    biomart_mirror <- preset$biomart_mirror
  }
  runs <- list()
  summary_rows <- list()
  coverage_rows <- list()

  for (engine in engines) {
    if (isTRUE(verbose)) {
      message(sprintf("Benchmarking annotation engine: %s", engine))
    }

    run <- tryCatch(
      annotate_expr(
        expr = expr,
        meta = meta,
        species = species,
        annotation_preset = annotation_preset,
        annotation_engine = engine,
        fields = fields,
        strip_version = strip_version,
        biomart_version = biomart_version,
        biomart_host = biomart_host,
        biomart_mirror = biomart_mirror,
        assay_name = assay_name,
        gene_id_col = gene_id_col,
        sample_col = sample_col,
        verbose = verbose
      ),
      error = function(e) e
    )

    runs[[engine]] <- run
    if (inherits(run, "error")) {
      summary_rows[[engine]] <- data.frame(
        engine = engine,
        status = "error",
        species = if (identical(species, "auto")) NA_character_ else species,
        annotated_genes = NA_integer_,
        total_genes = NA_integer_,
        ambiguous_gene_fields = NA_integer_,
        symbol_coverage = NA_real_,
        gene_name_coverage = NA_real_,
        message = conditionMessage(run),
        stringsAsFactors = FALSE
      )
      next
    }

    metrics <- annotation_summary_metrics(
      annotation = run$annotation,
      report = run$report,
      ambiguity_report = run$ambiguity_report
    )
    summary_rows[[engine]] <- cbind(
      data.frame(
        engine = engine,
        status = "ok",
        species = run$params$species,
        message = NA_character_,
        stringsAsFactors = FALSE
      ),
      metrics,
      stringsAsFactors = FALSE
    )

    coverage_rows[[engine]] <- data.frame(
      engine = engine,
      run$report,
      stringsAsFactors = FALSE
    )
  }

  summary <- do.call(rbind, summary_rows)
  coverage <- if (length(coverage_rows) == 0L) {
    data.frame(engine = character(0), field = character(0), annotation_rate = numeric(0), stringsAsFactors = FALSE)
  } else {
    do.call(rbind, coverage_rows)
  }

  write_if_requested(summary, output_file = output_file)
  write_if_requested(coverage, output_file = coverage_file)

  new_expranno_benchmark(
    summary = summary,
    coverage = coverage,
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

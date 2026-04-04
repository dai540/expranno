#' Run the full expranno workflow
#'
#' This wrapper annotates genes, writes `expr_anno.csv`, writes provenance and
#' ambiguity reports, merges expression with metadata, optionally runs
#' deconvolution and signature scoring, optionally benchmarks annotation
#' engines, and returns all outputs in one structured object.
#'
#' @param expr Expression table with `gene_id` in the first column, or a
#'   `SummarizedExperiment`-like object.
#' @param meta Metadata table with `sample` in the first column. Leave `NULL`
#'   when `expr` is a `SummarizedExperiment`.
#' @param species Either `"auto"`, `"human"`, or `"mouse"`.
#' @param annotation_preset Optional preset that fixes a reproducible
#'   annotation configuration. Supported values are `"human_v102"`,
#'   `"mouse_v102"`, `"human_tpm_v102"`, `"mouse_tpm_v102"`,
#'   `"human_count_v102"`, and `"mouse_count_v102"`.
#' @param annotation_engine Annotation backend strategy.
#' @param output_dir Output directory.
#' @param biomart_version Fixed Ensembl release used by the `biomaRt` backend.
#' @param biomart_host Optional explicit Ensembl host.
#' @param biomart_mirror Optional Ensembl mirror name.
#' @param assay_name Optional assay name when `expr` is a
#'   `SummarizedExperiment`.
#' @param gene_id_col Optional row-data column containing Ensembl IDs when
#'   `expr` is a `SummarizedExperiment`.
#' @param sample_col Metadata column to use as `sample` when `expr` is a
#'   `SummarizedExperiment`.
#' @param run_deconvolution Whether to run `immunedeconv`.
#' @param deconv_methods Either `"all_except_cibersort"` or an explicit vector.
#' @param expr_scale Expression scale used to choose duplicate-symbol handling.
#' @param duplicate_strategy Strategy used when multiple rows map to the same
#'   symbol. `"auto"` uses `"sum"` for counts and `"mean"` otherwise.
#' @param deconv_args Optional named list of additional arguments passed to
#'   [run_cell_deconvolution()], including method-specific values such as
#'   `indications`.
#' @param run_signature Whether to run GSVA or ssGSEA.
#' @param geneset_file Optional GMT path.
#' @param gene_sets Optional named list of gene sets.
#' @param signature_method One of `"gsva"`, `"ssgsea"`, or `"both"`.
#' @param signature_kcdf Kernel choice passed to [run_signature_analysis()].
#' @param signature_min_size Minimum gene-set size after ID matching.
#' @param signature_max_size Maximum gene-set size after ID matching.
#' @param gsva_args Optional named list of extra arguments passed to
#'   `GSVA::gsvaParam()`.
#' @param ssgsea_args Optional named list of extra arguments passed to
#'   `GSVA::ssgseaParam()`.
#' @param run_benchmark Whether to run [benchmark_annotation_engines()] on the
#'   same inputs.
#' @param benchmark_engines Annotation engines to benchmark when
#'   `run_benchmark = TRUE`.
#' @param run_validation Whether to run [validate_annotation_engines()] against
#'   a supplied truth table.
#' @param validation_truth Optional truth table keyed by Ensembl gene ID.
#' @param validation_truth_gene_col Column in `validation_truth` containing the
#'   Ensembl gene ID.
#' @param validation_fields Optional truth fields to validate, such as
#'   `symbol`.
#' @param save_session_info Whether to write `session_info.txt`.
#' @param verbose Whether to emit progress messages.
#'
#' @return An `expranno_result` object.
#' @export
run_expranno <- function(
    expr,
    meta = NULL,
    species = c("auto", "human", "mouse"),
    annotation_preset = NULL,
    annotation_engine = c("hybrid", "biomart", "orgdb", "ensdb", "none"),
    output_dir = ".",
    biomart_version = 102,
    biomart_host = NULL,
    biomart_mirror = NULL,
    assay_name = NULL,
    gene_id_col = NULL,
    sample_col = "sample",
    run_deconvolution = TRUE,
    deconv_methods = "all_except_cibersort",
    expr_scale = c("auto", "count", "abundance", "log"),
    duplicate_strategy = c("auto", "sum", "mean", "max", "first"),
    deconv_args = list(),
    run_signature = TRUE,
    geneset_file = NULL,
    gene_sets = NULL,
    signature_method = c("gsva", "ssgsea", "both"),
    signature_kcdf = c("auto", "Gaussian", "Poisson", "none"),
    signature_min_size = 1,
    signature_max_size = Inf,
    gsva_args = list(),
    ssgsea_args = list(),
    run_benchmark = FALSE,
    benchmark_engines = c("none", "biomart", "orgdb", "ensdb", "hybrid"),
    run_validation = FALSE,
    validation_truth = NULL,
    validation_truth_gene_col = "gene_id",
    validation_fields = NULL,
    save_session_info = TRUE,
    verbose = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  preset <- annotation_preset_defaults(annotation_preset)
  species <- if (is.null(preset)) match.arg(species) else preset$species
  annotation_engine <- if (is.null(preset)) match.arg(annotation_engine) else preset$annotation_engine
  expr_scale <- match.arg(expr_scale)
  duplicate_strategy <- match.arg(duplicate_strategy)
  signature_method <- match.arg(signature_method)
  signature_kcdf <- match.arg(signature_kcdf)
  if (!is.null(preset)) {
    biomart_version <- preset$biomart_version
    biomart_host <- preset$biomart_host
    biomart_mirror <- preset$biomart_mirror
    if (identical(expr_scale, "auto")) {
      expr_scale <- preset$expr_scale
    }
    if (identical(duplicate_strategy, "auto")) {
      duplicate_strategy <- preset$duplicate_strategy
    }
  }

  annotation <- annotate_expr(
    expr = expr,
    meta = meta,
    species = species,
    annotation_preset = annotation_preset,
    annotation_engine = annotation_engine,
    biomart_version = biomart_version,
    biomart_host = biomart_host,
    biomart_mirror = biomart_mirror,
    assay_name = assay_name,
    gene_id_col = gene_id_col,
    sample_col = sample_col,
    output_file = file.path(output_dir, "expr_anno.csv"),
    report_file = file.path(output_dir, "annotation_report.csv"),
    ambiguity_file = file.path(output_dir, "annotation_ambiguity.csv"),
    provenance_file = file.path(output_dir, "annotation_provenance.csv"),
    verbose = verbose
  )

  expr_meta_merged <- merge_expr_meta(
    expr_anno = annotation$expr_anno,
    meta = annotation$meta_checked,
    output_file = file.path(output_dir, "expr_meta_merged.csv")
  )

  deconvolution <- NULL
  if (isTRUE(run_deconvolution)) {
    deconvolution <- do.call(
      run_cell_deconvolution,
      c(
        list(
          expr_anno = annotation$expr_anno,
          species = annotation$params$species,
          methods = deconv_methods,
          expr_scale = expr_scale,
          duplicate_strategy = duplicate_strategy,
          output_dir = output_dir,
          verbose = verbose
        ),
        deconv_args
      )
    )
  }

  signatures <- NULL
  if (isTRUE(run_signature)) {
    signatures <- run_signature_analysis(
      expr_anno = annotation$expr_anno,
      geneset_file = geneset_file,
      gene_sets = gene_sets,
      method = signature_method,
      expr_scale = expr_scale,
      duplicate_strategy = duplicate_strategy,
      kcdf = signature_kcdf,
      min_size = signature_min_size,
      max_size = signature_max_size,
      gsva_args = gsva_args,
      ssgsea_args = ssgsea_args,
      output_dir = output_dir
    )
  }

  benchmark <- NULL
  if (isTRUE(run_benchmark)) {
    benchmark <- benchmark_annotation_engines(
      expr = expr,
      meta = meta,
      species = species,
      annotation_preset = annotation_preset,
      engines = benchmark_engines,
      strip_version = annotation$params$strip_version,
      biomart_version = biomart_version,
      biomart_host = biomart_host,
      biomart_mirror = biomart_mirror,
      assay_name = assay_name,
      gene_id_col = gene_id_col,
      sample_col = sample_col,
      output_file = file.path(output_dir, "annotation_benchmark_summary.csv"),
      coverage_file = file.path(output_dir, "annotation_benchmark_coverage.csv"),
      verbose = verbose
    )
  }

  validation <- NULL
  if (isTRUE(run_validation)) {
    if (is.null(validation_truth)) {
      stop("Provide `validation_truth` when `run_validation = TRUE`.", call. = FALSE)
    }
    validation_engines <- if (isTRUE(run_benchmark)) {
      unique(c(annotation_engine, benchmark_engines[benchmark_engines != "none"]))
    } else {
      annotation_engine
    }
    validation <- validate_annotation_engines(
      expr = expr,
      meta = meta,
      truth = validation_truth,
      truth_gene_col = validation_truth_gene_col,
      species = species,
      annotation_preset = annotation_preset,
      engines = validation_engines,
      fields = validation_fields,
      strip_version = annotation$params$strip_version,
      biomart_version = biomart_version,
      biomart_host = biomart_host,
      biomart_mirror = biomart_mirror,
      assay_name = assay_name,
      gene_id_col = gene_id_col,
      sample_col = sample_col,
      output_file = file.path(output_dir, "annotation_validation_summary.csv"),
      detail_file = file.path(output_dir, "annotation_validation_detail.csv"),
      verbose = verbose
    )
  }

  session_info <- NULL
  if (isTRUE(save_session_info)) {
    session_info <- session_info_text()
    write_text_if_requested(
      session_info,
      output_file = file.path(output_dir, "session_info.txt")
    )
  }

  files <- list.files(output_dir, full.names = TRUE)
  new_expranno_result(
    annotation = annotation,
    expr_meta_merged = expr_meta_merged,
    deconvolution = deconvolution,
    signatures = signatures,
    files = files,
    benchmark = benchmark,
    validation = validation,
    session_info = session_info
  )
}

#' Run the full expranno workflow
#'
#' This wrapper annotates genes, writes `expr_anno.csv`, merges expression with
#' metadata, optionally runs immune deconvolution and signature scoring, and
#' returns all outputs in one structured object.
#'
#' @param expr Expression table with `gene_id` in the first column.
#' @param meta Metadata table with `sample` in the first column.
#' @param species Either `"auto"`, `"human"`, or `"mouse"`.
#' @param annotation_engine Annotation backend strategy.
#' @param output_dir Output directory.
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
#' @param verbose Whether to emit progress messages.
#'
#' @return An `expranno_result` object.
#' @export
run_expranno <- function(
    expr,
    meta,
    species = c("auto", "human", "mouse"),
    annotation_engine = c("hybrid", "biomart", "orgdb", "ensdb", "none"),
    output_dir = ".",
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
    verbose = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  expr_scale <- match.arg(expr_scale)
  duplicate_strategy <- match.arg(duplicate_strategy)
  signature_kcdf <- match.arg(signature_kcdf)

  annotation <- annotate_expr(
    expr = expr,
    meta = meta,
    species = match.arg(species),
    annotation_engine = match.arg(annotation_engine),
    output_file = file.path(output_dir, "expr_anno.csv"),
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
      method = match.arg(signature_method),
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

  files <- list.files(output_dir, full.names = TRUE)
  new_expranno_result(
    annotation = annotation,
    expr_meta_merged = expr_meta_merged,
    deconvolution = deconvolution,
    signatures = signatures,
    files = files
  )
}

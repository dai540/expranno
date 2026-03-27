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
#' @param run_signature Whether to run GSVA or ssGSEA.
#' @param geneset_file Optional GMT path.
#' @param gene_sets Optional named list of gene sets.
#' @param signature_method One of `"gsva"`, `"ssgsea"`, or `"both"`.
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
    run_signature = TRUE,
    geneset_file = NULL,
    gene_sets = NULL,
    signature_method = c("gsva", "ssgsea", "both"),
    verbose = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
    deconvolution <- run_cell_deconvolution(
      expr_anno = annotation$expr_anno,
      species = annotation$params$species,
      methods = deconv_methods,
      output_dir = output_dir,
      verbose = verbose
    )
  }

  signatures <- NULL
  if (isTRUE(run_signature)) {
    signatures <- run_signature_analysis(
      expr_anno = annotation$expr_anno,
      geneset_file = geneset_file,
      gene_sets = gene_sets,
      method = match.arg(signature_method),
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

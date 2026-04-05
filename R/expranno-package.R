#' expranno: annotation-first RNA-seq workflows
#'
#' `expranno` annotates human and mouse expression matrices, aligns sample
#' metadata, optionally runs immune deconvolution with `immunedeconv`, and
#' computes GSVA or ssGSEA signature scores from user-supplied gene sets.
#'
#' The package uses a coverage-first annotation strategy that can combine
#' `biomaRt`, `org.*.eg.db`, and optional `EnsDb` packages, and it returns
#' a compact annotation coverage report to support downstream quality checks.
#'
#' Built-in helpers such as `list_annotation_presets()`,
#' `example_annotation_truth()`, and `as_expranno_se()` make it easier to
#' standardize fixed human or mouse workflows, reproduce validation runs,
#' and move results back into Bioconductor containers.
#'
#' @importFrom stats aggregate
#' @importFrom stats setNames
#' @importFrom utils capture.output
#' @importFrom utils write.csv
#' @keywords internal
"_PACKAGE"

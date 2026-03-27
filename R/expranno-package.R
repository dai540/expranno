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
#' @importFrom stats aggregate
#' @importFrom utils write.csv
#' @keywords internal
"_PACKAGE"

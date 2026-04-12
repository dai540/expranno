#' expranno: minimal annotation-first RNA-seq workflows
#'
#' `expranno` provides a compact downstream workflow for human and mouse
#' expression matrices whose first column is `gene_id` and metadata tables whose
#' first column is `sample`.
#'
#' The package focuses on:
#'
#' - gene annotation
#' - expression and metadata integration
#' - optional deconvolution
#' - optional GSVA or ssGSEA signature scoring
#'
#' Heavy backends are kept in `Suggests` so the installed package stays small.
#'
#' @importFrom stats aggregate
#' @importFrom utils write.csv
#' @keywords internal
"_PACKAGE"


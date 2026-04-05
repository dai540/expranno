#' Example annotation truth tables
#'
#' Loads small bundled truth tables for human or mouse annotation validation.
#' These tables are intended for examples, tests, and reproducible benchmark
#' demonstrations with `validate_annotation_engines()`.
#'
#' @param species Either `"human"` or `"mouse"`.
#'
#' @return A `data.frame` with `gene_id`, `symbol`, `gene_name`, and `biotype`.
#' @export
example_annotation_truth <- function(species = c("human", "mouse")) {
  species <- match.arg(species)
  path <- system.file(
    "extdata",
    sprintf("%s_annotation_truth_v102.csv", species),
    package = "expranno"
  )
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

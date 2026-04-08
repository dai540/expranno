#' Example expression and metadata inputs
#'
#' Returns a small human or mouse RNA-seq example that follows the package
#' input contract. The toy matrix is suitable for examples, tests, and
#' documentation.
#'
#' @param species Either `"human"` or `"mouse"`.
#' @return A named list with `expr` and `meta`.
#' @export
example_expranno_data <- function(species = c("human", "mouse")) {
  species <- match.arg(species)

  expr <- switch(
    species,
    human = data.frame(
      gene_id = c("ENSG00000141510.17", "ENSG00000146648.18", "ENSG00000012048.23"),
      sample_a = c(120, 80, 25),
      sample_b = c(140, 77, 30),
      sample_c = c(118, 91, 21),
      stringsAsFactors = FALSE
    ),
    mouse = data.frame(
      gene_id = c("ENSMUSG00000059552.8", "ENSMUSG00000020122.15", "ENSMUSG00000017167.16"),
      sample_a = c(220, 150, 44),
      sample_b = c(205, 163, 47),
      sample_c = c(214, 158, 40),
      stringsAsFactors = FALSE
    )
  )

  meta <- data.frame(
    sample = c("sample_a", "sample_b", "sample_c"),
    group = c("case", "control", "case"),
    batch = c("b1", "b1", "b2"),
    species = species,
    stringsAsFactors = FALSE
  )

  list(expr = expr, meta = meta)
}

#' Example annotation truth tables
#'
#' Loads small bundled truth tables for human or mouse annotation validation.
#' These tables are intended for examples, tests, and reproducible benchmark
#' demonstrations with [validate_annotation_engines()].
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

#' List built-in annotation presets
#'
#' Returns the fixed annotation presets that ship with `expranno`. These presets
#' are intended to make repeated human and mouse annotation runs easier to
#' standardize across projects.
#'
#' @return A `data.frame` describing the available presets, including
#'   expression-scale intent, backend order, and bundled truth-table resource.
#' @export
list_annotation_presets <- function() {
  data.frame(
    annotation_preset = available_annotation_presets(),
    species = c("human", "mouse", "human", "mouse", "human", "mouse"),
    recommended_input = c("any", "any", "TPM-like abundance", "TPM-like abundance", "raw counts", "raw counts"),
    expr_scale = c("auto", "auto", "abundance", "abundance", "count", "count"),
    annotation_engine = rep("hybrid", 6L),
    strip_version = rep(TRUE, 6L),
    biomart_version = rep(102L, 6L),
    symbol_priority = c(
      "hgnc_symbol -> external_gene_name",
      "mgi_symbol -> external_gene_name",
      "hgnc_symbol -> external_gene_name",
      "mgi_symbol -> external_gene_name",
      "hgnc_symbol -> external_gene_name",
      "mgi_symbol -> external_gene_name"
    ),
    fallback_order = c(
      "biomaRt -> org.Hs.eg.db -> EnsDb.Hsapiens.v86",
      "biomaRt -> org.Mm.eg.db -> EnsDb.Mmusculus.v79",
      "biomaRt -> org.Hs.eg.db -> EnsDb.Hsapiens.v86",
      "biomaRt -> org.Mm.eg.db -> EnsDb.Mmusculus.v79",
      "biomaRt -> org.Hs.eg.db -> EnsDb.Hsapiens.v86",
      "biomaRt -> org.Mm.eg.db -> EnsDb.Mmusculus.v79"
    ),
    bundled_truth = c(
      "example_annotation_truth('human')",
      "example_annotation_truth('mouse')",
      "example_annotation_truth('human')",
      "example_annotation_truth('mouse')",
      "example_annotation_truth('human')",
      "example_annotation_truth('mouse')"
    ),
    stringsAsFactors = FALSE
  )
}

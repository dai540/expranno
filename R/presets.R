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

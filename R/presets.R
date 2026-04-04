#' List built-in annotation presets
#'
#' Returns the fixed annotation presets that ship with `expranno`. These presets
#' are intended to make repeated human and mouse annotation runs easier to
#' standardize across projects.
#'
#' @return A `data.frame` describing the available presets.
#' @export
list_annotation_presets <- function() {
  data.frame(
    annotation_preset = available_annotation_presets(),
    species = c("human", "mouse", "human", "mouse", "human", "mouse"),
    expr_scale = c("auto", "auto", "abundance", "abundance", "count", "count"),
    annotation_engine = rep("hybrid", 6L),
    strip_version = rep(TRUE, 6L),
    biomart_version = rep(102L, 6L),
    stringsAsFactors = FALSE
  )
}

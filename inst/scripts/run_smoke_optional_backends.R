.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

library(expranno)
library(immunedeconv)

make_real_like_expr_anno <- function(expr_mat) {
  data.frame(
    gene_id_raw = rownames(expr_mat),
    gene_id = rownames(expr_mat),
    symbol = rownames(expr_mat),
    expr_mat,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

normalize_output_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args) >= 1L) args[[1]] else file.path(getwd(), "smoke-outputs")
output_root <- normalize_output_dir(output_root)

human_dir <- file.path(output_root, "human-wrapper")
human_deconv_dir <- file.path(output_root, "human-deconvolution")
mouse_deconv_dir <- file.path(output_root, "mouse-deconvolution")
dir.create(human_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(human_deconv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mouse_deconv_dir, recursive = TRUE, showWarnings = FALSE)

human_demo <- example_expranno_data("human")
human_run <- run_expranno(
  expr = human_demo$expr,
  meta = human_demo$meta,
  species = "human",
  annotation_engine = "hybrid",
  expr_scale = "count",
  duplicate_strategy = "sum",
  output_dir = human_dir,
  run_deconvolution = FALSE,
  run_signature = TRUE,
  geneset_file = system.file("extdata", "hallmark_demo.gmt", package = "expranno"),
  signature_method = "both",
  signature_kcdf = "Poisson",
  verbose = FALSE
)

data("dataset_racle", package = "immunedeconv")
human_expr_anno <- make_real_like_expr_anno(dataset_racle$expr_mat)
human_deconv <- run_cell_deconvolution(
  human_expr_anno,
  species = "human",
  methods = c("epic", "mcp_counter", "quantiseq", "xcell"),
  expr_scale = "abundance",
  duplicate_strategy = "first",
  output_dir = human_deconv_dir,
  verbose = FALSE
)

data("dataset_petitprez", package = "immunedeconv")
mouse_expr_anno <- make_real_like_expr_anno(dataset_petitprez$expr_mat)
mouse_deconv <- run_cell_deconvolution(
  mouse_expr_anno,
  species = "mouse",
  methods = c("mmcp_counter", "base", "dcq"),
  expr_scale = "count",
  duplicate_strategy = "first",
  output_dir = mouse_deconv_dir,
  verbose = FALSE
)

summary_table <- data.frame(
  component = c(
    "human_annotation",
    "human_signature_gsva",
    "human_signature_ssgsea",
    "human_deconv_epic",
    "human_deconv_mcp_counter",
    "human_deconv_quantiseq",
    "human_deconv_xcell",
    "mouse_deconv_mmcp_counter",
    "mouse_deconv_base",
    "mouse_deconv_dcq"
  ),
  rows = c(
    nrow(human_run$annotation$expr_anno),
    nrow(human_run$signatures$gsva),
    nrow(human_run$signatures$ssgsea),
    nrow(human_deconv$epic),
    nrow(human_deconv$mcp_counter),
    nrow(human_deconv$quantiseq),
    nrow(human_deconv$xcell),
    nrow(mouse_deconv$mmcp_counter),
    nrow(mouse_deconv$base),
    nrow(mouse_deconv$dcq)
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(summary_table, file = file.path(output_root, "smoke_summary.csv"), row.names = FALSE)
print(summary_table)

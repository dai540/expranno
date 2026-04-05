.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

library(expranno)

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
human_truth <- example_annotation_truth("human")[, c("gene_id", "symbol"), drop = FALSE]

human_annotation <- annotate_expr(
  expr = human_demo$expr,
  meta = human_demo$meta,
  annotation_preset = "human_count_v102",
  output_file = file.path(human_dir, "expr_anno.csv"),
  report_file = file.path(human_dir, "annotation_report.csv"),
  ambiguity_file = file.path(human_dir, "annotation_ambiguity.csv"),
  provenance_file = file.path(human_dir, "annotation_provenance.csv"),
  verbose = FALSE
)

summary_components <- data.frame(
  component = c(
    "human_annotation",
    "human_merged"
  ),
  rows = c(
    nrow(human_annotation$expr_anno),
    nrow(
      merge_expr_meta(
        human_annotation$expr_anno,
        human_demo$meta,
        output_file = file.path(human_dir, "expr_meta_merged.csv")
      )
    )
  ),
  stringsAsFactors = FALSE
)

notes <- data.frame(
  component = character(),
  status = character(),
  detail = character(),
  stringsAsFactors = FALSE
)

validation_error <- NULL
validation <- tryCatch(
  validate_annotation_engines(
    expr = human_demo$expr,
    meta = human_demo$meta,
    truth = human_truth,
    species = "human",
    annotation_preset = "human_count_v102",
    engines = "hybrid",
    fields = "symbol",
    output_file = file.path(human_dir, "annotation_validation_summary.csv"),
    detail_file = file.path(human_dir, "annotation_validation_detail.csv"),
    verbose = FALSE
  ),
  error = function(e) {
    validation_error <<- conditionMessage(e)
    NULL
  }
)

if (is.null(validation)) {
  notes <- rbind(
    notes,
    data.frame(
      component = "validation",
      status = "error",
      detail = validation_error,
      stringsAsFactors = FALSE
    )
  )
} else {
  summary_components <- rbind(
    summary_components,
    data.frame(
      component = "human_validation_summary",
      rows = nrow(validation$summary),
      stringsAsFactors = FALSE
    )
  )
  notes <- rbind(
    notes,
    data.frame(
      component = "validation",
      status = "ran",
      detail = "annotation validation completed",
      stringsAsFactors = FALSE
    )
  )
}

signature_error <- NULL
has_signature_genes <- any(
  !is.na(human_annotation$expr_anno$symbol) &
    human_annotation$expr_anno$symbol != ""
)
signatures <- NULL

if (isTRUE(has_signature_genes)) {
  signatures <- tryCatch(
    run_signature_analysis(
      expr_anno = human_annotation$expr_anno,
      geneset_file = system.file("extdata", "hallmark_demo.gmt", package = "expranno"),
      method = "both",
      expr_scale = "count",
      duplicate_strategy = "sum",
      kcdf = "Poisson",
      output_dir = human_dir
    ),
    error = function(e) {
      signature_error <<- conditionMessage(e)
      NULL
    }
  )
}

if (!isTRUE(has_signature_genes)) {
  notes <- rbind(
    notes,
    data.frame(
      component = "signature",
      status = "skipped",
      detail = "no annotated symbols available for signature scoring",
      stringsAsFactors = FALSE
    )
  )
} else if (is.null(signatures)) {
  notes <- rbind(
    notes,
    data.frame(
      component = "signature",
      status = "error",
      detail = signature_error,
      stringsAsFactors = FALSE
    )
  )
} else {
  summary_components <- rbind(
    summary_components,
    data.frame(
      component = c("human_signature_gsva", "human_signature_ssgsea"),
      rows = c(nrow(signatures$gsva), nrow(signatures$ssgsea)),
      stringsAsFactors = FALSE
    )
  )
  notes <- rbind(
    notes,
    data.frame(
      component = "signature",
      status = "ran",
      detail = "signature scoring completed",
      stringsAsFactors = FALSE
    )
  )
}

writeLines(capture.output(sessionInfo()), con = file.path(human_dir, "session_info.txt"))

if (requireNamespace("immunedeconv", quietly = TRUE)) {
  deconv_error <- NULL
  tryCatch(
    {
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

      summary_components <- rbind(
        summary_components,
        data.frame(
          component = c(
            "human_deconv_epic",
            "human_deconv_mcp_counter",
            "human_deconv_quantiseq",
            "human_deconv_xcell",
            "mouse_deconv_mmcp_counter",
            "mouse_deconv_base",
            "mouse_deconv_dcq"
          ),
          rows = c(
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
      )

      notes <- data.frame(
        component = "deconvolution",
        status = "ran",
        detail = "immunedeconv installed",
        stringsAsFactors = FALSE
      )
    },
    error = function(e) {
      deconv_error <<- conditionMessage(e)
    }
  )

  if (!is.null(deconv_error)) {
    notes <- rbind(
      notes,
      data.frame(
        component = "deconvolution",
        status = "error",
        detail = deconv_error,
        stringsAsFactors = FALSE
      )
    )
  }
} else {
  notes <- rbind(
    notes,
    data.frame(
      component = "deconvolution",
      status = "skipped",
      detail = "immunedeconv not installed",
      stringsAsFactors = FALSE
    )
  )
}

utils::write.csv(summary_components, file = file.path(output_root, "smoke_summary.csv"), row.names = FALSE)
utils::write.csv(notes, file = file.path(output_root, "smoke_notes.csv"), row.names = FALSE)
print(summary_components)
print(notes)

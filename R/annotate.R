biomart_dataset <- function(species) {
  if (identical(species, "human")) "hsapiens_gene_ensembl" else "mmusculus_gene_ensembl"
}

biomart_attributes <- function(species) {
  if (identical(species, "human")) {
    c("ensembl_gene_id", "hgnc_symbol", "external_gene_name")
  } else {
    c("ensembl_gene_id", "mgi_symbol", "external_gene_name")
  }
}

annotate_with_biomart <- function(gene_ids, species, version) {
  require_namespace("biomaRt", reason = "biomaRt-based annotation")

  mart <- biomaRt::useEnsembl(
    biomart = "ensembl",
    dataset = biomart_dataset(species),
    version = version
  )
  raw <- biomaRt::getBM(
    attributes = biomart_attributes(species),
    filters = "ensembl_gene_id",
    values = unique(gene_ids),
    mart = mart
  )
  if (nrow(raw) == 0L) {
    return(data.frame(gene_id = character(), symbol = character(), gene_name = character(), stringsAsFactors = FALSE))
  }

  symbol_col <- if (identical(species, "human")) "hgnc_symbol" else "mgi_symbol"
  raw$symbol <- raw[[symbol_col]]
  raw$symbol[raw$symbol == ""] <- NA_character_
  raw$external_gene_name[raw$external_gene_name == ""] <- NA_character_
  raw$symbol <- ifelse(is.na(raw$symbol), raw$external_gene_name, raw$symbol)

  out <- raw[, c("ensembl_gene_id", "symbol", "external_gene_name"), drop = FALSE]
  names(out) <- c("gene_id", "symbol", "gene_name")
  out <- out[!duplicated(out$gene_id), , drop = FALSE]
  out
}

annotate_with_orgdb <- function(gene_ids, species) {
  require_namespace("AnnotationDbi", reason = "OrgDb-based annotation")

  org_pkg <- if (identical(species, "human")) "org.Hs.eg.db" else "org.Mm.eg.db"
  require_namespace(org_pkg, reason = "OrgDb-based annotation")
  orgdb <- getExportedValue(org_pkg, org_pkg)

  raw <- AnnotationDbi::select(
    x = orgdb,
    keys = unique(gene_ids),
    keytype = "ENSEMBL",
    columns = c("SYMBOL", "GENENAME")
  )
  if (nrow(raw) == 0L) {
    return(data.frame(gene_id = character(), symbol = character(), gene_name = character(), stringsAsFactors = FALSE))
  }

  out <- raw[, c("ENSEMBL", "SYMBOL", "GENENAME"), drop = FALSE]
  names(out) <- c("gene_id", "symbol", "gene_name")
  out <- out[!duplicated(out$gene_id), , drop = FALSE]
  out
}

combine_annotation_sources <- function(gene_ids, biomart = NULL, orgdb = NULL) {
  out <- data.frame(
    gene_id = unique(gene_ids),
    symbol = NA_character_,
    gene_name = NA_character_,
    annotation_source = NA_character_,
    stringsAsFactors = FALSE
  )

  if (!is.null(biomart) && nrow(biomart) > 0L) {
    idx <- match(out$gene_id, biomart$gene_id)
    has <- !is.na(idx)
    out$symbol[has] <- biomart$symbol[idx[has]]
    out$gene_name[has] <- biomart$gene_name[idx[has]]
    out$annotation_source[has] <- "biomart"
  }

  if (!is.null(orgdb) && nrow(orgdb) > 0L) {
    idx <- match(out$gene_id, orgdb$gene_id)
    has <- !is.na(idx)
    fill_symbol <- has & (is.na(out$symbol) | !nzchar(out$symbol))
    fill_name <- has & (is.na(out$gene_name) | !nzchar(out$gene_name))
    out$symbol[fill_symbol] <- orgdb$symbol[idx[fill_symbol]]
    out$gene_name[fill_name] <- orgdb$gene_name[idx[fill_name]]
    out$annotation_source[fill_symbol | fill_name] <- ifelse(
      is.na(out$annotation_source[fill_symbol | fill_name]),
      "orgdb",
      "hybrid"
    )
  }

  out$annotation_status <- ifelse(
    (!is.na(out$symbol) & nzchar(out$symbol)) | (!is.na(out$gene_name) & nzchar(out$gene_name)),
    "annotated",
    "unmapped"
  )
  out
}

annotation_report <- function(expr_anno, species, engine, version) {
  annotated <- sum(expr_anno$annotation_status == "annotated")
  total <- nrow(expr_anno)
  data.frame(
    species = species,
    annotation_engine = engine,
    biomart_version = if (identical(engine, "none")) NA_integer_ else version,
    rows = total,
    unique_gene_ids = length(unique(expr_anno$gene_id)),
    annotated_rows = annotated,
    unannotated_rows = total - annotated,
    annotation_rate = if (total == 0L) NA_real_ else annotated / total,
    stringsAsFactors = FALSE
  )
}

#' Annotate an expression matrix
#'
#' @param expr Expression table with `gene_id` in the first column.
#' @param species `"auto"`, `"human"`, or `"mouse"`.
#' @param annotation_engine `"hybrid"`, `"biomart"`, `"orgdb"`, or `"none"`.
#' @param strip_version Whether to remove Ensembl version suffixes.
#' @param biomart_version Ensembl release used when `biomaRt` is active.
#' @param output_dir Optional output directory.
#'
#' @return A list with `expr_anno`, `report`, and `species`.
#' @export
annotate_expr <- function(
    expr,
    species = c("auto", "human", "mouse"),
    annotation_engine = c("hybrid", "biomart", "orgdb", "none"),
    strip_version = TRUE,
    biomart_version = 102,
    output_dir = NULL) {
  validate_expr(expr)
  species <- match.arg(species)
  annotation_engine <- match.arg(annotation_engine)

  gene_id_raw <- as.character(expr[[1]])
  gene_id <- if (isTRUE(strip_version)) sub("\\.\\d+$", "", gene_id_raw) else gene_id_raw
  species <- if (identical(species, "auto")) detect_species(gene_id) else species

  biomart <- NULL
  orgdb <- NULL
  if (annotation_engine %in% c("hybrid", "biomart")) {
    biomart <- annotate_with_biomart(gene_id, species = species, version = biomart_version)
  }
  if (annotation_engine %in% c("hybrid", "orgdb")) {
    orgdb <- annotate_with_orgdb(gene_id, species = species)
  }

  annotation <- if (identical(annotation_engine, "none")) {
    data.frame(
      gene_id = unique(gene_id),
      symbol = NA_character_,
      gene_name = NA_character_,
      annotation_source = "none",
      annotation_status = "unmapped",
      stringsAsFactors = FALSE
    )
  } else {
    combine_annotation_sources(gene_id, biomart = biomart, orgdb = orgdb)
  }

  idx <- match(gene_id, annotation$gene_id)
  expr_anno <- data.frame(
    gene_id_raw = gene_id_raw,
    gene_id = gene_id,
    symbol = annotation$symbol[idx],
    gene_name = annotation$gene_name[idx],
    annotation_source = annotation$annotation_source[idx],
    annotation_status = annotation$annotation_status[idx],
    expr[, -1, drop = FALSE],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  report <- annotation_report(
    expr_anno = expr_anno,
    species = species,
    engine = annotation_engine,
    version = biomart_version
  )

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write_if_requested(expr_anno, file.path(output_dir, "expr_anno.csv"))
    write_if_requested(report, file.path(output_dir, "annotation_report.csv"))
  }

  structure(
    list(expr_anno = expr_anno, report = report, species = species),
    class = "expranno_annotation"
  )
}


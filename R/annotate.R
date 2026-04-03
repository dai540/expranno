query_biomart_annotation <- function(gene_ids, species) {
  require_namespace("biomaRt", reason = "hybrid annotation")

  dataset <- switch(
    species,
    human = "hsapiens_gene_ensembl",
    mouse = "mmusculus_gene_ensembl",
    stop("Unsupported species for biomaRt.", call. = FALSE)
  )

  mart <- biomaRt::useEnsembl(biomart = "genes", dataset = dataset)
  raw <- biomaRt::getBM(
    attributes = c(
      "ensembl_gene_id",
      "external_gene_name",
      "description",
      "entrezgene_id",
      "gene_biotype",
      "chromosome_name",
      "start_position",
      "end_position",
      "strand"
    ),
    filters = "ensembl_gene_id",
    values = unique(gene_ids),
    mart = mart
  )

  out <- data.frame(
    gene_id = raw$ensembl_gene_id,
    symbol = raw$external_gene_name,
    gene_name = raw$description,
    entrez_id = as.character(raw$entrezgene_id),
    biotype = raw$gene_biotype,
    chromosome = raw$chromosome_name,
    start = raw$start_position,
    end = raw$end_position,
    strand = raw$strand,
    stringsAsFactors = FALSE
  )
  out[!duplicated(out$gene_id), , drop = FALSE]
}

query_orgdb_annotation <- function(gene_ids, species) {
  require_namespace("AnnotationDbi", reason = "OrgDb annotation")
  org_pkg <- switch(
    species,
    human = "org.Hs.eg.db",
    mouse = "org.Mm.eg.db",
    stop("Unsupported species for OrgDb.", call. = FALSE)
  )
  require_namespace(org_pkg, reason = "OrgDb annotation")

  db <- getExportedValue(org_pkg, org_pkg)

  data.frame(
    gene_id = gene_ids,
    symbol = unname(AnnotationDbi::mapIds(db, keys = gene_ids, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")),
    gene_name = unname(AnnotationDbi::mapIds(db, keys = gene_ids, keytype = "ENSEMBL", column = "GENENAME", multiVals = "first")),
    entrez_id = unname(AnnotationDbi::mapIds(db, keys = gene_ids, keytype = "ENSEMBL", column = "ENTREZID", multiVals = "first")),
    stringsAsFactors = FALSE
  )
}

query_ensdb_annotation <- function(gene_ids, species) {
  require_namespace("ensembldb", reason = "EnsDb annotation")
  require_namespace("AnnotationFilter", reason = "EnsDb annotation")

  ens_pkg <- switch(
    species,
    human = "EnsDb.Hsapiens.v86",
    mouse = "EnsDb.Mmusculus.v79",
    stop("Unsupported species for EnsDb.", call. = FALSE)
  )
  require_namespace(ens_pkg, reason = "EnsDb annotation")

  db_name <- switch(
    species,
    human = "EnsDb.Hsapiens.v86",
    mouse = "EnsDb.Mmusculus.v79"
  )
  db <- getExportedValue(ens_pkg, db_name)
  raw <- ensembldb::genes(
    db,
    filter = AnnotationFilter::GeneIdFilter(unique(gene_ids)),
    return.type = "DataFrame"
  )

  out <- as.data.frame(raw)
  rename_map <- c(
    gene_id = "gene_id",
    symbol = "gene_name",
    gene_name = "description",
    biotype = "gene_biotype",
    chromosome = "seq_name",
    start = "gene_seq_start",
    end = "gene_seq_end",
    strand = "seq_strand"
  )
  keep <- intersect(unname(rename_map), names(out))
  out <- out[, keep, drop = FALSE]
  names(out) <- names(rename_map)[match(names(out), rename_map)]
  out[!duplicated(out$gene_id), , drop = FALSE]
}

build_annotation <- function(gene_info, species, annotation_engine = c("hybrid", "biomart", "orgdb", "ensdb", "none"), fields = NULL, verbose = TRUE) {
  annotation_engine <- match.arg(annotation_engine)
  if (is.null(fields)) {
    fields <- default_annotation_fields()
  }

  base <- gene_info
  for (field in fields) {
    base[[field]] <- NA
  }
  base$annotation_source <- NA_character_
  base$annotation_status <- "unannotated"

  if (annotation_engine == "none") {
    return(base)
  }

  sources <- switch(
    annotation_engine,
    hybrid = c("biomart", "orgdb", "ensdb"),
    biomart = "biomart",
    orgdb = "orgdb",
    ensdb = "ensdb"
  )

  source_results <- list()
  for (source in sources) {
    source_results[[source]] <- tryCatch(
      {
        if (verbose) {
          message(sprintf("Querying %s annotation backend...", source))
        }
        switch(
          source,
          biomart = query_biomart_annotation(base$gene_id, species),
          orgdb = query_orgdb_annotation(base$gene_id, species),
          ensdb = query_ensdb_annotation(base$gene_id, species)
        )
      },
      error = function(e) {
        if (verbose) {
          message(sprintf("Skipping %s backend: %s", source, conditionMessage(e)))
        }
        NULL
      }
    )
  }

  for (source in names(source_results)) {
    incoming <- source_results[[source]]
    if (is.null(incoming)) {
      next
    }
    for (field in fields) {
      before <- base[[field]]
      base <- merge_fill(base, incoming, field)
      after <- base[[field]]
      filled <- (is.na(before) | before == "") & (!is.na(after) & after != "")
      source_col <- paste0(field, "_source")
      if (!source_col %in% names(base)) {
        base[[source_col]] <- NA_character_
      }
      base[[source_col]][filled] <- source
    }
  }

  any_annotated <- Reduce(
    `|`,
    lapply(fields, function(field) !is.na(base[[field]]) & base[[field]] != "")
  )
  base$annotation_status[any_annotated] <- "annotated"
  base$annotation_source[any_annotated] <- vapply(
    which(any_annotated),
    function(i) {
      filled_sources <- unlist(base[i, grep("_source$", names(base)), drop = FALSE], use.names = FALSE)
      filled_sources <- unique(filled_sources[!is.na(filled_sources)])
      if (length(filled_sources) == 0L) NA_character_ else paste(filled_sources, collapse = ";")
    },
    character(1)
  )

  base
}

#' Annotate an expression matrix
#'
#' Validates `expr` and `meta`, normalizes Ensembl IDs, infers or accepts the
#' species, runs a coverage-first annotation workflow, binds annotation columns
#' to the expression table, and optionally writes `expr_anno.csv`.
#'
#' @param expr Expression table with `gene_id` in the first column.
#' @param meta Metadata table with `sample` in the first column.
#' @param species Either `"auto"`, `"human"`, or `"mouse"`.
#' @param annotation_engine Annotation backend strategy. Use `"hybrid"` to try
#'   `biomaRt`, `orgdb`, and `EnsDb` in sequence.
#' @param fields Optional annotation fields to keep.
#' @param strip_version Whether to remove Ensembl version suffixes.
#' @param output_file Optional CSV path for the annotated matrix.
#' @param verbose Whether to emit progress messages.
#'
#' @return An `expranno_annotation` object.
#' @export
annotate_expr <- function(
    expr,
    meta,
    species = c("auto", "human", "mouse"),
    annotation_engine = c("hybrid", "biomart", "orgdb", "ensdb", "none"),
    fields = NULL,
    strip_version = TRUE,
    output_file = NULL,
    verbose = TRUE) {
  validate_expr(expr)
  validate_meta(meta)

  gene_info <- normalize_gene_ids(expr$gene_id, strip_version = strip_version)
  resolved_species <- resolve_species(species, gene_info$gene_id)
  meta_checked <- match_meta_to_expr(expr, meta)

  annotation <- build_annotation(
    gene_info = gene_info,
    species = resolved_species,
    annotation_engine = match.arg(annotation_engine),
    fields = fields,
    verbose = verbose
  )

  expr_anno <- cbind(annotation, expr[, -1, drop = FALSE], stringsAsFactors = FALSE)
  report <- annotation_coverage_report(annotation, fields %||% default_annotation_fields())

  write_if_requested(expr_anno, output_file = output_file)

  new_expranno_annotation(
    expr_anno = expr_anno,
    annotation = annotation,
    meta_checked = meta_checked,
    report = report,
    params = list(
      species = resolved_species,
      annotation_engine = annotation_engine,
      strip_version = strip_version
    )
  )
}

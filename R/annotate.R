biomart_dataset_for_species <- function(species) {
  switch(
    species,
    human = "hsapiens_gene_ensembl",
    mouse = "mmusculus_gene_ensembl",
    stop("Unsupported species for biomaRt.", call. = FALSE)
  )
}

biomart_symbol_attribute_for_species <- function(species) {
  switch(
    species,
    human = "hgnc_symbol",
    mouse = "mgi_symbol",
    stop("Unsupported species for biomaRt.", call. = FALSE)
  )
}

orgdb_package_for_species <- function(species) {
  switch(
    species,
    human = "org.Hs.eg.db",
    mouse = "org.Mm.eg.db",
    stop("Unsupported species for OrgDb.", call. = FALSE)
  )
}

ensdb_package_for_species <- function(species) {
  switch(
    species,
    human = "EnsDb.Hsapiens.v86",
    mouse = "EnsDb.Mmusculus.v79",
    stop("Unsupported species for EnsDb.", call. = FALSE)
  )
}

empty_backend_annotation <- function(gene_ids, fields, source) {
  out <- data.frame(gene_id = unique(gene_ids), stringsAsFactors = FALSE)
  for (field in fields) {
    out[[field]] <- NA_character_
    out[[paste0(field, "_candidates")]] <- NA_character_
    out[[paste0(field, "_candidate_count")]] <- 0L
    out[[paste0(field, "_is_ambiguous")]] <- FALSE
    out[[paste0(field, "_source")]] <- NA_character_
  }
  out
}

new_backend_provenance <- function(
    source,
    species,
    status,
    genes_queried,
    rows_returned = 0L,
    dataset = NA_character_,
    package_name = NA_character_,
    backend_release = NA_character_,
    host = NA_character_,
    mirror = NA_character_,
    message = NA_character_) {
  data.frame(
    source = source,
    species = species,
    dataset = dataset,
    package_name = package_name,
    backend_release = backend_release,
    host = host,
    mirror = mirror,
    genes_queried = as.integer(genes_queried),
    rows_returned = as.integer(rows_returned),
    annotation_date = as.character(Sys.Date()),
    status = status,
    message = message,
    stringsAsFactors = FALSE
  )
}

collapse_backend_records <- function(raw, gene_ids, field_map, source) {
  fields <- names(field_map)
  out <- data.frame(gene_id = unique(gene_ids), stringsAsFactors = FALSE)

  if (nrow(raw) == 0L) {
    return(empty_backend_annotation(gene_ids = gene_ids, fields = fields, source = source))
  }

  for (field in fields) {
    column <- field_map[[field]]
    values_by_gene <- split(raw[[column]], raw$gene_id)
    collapsed <- lapply(out$gene_id, function(gene_id) {
      collapse_annotation_field(values_by_gene[[gene_id]] %||% character(0))
    })

    out[[field]] <- vapply(collapsed, `[[`, character(1), "value")
    out[[paste0(field, "_candidates")]] <- vapply(collapsed, `[[`, character(1), "candidates")
    out[[paste0(field, "_candidate_count")]] <- vapply(collapsed, `[[`, integer(1), "candidate_count")
    out[[paste0(field, "_is_ambiguous")]] <- vapply(collapsed, `[[`, logical(1), "is_ambiguous")
    out[[paste0(field, "_source")]] <- ifelse(
      !is.na(out[[field]]) & nzchar(out[[field]]),
      source,
      NA_character_
    )
  }

  out
}

query_biomart_annotation <- function(
    gene_ids,
    species,
    fields,
    biomart_version = 102,
    biomart_host = NULL,
    biomart_mirror = NULL) {
  require_namespace("biomaRt", reason = "hybrid annotation")

  dataset <- biomart_dataset_for_species(species)
  symbol_attr <- biomart_symbol_attribute_for_species(species)
  args <- list(
    biomart = "genes",
    dataset = dataset
  )
  if (!is.null(biomart_host) && nzchar(biomart_host)) {
    args$host <- biomart_host
  } else if (!is.null(biomart_version) && !is.na(biomart_version)) {
    args$version <- biomart_version
  } else if (!is.null(biomart_mirror) && nzchar(biomart_mirror)) {
    args$mirror <- biomart_mirror
  }

  mart <- do.call(biomaRt::useEnsembl, args)
  raw <- biomaRt::getBM(
    attributes = unique(c(
      "ensembl_gene_id",
      symbol_attr,
      "external_gene_name",
      "description",
      "entrezgene_id",
      "gene_biotype",
      "chromosome_name",
      "start_position",
      "end_position",
      "strand"
    )),
    filters = "ensembl_gene_id",
    values = unique(gene_ids),
    mart = mart
  )

  raw$symbol <- clean_annotation_value(raw[[symbol_attr]])
  fallback_symbol <- clean_annotation_value(raw$external_gene_name)
  take_fallback <- is.na(raw$symbol)
  raw$symbol[take_fallback] <- fallback_symbol[take_fallback]
  raw$gene_name <- clean_annotation_value(sub(" \\[Source:.*$", "", raw$description))
  raw$entrez_id <- clean_annotation_value(raw$entrezgene_id)
  raw$biotype <- clean_annotation_value(raw$gene_biotype)
  raw$chromosome <- clean_annotation_value(raw$chromosome_name)
  raw$start <- clean_annotation_value(raw$start_position)
  raw$end <- clean_annotation_value(raw$end_position)
  raw$strand <- clean_annotation_value(raw$strand)

  collapsed <- collapse_backend_records(
    raw = data.frame(
      gene_id = raw$ensembl_gene_id,
      symbol = raw$symbol,
      gene_name = raw$gene_name,
      entrez_id = raw$entrez_id,
      biotype = raw$biotype,
      chromosome = raw$chromosome,
      start = raw$start,
      end = raw$end,
      strand = raw$strand,
      stringsAsFactors = FALSE
    ),
    gene_ids = gene_ids,
    field_map = stats::setNames(fields, fields),
    source = "biomart"
  )

  list(
    data = collapsed,
    provenance = new_backend_provenance(
      source = "biomart",
      species = species,
      status = "ok",
      genes_queried = length(unique(gene_ids)),
      rows_returned = nrow(raw),
      dataset = dataset,
      package_name = "biomaRt",
      backend_release = if (!is.null(biomart_version) && !is.na(biomart_version)) {
        sprintf("Ensembl v%s", biomart_version)
      } else {
        "Ensembl current"
      },
      host = biomart_host %||% NA_character_,
      mirror = biomart_mirror %||% NA_character_
    )
  )
}

query_orgdb_annotation <- function(gene_ids, species, fields) {
  require_namespace("AnnotationDbi", reason = "OrgDb annotation")
  org_pkg <- orgdb_package_for_species(species)
  require_namespace(org_pkg, reason = "OrgDb annotation")

  db <- getExportedValue(org_pkg, org_pkg)
  raw <- AnnotationDbi::select(
    db,
    keys = unique(gene_ids),
    keytype = "ENSEMBL",
    columns = c("SYMBOL", "GENENAME", "ENTREZID")
  )

  collapsed <- collapse_backend_records(
    raw = data.frame(
      gene_id = raw$ENSEMBL,
      symbol = clean_annotation_value(raw$SYMBOL),
      gene_name = clean_annotation_value(raw$GENENAME),
      entrez_id = clean_annotation_value(raw$ENTREZID),
      stringsAsFactors = FALSE
    ),
    gene_ids = gene_ids,
    field_map = c(symbol = "symbol", gene_name = "gene_name", entrez_id = "entrez_id"),
    source = "orgdb"
  )
  for (field in setdiff(fields, c("symbol", "gene_name", "entrez_id"))) {
    collapsed[[field]] <- NA_character_
    collapsed[[paste0(field, "_candidates")]] <- NA_character_
    collapsed[[paste0(field, "_candidate_count")]] <- 0L
    collapsed[[paste0(field, "_is_ambiguous")]] <- FALSE
    collapsed[[paste0(field, "_source")]] <- NA_character_
  }
  collapsed <- collapsed[, c("gene_id", unlist(lapply(fields, function(field) {
    c(
      field,
      paste0(field, "_candidates"),
      paste0(field, "_candidate_count"),
      paste0(field, "_is_ambiguous"),
      paste0(field, "_source")
    )
  }))), drop = FALSE]

  list(
    data = collapsed,
    provenance = new_backend_provenance(
      source = "orgdb",
      species = species,
      status = "ok",
      genes_queried = length(unique(gene_ids)),
      rows_returned = nrow(raw),
      dataset = org_pkg,
      package_name = org_pkg,
      backend_release = sprintf("%s %s", org_pkg, as.character(utils::packageVersion(org_pkg)))
    )
  )
}

query_ensdb_annotation <- function(gene_ids, species, fields) {
  require_namespace("ensembldb", reason = "EnsDb annotation")
  require_namespace("AnnotationFilter", reason = "EnsDb annotation")

  ens_pkg <- ensdb_package_for_species(species)
  require_namespace(ens_pkg, reason = "EnsDb annotation")

  db <- getExportedValue(ens_pkg, ens_pkg)
  raw <- ensembldb::genes(
    db,
    filter = AnnotationFilter::GeneIdFilter(unique(gene_ids)),
    return.type = "DataFrame"
  )
  raw <- as.data.frame(raw)

  collapsed <- collapse_backend_records(
    raw = data.frame(
      gene_id = clean_annotation_value(raw$gene_id),
      symbol = clean_annotation_value(raw$gene_name),
      gene_name = clean_annotation_value(raw$description),
      biotype = clean_annotation_value(raw$gene_biotype),
      chromosome = clean_annotation_value(raw$seq_name),
      start = clean_annotation_value(raw$gene_seq_start),
      end = clean_annotation_value(raw$gene_seq_end),
      strand = clean_annotation_value(raw$seq_strand),
      stringsAsFactors = FALSE
    ),
    gene_ids = gene_ids,
    field_map = c(
      symbol = "symbol",
      gene_name = "gene_name",
      biotype = "biotype",
      chromosome = "chromosome",
      start = "start",
      end = "end",
      strand = "strand"
    ),
    source = "ensdb"
  )
  if (!"entrez_id" %in% fields) {
    collapsed$entrez_id <- NULL
  }
  for (field in setdiff(fields, c("symbol", "gene_name", "biotype", "chromosome", "start", "end", "strand"))) {
    collapsed[[field]] <- NA_character_
    collapsed[[paste0(field, "_candidates")]] <- NA_character_
    collapsed[[paste0(field, "_candidate_count")]] <- 0L
    collapsed[[paste0(field, "_is_ambiguous")]] <- FALSE
    collapsed[[paste0(field, "_source")]] <- NA_character_
  }
  collapsed <- collapsed[, c("gene_id", unlist(lapply(fields, function(field) {
    c(
      field,
      paste0(field, "_candidates"),
      paste0(field, "_candidate_count"),
      paste0(field, "_is_ambiguous"),
      paste0(field, "_source")
    )
  }))), drop = FALSE]

  list(
    data = collapsed,
    provenance = new_backend_provenance(
      source = "ensdb",
      species = species,
      status = "ok",
      genes_queried = length(unique(gene_ids)),
      rows_returned = nrow(raw),
      dataset = ens_pkg,
      package_name = ens_pkg,
      backend_release = sprintf("%s %s", ens_pkg, as.character(utils::packageVersion(ens_pkg)))
    )
  )
}

annotation_source_order <- function(annotation_engine, fields) {
  if (annotation_engine == "hybrid") {
    orders <- list(
      symbol = c("biomart", "orgdb", "ensdb"),
      gene_name = c("biomart", "orgdb", "ensdb"),
      entrez_id = c("biomart", "orgdb", "ensdb"),
      biotype = c("biomart", "ensdb", "orgdb"),
      chromosome = c("biomart", "ensdb", "orgdb"),
      start = c("biomart", "ensdb", "orgdb"),
      end = c("biomart", "ensdb", "orgdb"),
      strand = c("biomart", "ensdb", "orgdb")
    )
    return(orders[fields])
  }

  setNames(rep(list(annotation_engine), length(fields)), fields)
}

backend_metadata_from_sources <- function(provenance, sources, field) {
  if (length(sources) == 0L) {
    return(NA_character_)
  }
  values <- provenance[[field]][match(sources, provenance$source)]
  values <- clean_annotation_value(values)
  values <- unique(values[!is.na(values)])
  if (length(values) == 0L) {
    NA_character_
  } else {
    paste(values, collapse = ";")
  }
}

build_annotation <- function(
    gene_info,
    species,
    annotation_engine = c("hybrid", "biomart", "orgdb", "ensdb", "none"),
    fields = NULL,
    verbose = TRUE,
    biomart_version = 102,
    biomart_host = NULL,
    biomart_mirror = NULL) {
  annotation_engine <- match.arg(annotation_engine)
  if (is.null(fields)) {
    fields <- default_annotation_fields()
  }

  base <- gene_info
  for (field in fields) {
    base[[field]] <- NA_character_
    base[[paste0(field, "_candidates")]] <- NA_character_
    base[[paste0(field, "_candidate_count")]] <- 0L
    base[[paste0(field, "_is_ambiguous")]] <- FALSE
    base[[paste0(field, "_source")]] <- NA_character_
  }
  base$annotation_source <- NA_character_
  base$annotation_status <- "unannotated"
  base$annotation_backend_release <- NA_character_
  base$annotation_backend_host <- NA_character_
  base$annotation_backend_mirror <- NA_character_
  base$annotation_date <- as.character(Sys.Date())

  if (annotation_engine == "none") {
    provenance <- new_backend_provenance(
      source = "none",
      species = species,
      status = "skipped",
      genes_queried = nrow(base),
      rows_returned = 0L,
      dataset = "none",
      package_name = "none",
      message = "Annotation lookup skipped by user request."
    )
    return(list(annotation = base, provenance = provenance))
  }

  sources <- switch(
    annotation_engine,
    hybrid = c("biomart", "orgdb", "ensdb"),
    biomart = "biomart",
    orgdb = "orgdb",
    ensdb = "ensdb"
  )

  source_results <- list()
  provenance_rows <- list()
  for (source in sources) {
    source_result <- tryCatch(
      {
        if (isTRUE(verbose)) {
          message(sprintf("Querying %s annotation backend...", source))
        }
        switch(
          source,
          biomart = query_biomart_annotation(
            gene_ids = base$gene_id,
            species = species,
            fields = fields,
            biomart_version = biomart_version,
            biomart_host = biomart_host,
            biomart_mirror = biomart_mirror
          ),
          orgdb = query_orgdb_annotation(
            gene_ids = base$gene_id,
            species = species,
            fields = fields
          ),
          ensdb = query_ensdb_annotation(
            gene_ids = base$gene_id,
            species = species,
            fields = fields
          )
        )
      },
      error = function(e) {
        if (isTRUE(verbose)) {
          message(sprintf("Skipping %s backend: %s", source, conditionMessage(e)))
        }
        list(
          data = empty_backend_annotation(gene_ids = base$gene_id, fields = fields, source = source),
          provenance = new_backend_provenance(
            source = source,
            species = species,
            status = "error",
            genes_queried = length(unique(base$gene_id)),
            rows_returned = 0L,
            dataset = switch(
              source,
              biomart = biomart_dataset_for_species(species),
              orgdb = orgdb_package_for_species(species),
              ensdb = ensdb_package_for_species(species)
            ),
            package_name = switch(
              source,
              biomart = "biomaRt",
              orgdb = orgdb_package_for_species(species),
              ensdb = ensdb_package_for_species(species)
            ),
            backend_release = switch(
              source,
              biomart = if (!is.null(biomart_version) && !is.na(biomart_version)) {
                sprintf("Ensembl v%s", biomart_version)
              } else {
                "Ensembl current"
              },
              orgdb = orgdb_package_for_species(species),
              ensdb = ensdb_package_for_species(species)
            ),
            host = biomart_host %||% NA_character_,
            mirror = biomart_mirror %||% NA_character_,
            message = conditionMessage(e)
          )
        )
      }
    )

    source_results[[source]] <- source_result$data
    provenance_rows[[source]] <- source_result$provenance
  }

  field_orders <- annotation_source_order(annotation_engine = annotation_engine, fields = fields)

  for (field in fields) {
    candidates_by_source <- list()
    chosen_value <- rep(NA_character_, nrow(base))
    chosen_source <- rep(NA_character_, nrow(base))

    for (source in field_orders[[field]]) {
      incoming <- source_results[[source]]
      idx <- match(base$gene_id, incoming$gene_id)
      field_values <- clean_annotation_value(incoming[[field]][idx])
      candidates_by_source[[source]] <- incoming[[paste0(field, "_candidates")]][idx]

      fill <- is.na(chosen_value) & !is.na(field_values)
      chosen_value[fill] <- field_values[fill]
      chosen_source[fill] <- source
    }

    candidate_summary <- lapply(seq_len(nrow(base)), function(i) {
      collapse_candidate_strings(vapply(candidates_by_source, `[`, character(1), i))
    })

    base[[field]] <- chosen_value
    base[[paste0(field, "_source")]] <- chosen_source
    base[[paste0(field, "_candidates")]] <- vapply(candidate_summary, `[[`, character(1), "candidates")
    base[[paste0(field, "_candidate_count")]] <- vapply(candidate_summary, `[[`, integer(1), "candidate_count")
    base[[paste0(field, "_is_ambiguous")]] <- vapply(candidate_summary, `[[`, logical(1), "is_ambiguous")
  }

  any_annotated <- Reduce(
    `|`,
    lapply(fields, function(field) !is.na(base[[field]]) & nzchar(base[[field]]))
  )
  source_cols <- paste0(fields, "_source")

  base$annotation_status[any_annotated] <- "annotated"
  base$annotation_source[any_annotated] <- vapply(
    which(any_annotated),
    function(i) {
      sources_i <- unique(unlist(base[i, source_cols, drop = FALSE], use.names = FALSE))
      sources_i <- sources_i[!is.na(sources_i)]
      if (length(sources_i) == 0L) {
        NA_character_
      } else {
        paste(sources_i, collapse = ";")
      }
    },
    character(1)
  )

  provenance <- do.call(rbind, provenance_rows)
  base$annotation_backend_release <- vapply(seq_len(nrow(base)), function(i) {
    backend_metadata_from_sources(
      provenance = provenance,
      sources = candidate_values_from_string(base$annotation_source[[i]]),
      field = "backend_release"
    )
  }, character(1))
  base$annotation_backend_host <- vapply(seq_len(nrow(base)), function(i) {
    backend_metadata_from_sources(
      provenance = provenance,
      sources = candidate_values_from_string(base$annotation_source[[i]]),
      field = "host"
    )
  }, character(1))
  base$annotation_backend_mirror <- vapply(seq_len(nrow(base)), function(i) {
    backend_metadata_from_sources(
      provenance = provenance,
      sources = candidate_values_from_string(base$annotation_source[[i]]),
      field = "mirror"
    )
  }, character(1))
  base$annotation_date <- vapply(seq_len(nrow(base)), function(i) {
    backend_metadata_from_sources(
      provenance = provenance,
      sources = candidate_values_from_string(base$annotation_source[[i]]),
      field = "annotation_date"
    )
  }, character(1))

  list(annotation = base, provenance = provenance)
}

#' Annotate an expression matrix
#'
#' Validates input, normalizes Ensembl IDs, infers or accepts the species, runs
#' a coverage-first annotation workflow, binds annotation columns to the
#' expression table, and optionally writes `expr_anno.csv` plus provenance and
#' ambiguity reports.
#'
#' @param expr Expression table with `gene_id` in the first column, or a
#'   `SummarizedExperiment`-like object.
#' @param meta Metadata table with `sample` in the first column. Leave `NULL`
#'   when `expr` is a `SummarizedExperiment` and metadata should be taken from
#'   `colData`.
#' @param species Either `"auto"`, `"human"`, or `"mouse"`.
#' @param annotation_engine Annotation backend strategy. Use `"hybrid"` to try
#'   `biomaRt`, `orgdb`, and `EnsDb` in sequence.
#' @param fields Optional annotation fields to keep.
#' @param strip_version Whether to remove Ensembl version suffixes.
#' @param biomart_version Fixed Ensembl release used by the `biomaRt` backend.
#' @param biomart_host Optional explicit Ensembl host.
#' @param biomart_mirror Optional Ensembl mirror name.
#' @param assay_name Optional assay name when `expr` is a
#'   `SummarizedExperiment`.
#' @param gene_id_col Optional row-data column containing Ensembl IDs when
#'   `expr` is a `SummarizedExperiment`.
#' @param sample_col Metadata column to use as `sample` when `expr` is a
#'   `SummarizedExperiment`.
#' @param output_file Optional CSV path for the annotated matrix.
#' @param report_file Optional CSV path for the coverage report.
#' @param ambiguity_file Optional CSV path for the ambiguity report.
#' @param provenance_file Optional CSV path for backend provenance.
#' @param verbose Whether to emit progress messages.
#'
#' @return An `expranno_annotation` object.
#' @export
annotate_expr <- function(
    expr,
    meta = NULL,
    species = c("auto", "human", "mouse"),
    annotation_engine = c("hybrid", "biomart", "orgdb", "ensdb", "none"),
    fields = NULL,
    strip_version = TRUE,
    biomart_version = 102,
    biomart_host = NULL,
    biomart_mirror = NULL,
    assay_name = NULL,
    gene_id_col = NULL,
    sample_col = "sample",
    output_file = NULL,
    report_file = NULL,
    ambiguity_file = NULL,
    provenance_file = NULL,
    verbose = TRUE) {
  species <- match.arg(species)
  annotation_engine <- match.arg(annotation_engine)
  inputs <- coerce_expranno_inputs(
    expr = expr,
    meta = meta,
    assay_name = assay_name,
    gene_id_col = gene_id_col,
    sample_col = sample_col
  )
  expr <- inputs$expr
  meta <- inputs$meta

  validate_expr(expr)
  validate_meta(meta)

  fields <- fields %||% default_annotation_fields()
  gene_info <- normalize_gene_ids(expr$gene_id, strip_version = strip_version)
  resolved_species <- resolve_species(species, gene_info$gene_id)
  meta_checked <- match_meta_to_expr(expr, meta)

  built <- build_annotation(
    gene_info = gene_info,
    species = resolved_species,
    annotation_engine = annotation_engine,
    fields = fields,
    verbose = verbose,
    biomart_version = biomart_version,
    biomart_host = biomart_host,
    biomart_mirror = biomart_mirror
  )

  annotation <- built$annotation
  provenance <- built$provenance
  expr_anno <- data.frame(
    annotation,
    expr[, -1, drop = FALSE],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  report <- annotation_coverage_report(annotation, fields)
  ambiguity_report <- annotation_ambiguity_report(annotation, fields)

  write_if_requested(expr_anno, output_file = output_file)
  write_if_requested(report, output_file = report_file)
  write_if_requested(ambiguity_report, output_file = ambiguity_file)
  write_if_requested(provenance, output_file = provenance_file)

  new_expranno_annotation(
    expr_anno = expr_anno,
    annotation = annotation,
    meta_checked = meta_checked,
    report = report,
    ambiguity_report = ambiguity_report,
    provenance = provenance,
    params = list(
      species = resolved_species,
      annotation_engine = annotation_engine,
      strip_version = strip_version,
      biomart_version = biomart_version,
      biomart_host = biomart_host,
      biomart_mirror = biomart_mirror,
      input_type = inputs$input_type
    )
  )
}

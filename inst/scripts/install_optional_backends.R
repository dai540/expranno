expranno_install_optional_backends <- function(
    profile = c("full", "core"),
    allow_failures = TRUE) {
  profile <- match.arg(profile)

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }

  options(repos = BiocManager::repositories())

  core_packages <- c(
    "AnnotationDbi",
    "AnnotationFilter",
    "annotate",
    "Biobase",
    "BiocParallel",
    "biomaRt",
    "ensembldb",
    "EnsDb.Hsapiens.v86",
    "EnsDb.Mmusculus.v79",
    "genefilter",
    "GSVA",
    "preprocessCore",
    "quantiseqr",
    "SummarizedExperiment",
    "S4Vectors",
    "SingleCellExperiment",
    "sva",
    "org.Hs.eg.db",
    "org.Mm.eg.db"
  )

  BiocManager::install(
    core_packages,
    ask = FALSE,
    update = FALSE
  )

  if (identical(profile, "core")) {
    return(invisible(TRUE))
  }

  install_optional_github <- function(repo) {
    if (isTRUE(allow_failures)) {
      tryCatch(
        remotes::install_github(repo, upgrade = "never"),
        error = function(e) {
          message(sprintf("Skipping optional GitHub backend %s: %s", repo, conditionMessage(e)))
        }
      )
    } else {
      remotes::install_github(repo, upgrade = "never")
    }
  }

  install_optional_github("dviraran/xCell")
  install_optional_github("GfellerLab/EPIC")
  install_optional_github("grst/MCPcounter")
  install_optional_github("grst/mMCPcounter")
  install_optional_github("cansysbio/ConsensusTME")

  if (!requireNamespace("immunedeconv", quietly = TRUE)) {
    if (isTRUE(allow_failures)) {
      tryCatch(
        remotes::install_github(
          "omnideconv/immunedeconv",
          upgrade = "never",
          dependencies = FALSE
        ),
        error = function(e) {
          message(sprintf("Skipping immunedeconv installation: %s", conditionMessage(e)))
        }
      )
    } else {
      remotes::install_github(
        "omnideconv/immunedeconv",
        upgrade = "never",
        dependencies = FALSE
      )
    }
  }

  invisible(TRUE)
}

if (identical(environmentName(environment()), "R_GlobalEnv") && !interactive()) {
  expranno_install_optional_backends(profile = "full", allow_failures = TRUE)
}

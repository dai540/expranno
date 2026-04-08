# expranno 2.5.2

- Reduced source-level file fragmentation by merging small helper modules into
  `R/data-helpers.R`, `R/bioconductor.R`, and `R/io.R`.
- Removed remaining non-package repository clutter so the tracked source tree
  stays closer to a minimal R package layout.
- Updated package metadata, README, and generated documentation pointers to
  match the simplified layout.

# expranno 2.5.1

- Reworked the source `README.md` and `_pkgdown.yml` so the package
  presentation matches the simpler `heteff`-style structure used across
  the repository set.

# expranno 2.5.0

- Added CI-level release metadata checks so `DESCRIPTION`, `NEWS.md`,
  `CITATION.cff`, `inst/CITATION`, and the README tarball example stay
  aligned.
- Added a required `optional-backend-core` GitHub Actions job that
  installs stable optional annotation and signature backends on `main`
  and runs the corresponding tests.
- Simplified the repository layout by removing smoke-only scripts,
  Docker helpers, pkgdown custom asset folders, and GitHub issue / PR
  template files that were not required for the package itself.

# expranno 2.3.1

- Fixed the pkgdown GitHub Actions workflow by allowing `pkgdown` to
  install the local package before building reference examples.

# expranno 2.3.0

- Synchronized the public version across `DESCRIPTION`, `NEWS.md`,
  `README.md`, citations, and the next pkgdown deployment.
- Slimmed the repository by dropping generated `docs/` from `main` and
  relying on the pkgdown workflow to rebuild and publish the site.
- Removed redundant figure assets and simplified the pkgdown workflow so
  the source tree stays smaller and easier to maintain.

# expranno 2.2.2

- Installed `SummarizedExperiment`, `S4Vectors`, and
  `SingleCellExperiment` in GitHub Actions so Bioconductor-output tests
  and vignettes pass consistently across the CI matrix.

# expranno 2.2.1

- Guarded the Bioconductor-output vignette example in `getting-started.Rmd`
  so `pkgdown`, `R CMD check`, and release builds no longer fail when
  `SummarizedExperiment` is not installed.

# expranno 2.2.0

- Added `example_annotation_truth()` plus bundled human and mouse truth
  tables so validation examples no longer depend on ad hoc local files.
- Added `as_expranno_se()` to convert annotated outputs back into a
  `SummarizedExperiment` with annotation in `rowData` and sample-level
  outputs in `colData`.
- Expanded `list_annotation_presets()` with recommended input scale,
  backend cascade, symbol priority, and bundled validation resource
  columns.
- Added a preset reference vignette and surfaced preset/truth/Bioconductor
  helpers more clearly in the README and pkgdown site.
- Added a bundled optional-backend installation helper and Docker recipe
  for more reproducible local backend validation.
- Strengthened public-package metadata and CI coverage for release
  workflows and cross-platform checks.

# expranno 2.1.0

- Reworked hybrid annotation around a fixed `biomaRt` release (`Ensembl v102`)
  with human/mouse-specific symbol recovery, backend provenance capture, and
  ambiguity reporting.
- Added explicit human and mouse annotation presets such as
  `human_tpm_v102` and `mouse_tpm_v102` for repeatable lab workflows.
- Added `annotation_report.csv`, `annotation_ambiguity.csv`,
  `annotation_provenance.csv`, and optional `session_info.txt` outputs to the
  main wrapper.
- Added `benchmark_annotation_engines()` so `none`, `biomart`, `orgdb`,
  `ensdb`, and `hybrid` can be compared on the same input.
- Added `validate_annotation_engines()` plus validation CSV outputs so
  annotation can be checked against truth tables rather than coverage alone.
- Added `SummarizedExperiment` / `SingleCellExperiment` input coercion through
  `as_expranno_input()`.
- Added a dedicated GitHub Actions workflow for optional-backend smoke tests.
- Improved signature scoring failures so empty annotation columns now raise a
  clear error instead of bubbling up as `no rows to aggregate`, and the
  optional smoke runner skips signature scoring when no symbols were annotated.

# expranno 2.0.2

- Fixed `EnsDb` fallback annotation by switching to
  `AnnotationFilter::GeneIdFilter()`.
- Fixed automatic deconvolution method discovery for current
  `immunedeconv` exports.
- Added skip-aware optional-backend smoke tests for hybrid annotation,
  real-data-like deconvolution, and end-to-end signature scoring.
- Added `inst/scripts/run_smoke_optional_backends.R` for local backend
  verification runs that write smoke-test outputs to disk.

# expranno 2.0.1

- Switched signature scoring to the current GSVA parameter-object API
  when available, with a legacy fallback for older GSVA releases.
- Added `expr_scale` and `duplicate_strategy` controls so duplicate
  symbols are no longer silently summed for all expression types.
- Added `deconv_args` support in the main wrapper and explicit handling
  of indication-specific `immunedeconv` methods such as `timer` and
  `consensus_tme`.
- Expanded tests and documentation around duplicate handling,
  deconvolution method requirements, and signature-scoring kernel
  choices.

# expranno 2.00

- Prepared GitHub-facing README, GitHub Actions, and pkgdown deployment
  for the public `dai540/expranno` repository.
- Redesigned the website workflow figure with a cleaner mermaid-based
  layout that separates inputs, core processing, outputs, and downstream
  analyses.
- Added package logo assets for the README and pkgdown site.
- Expanded the README and getting-started guide with clearer workflow
  framing, output interpretation, and post-run quality checks.
- Kept the human and mouse case studies centered on annotation coverage,
  merged outputs, and downstream deconvolution/signature workflows.

# expranno 1.00

- First public release candidate.
- Added hybrid human and mouse annotation workflow with `biomaRt`,
  `org.*.eg.db`, and optional `EnsDb` fallback.
- Added expression and metadata integration helpers.
- Added immune deconvolution wrappers around `immunedeconv`.
- Added GSVA and ssGSEA wrappers for user-supplied gene set databases.

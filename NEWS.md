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

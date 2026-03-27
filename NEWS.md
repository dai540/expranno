# expranno 2.00

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

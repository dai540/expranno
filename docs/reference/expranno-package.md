# expranno: annotation-first RNA-seq workflows

`expranno` annotates human and mouse expression matrices, aligns sample
metadata, optionally runs immune deconvolution with `immunedeconv`, and
computes GSVA or ssGSEA signature scores from user-supplied gene sets.

## Details

The package uses a coverage-first annotation strategy that can combine
`biomaRt`, `org.*.eg.db`, and optional `EnsDb` packages.

## Author

**Maintainer**: Daiki User <daiki@example.com> \[copyright holder\]

# expranno: annotation-first RNA-seq workflows

`expranno` annotates human and mouse expression matrices, aligns sample
metadata, optionally runs immune deconvolution with `immunedeconv`, and
computes GSVA or ssGSEA signature scores from user-supplied gene sets.

## Details

The package uses a coverage-first annotation strategy that can combine
`biomaRt`, `org.*.eg.db`, and optional `EnsDb` packages, and it returns
a compact annotation coverage report to support downstream quality
checks.

Built-in helpers such as
[`list_annotation_presets()`](https://dai540.github.io/expranno/reference/list_annotation_presets.md),
[`example_annotation_truth()`](https://dai540.github.io/expranno/reference/example_annotation_truth.md),
and
[`as_expranno_se()`](https://dai540.github.io/expranno/reference/as_expranno_se.md)
make it easier to standardize fixed human or mouse workflows, reproduce
validation runs, and move results back into Bioconductor containers.

## See also

Useful links:

- <https://github.com/dai540/expranno>

- <https://dai540.github.io/expranno/>

- Report bugs at <https://github.com/dai540/expranno/issues>

## Author

**Maintainer**: Dai Dai <daik54015@gmail.com>

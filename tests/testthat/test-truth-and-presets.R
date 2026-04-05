test_that("bundled truth tables expose the documented schema", {
  human_truth <- example_annotation_truth("human")
  mouse_truth <- example_annotation_truth("mouse")

  expect_identical(
    names(human_truth),
    c("gene_id", "symbol", "gene_name", "biotype")
  )
  expect_identical(names(mouse_truth), names(human_truth))
  expect_true(all(grepl("^ENSG", human_truth$gene_id)))
  expect_true(all(grepl("^ENSMUSG", mouse_truth$gene_id)))
})

test_that("preset table stays aligned with built-in validation helpers", {
  presets <- list_annotation_presets()

  expect_true(all(c(
    "annotation_preset", "recommended_input", "symbol_priority",
    "fallback_order", "bundled_truth"
  ) %in% names(presets)))
  expect_equal(nrow(presets), length(unique(presets$annotation_preset)))
  expect_true(all(grepl("^example_annotation_truth\\(", presets$bundled_truth)))
})

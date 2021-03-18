test_that("Loading through tabix should work", {
  results <-
    .load_tabix(here("inst/extdata/example_frag.bed.gz"), regions = "21")
  expect_true(
    nrow(results) == 1447909 &&
      nrow(results[is.na(mapq)]) &&
      sum(results$mapq, na.rm = TRUE) == 79990393 &&
      colnames(results) == c(
        "chrom",
        "start",
        "end",
        "name",
        "mapq",
        "strand",
        "cigar1",
        "cigar2"
      )
  )
})

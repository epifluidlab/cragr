test_that("Loading through tabix should work", {
  results <-
    load_fragments(here("inst/extdata/example_frag.bed.gz"), region = "21")
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


test_that("Loading whole fragment set should work", {
  results <-
    load_fragments(here("inst/extdata/example_frag.bed.gz"))
  expect_true(
    nrow(results) == 2969353 &&
      nrow(results[is.na(mapq)]) &&
      sum(results$mapq, na.rm = TRUE) == 164589277 &&
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


test_that("Filtering by chromosome should work", {
  dt <- load_fragments(here("inst/extdata/example_frag.small.bed.gz"))
  expect_true(nrow(.filter_by_chrom(dt, c("21", "22"))) == 94306)
  expect_true(nrow(.filter_by_chrom(dt, c("21", "22"), negate = TRUE)) == 0)
  expect_true(nrow(.filter_by_chrom(dt, c("1"))) == 0)
  expect_true(nrow(.filter_by_chrom(dt, c("21"))) == 52689)
  expect_true(nrow(.filter_by_chrom(dt, c("21"), negate = TRUE)) == 41617)
})


test_that("Filtering by genomic ranges should work", {
  granges <-
    tmp_gr <-
    GenomicRanges::GRanges(c("21:48119752-48119762", "22:16050481-16050485"))
  dt <-
    load_fragments(here("inst/extdata/example_frag.small.bed.gz"))
  dt_filtered <- .filter_by_granges(dt, granges, negate = FALSE)
  expect_equal(colnames(dt), colnames(dt_filtered))
  expect_equal(sapply(dt, typeof), sapply(dt_filtered, typeof))
  expect_equal(nrow(dt_filtered), 11L)
  expect_equal(sum(dt_filtered$mapq, na.rm = TRUE), 505L)

  dt_filtered <- .filter_by_granges(dt, granges, negate = TRUE)
  expect_equal(colnames(dt), colnames(dt_filtered))
  expect_equal(sapply(dt, typeof), sapply(dt_filtered, typeof))
  expect_equal(nrow(dt_filtered), 94295L)
  expect_equal(sum(dt_filtered$mapq, na.rm = TRUE), 4270833L)
})

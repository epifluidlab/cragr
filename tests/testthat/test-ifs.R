test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("Excluding regions works", {
  dt <-
    load_fragments(here("inst/extdata/example_frag.small.bed.gz"))
  region <-
    load_bed(here(
      "inst/extdata/wgEncodeDacMapabilityConsensusExcludable.hs37-1kg.bed"
    ))
  filtered <- .exclude_region(dt, region = NULL)
  expect_equal(nrow(filtered), 94306)
  filtered <- .exclude_region(dt, region = region, check_sort = TRUE)
  expect_equal(nrow(filtered), 81952)
})

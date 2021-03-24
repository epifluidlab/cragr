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


test_that("Constructing ovearll blacklist regions works", {
  black_list <- system.file("extdata", "wgEncodeDacMapabilityConsensusExcludable.hs37-1kg.bed", package = "cragr")
  maps <- system.file("extdata", "mappability.subset.bedGraph.gz", package = "cragr")

  regions <-
    suppressWarnings(
      .build_exclude_region(
        blacklist_region = black_list,
        mappability_region = maps,
        chrom = c("fake", "21")
      )
    )
  expect_equal(nrow(regions), 45568L)
  expect_equal(regions[, sum(end - start)], 6476167L)
  expect_equal(colnames(regions), c("chrom", "start", "end"))
  expect_equal(as.character(map_chr(regions, class)), c("factor", "integer", "integer"))

  regions <- suppressWarnings(.build_exclude_region(blacklist_region = NULL, mappability_region = maps, chrom = "19"))
  expect_null(regions)

  regions <- .build_exclude_region(blacklist_region = NULL, mappability_region = maps, chrom = "21")
  expect_equal(nrow(regions), 46592L)
  expect_equal(regions[, sum(end - start)], 6352800L)

  regions <- suppressWarnings(.build_exclude_region(blacklist_region = black_list, mappability_region = NULL, chrom = "21"))
  expect_equal(nrow(regions), 13L)
  expect_equal(regions[, sum(end - start)], 233249L)

  regions <- .build_exclude_region(blacklist_region = NULL, mappability_region = maps)
  expect_equal(nrow(regions), 104223L)
  expect_equal(regions[, sum(end - start)], 14455260L)

  expect_null(.build_exclude_region(blacklist_region = NULL, mappability_region = NULL))
})

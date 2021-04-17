#' Check if required binaries exist in PATH
#'
#' Some of cragr's functionalities require certain third-party binaries to exist
#' in PATH. For example, selectively loading cfDNA fragments from a certain
#' region under the hood invokes tabix to complete the task. This function checks
#' whether certain binaries exist.
#' @param binaries Character vector indicating what binaries to check. Default
#'   is tabix.
#' @param verbose Logical value indicating whether to output detailed messages.
#'   Default is FALSE.
#' @param stop_on_error Stop the script completely if the result is `FALSE`.
#' @return Logical value. TRUE if all required binaries can be found in PATH.
#' @export
check_binaries <- function(binaries = "tabix", verbose = FALSE, stop_on_fail = FALSE) {
  check_results <- binaries %>% purrr::map_lgl(function(binary) {
    # check if binary is in path
    if ("" == Sys.which(binary)) {
      if (verbose) {
        message(str_interp("Checking ${binary}: FAIL"))
      }
      return(FALSE)
    }
    else {
      if (verbose) {
        message(str_interp("Checking ${binary}: PASS"))
      }
      return(TRUE)
    }
  })

  check_results <- all(check_results)
  if (stop_on_fail)
    assertthat::assert_that(check_results, msg = paste0("Failed in locating ", paste(binaries, collapse = ", ")))
  else
    check_results
}



# To save time, pre-calculate low mappability regions
build_high_mappability_regions <- function(maps_file, chrom, genome, window_size, step_size, threshold) {
  stopifnot(length(chrom) == 1)
  assertthat::are_equal(window_size %% step_size, 0)
  chrom0 <- chrom

  chrom_sizes <- bedtorch::get_seqinfo(genome)
  chrom_sizes <-
    data.table(
      chrom = factor(seqnames(chrom_sizes), levels = as.character(seqnames(chrom_sizes))),
      start = 0L,
      end = GenomeInfoDb::seqlengths(chrom_sizes)
    ) %>%
    bedtorch::as.bedtorch_table(genome = genome)
  maps <- bedtorch::read_bed(maps_file, use_gr = FALSE, range = chrom, genome = genome)

  n <- window_size %/% step_size

  # Construct a continuous IFS score track, with filled 0s.

  scaffold <- chrom_sizes[chrom == chrom0]
  scaffold$start <- as.integer(maps$start[1] %/% step_size * step_size)
  scaffold$end <- as.integer(ceiling(scaffold$end / step_size) * step_size)

    # data.table::data.table(chrom = maps$chrom[1])[chrom_sizes, nomatch = 0, on = "chrom"]
  scaffold <- scaffold[, {
    start <- as.integer(seq(start, end, by = step_size) %>% head(-1))
    end <- as.integer(start + step_size)
    # gid <- 0:(size %/% step_size)
    # start <- as.integer(gid * step_size)
    # end <- as.integer((gid + 1) * step_size)
    list(chrom, start, end)
  }]

  data.table::setkey(scaffold, "chrom", "start", "end")

  # The original mappability scores may not be continuous, thus need scaffolding
  maps <- maps[scaffold]
  maps[, `:=`(
    end = as.integer(start + window_size),
    score = bedtorch::rollmean(score, k = n, na_pad = TRUE, align = "left")
  )]

  #
  maps[!is.na(score) & score >= threshold] %>%
    bedtorch::as.GenomicRanges()
}

# c(1:22, "X", "Y") %>% map(function(chrom) {
#   cat(str_interp("Processing chr${chrom}\n"))
#   build_high_mappability_regions("wgEncodeDukeMapabilityUniqueness35bp.hs37-1kg.20bp.bedGraph.gz",
#                                  chrom, "hs37-1kg", 200, 20, 0.9)
# })



log_mem <- function(label = "Unknown") {
  if (requireNamespace("lobstr")) {
    mem <- as.numeric(lobstr::mem_used()) / 1024**2
    logging::logdebug(str_interp("[${label}] Memory used: ${mem} MB"))
  }
}

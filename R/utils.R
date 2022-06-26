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
    } else {
      if (verbose) {
        message(str_interp("Checking ${binary}: PASS"))
      }
      return(TRUE)
    }
  })

  check_results <- all(check_results)
  if (stop_on_fail) {
    assertthat::assert_that(check_results, msg = paste0("Failed in locating ", paste(binaries, collapse = ", ")))
  } else {
    check_results
  }
}


get_style <- function(gr) {
  gr <- head(gr)
  new_seqinfo <-
    GenomeInfoDb::Seqinfo(seqnames = unique(as.character(seqnames(gr))))
  seqlevels(gr) <- seqlevels(new_seqinfo)
  seqinfo(gr) <- new_seqinfo

  GenomeInfoDb::seqlevelsStyle(gr)[1]
}


# Set style to either "NCBI" or "UCSC"
# The method in GenomeInfoDb may modify the search path :-()
set_style <- function(gr, style) {
  assert_that(style %in% c("NCBI", "UCSC"))

  path1 <- search()
  on.exit({
    path2 <- search()
    setdiff(path2, path1) %>% walk(function(x) {
      detach(
        name = x,
        force = TRUE,
        character.only = TRUE
      )
    })
  })
  original_style <- GenomeInfoDb::seqlevelsStyle(gr)[1]


  original_seqinfo <- GenomeInfoDb::seqinfo(gr)

  if (is_true(original_style == style)) {
    return(gr)
  }

  if (any(is.na(original_seqinfo@genome))) {
    GenomeInfoDb::seqlevelsStyle(gr) <- style
    return(gr)
  }

  original_genome <- original_seqinfo@genome[1]
  assert_that(original_genome %in% c("GRCh37", "GRCh38", "hg19", "hg38"))

  # Remove existing genomes and then change style
  new_seqinfo <-
    GenomeInfoDb::Seqinfo(seqnames = unique(as.character(seqnames(gr))))
  seqlevels(gr) <- seqlevels(new_seqinfo)
  seqinfo(gr) <- new_seqinfo
  GenomeInfoDb::seqlevelsStyle(gr) <- style

  # Set the final seqinfo
  if (is_true(original_genome == "GRCh37")) {
    new_genome <- "hg19"
  } else if (is_true(original_genome == "GRCh38")) {
    new_genome <- "hg38"
  } else if (is_true(original_genome == "hg19")) {
    new_genome <- "GRCh37"
  } else if (is_true(original_genome == "hg38")) {
    new_genome <- "GRCh38"
  } else {
    stop("Unknown genome: ", original_genome)
  }
  new_seqinfo <- bedtorch::get_seqinfo(new_genome)
  seqlevels(gr) <- seqlevels(new_seqinfo)
  seqinfo(gr) <- new_seqinfo

  return(gr)
}



# # To save time, pre-calculate low mappability regions
# build_high_mappability_regions <- function(maps_file, chrom, genome, window_size, step_size, threshold) {
#   stopifnot(length(chrom) == 1)
#   assertthat::are_equal(window_size %% step_size, 0)
#   chrom0 <- chrom

#   chrom_sizes <- bedtorch::get_seqinfo(genome)
#   chrom_sizes <-
#     data.table(
#       chrom = factor(seqnames(chrom_sizes), levels = as.character(seqnames(chrom_sizes))),
#       start = 0L,
#       end = GenomeInfoDb::seqlengths(chrom_sizes)
#     ) %>%
#     bedtorch::as.bedtorch_table(genome = genome)
#   maps <- bedtorch::read_bed(maps_file, use_gr = FALSE, range = chrom) #, genome = genome)

#   n <- window_size %/% step_size

#   # Construct a continuous IFS score track, with filled 0s.

#   scaffold <- chrom_sizes[chrom == chrom0]
#   scaffold$start <- as.integer(maps$start[1] %/% step_size * step_size)
#   scaffold$end <- as.integer(ceiling(scaffold$end / step_size) * step_size)

#   # data.table::data.table(chrom = maps$chrom[1])[chrom_sizes, nomatch = 0, on = "chrom"]
#   scaffold <- scaffold[, {
#     start <- as.integer(seq(start, end, by = step_size) %>% head(-1))
#     end <- as.integer(start + step_size)
#     # gid <- 0:(size %/% step_size)
#     # start <- as.integer(gid * step_size)
#     # end <- as.integer((gid + 1) * step_size)
#     list(chrom, start, end)
#   }]

#   data.table::setkey(scaffold, "chrom", "start", "end")

#   # The original mappability scores may not be continuous, thus need scaffolding
#   maps <- maps[scaffold]
#   maps[, `:=`(
#     end = as.integer(start + window_size),
#     score = bedtorch::rollmean(score, k = n, na_pad = TRUE, align = "left")
#   )]

#   #
#   maps[!is.na(score) & score >= threshold] %>%
#     bedtorch::as.GenomicRanges()
# }

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



#' Compare two hotspot files and return the Jaccard index
#'
#' Calculate Jaccard index using merged hotspots
#' @export
jaccard_index <- function(hotspot1, hotspot2, max_distance = 200L) {
  hotspot1 <- GenomicRanges::reduce(hotspot1, min.gapwidth = max_distance + 1)
  mcols(hotspot1) <- NULL
  hotspot2 <- GenomicRanges::reduce(hotspot2, min.gapwidth = max_distance + 1)
  mcols(hotspot2) <- NULL

  hits <- findOverlaps(hotspot1, hotspot2)
  hotspot1_common <- unique(queryHits(hits))
  hotspot1_unique <- setdiff(seq.int(length(hotspot1)), hotspot1_common)
  hotspot2_common <- unique(subjectHits(hits))
  hotspot2_unique <- setdiff(seq.int(length(hotspot2)), hotspot2_common)

  assertthat::are_equal(length(hotspot1_common), length(hotspot2_common))

  return(length(hotspot1_common) / (
    length(hotspot1_unique) + length(hotspot2_unique) + length(hotspot1_common)
  ))
}

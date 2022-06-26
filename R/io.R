#' Guess whether the seqinfo is hg19/hg38 or GRCh37/38 by guessing the style
#' @param genome Should be either GRCh37 or GRCh38
#' @export
assign_seqinfo <- function(gr, genome) {
  assert_that(!is_null(genome), genome %in% c("GRCh37", "GRCh38"))
  # Assign seqinfo
  genome_style <- get_style(gr)
  if (is_true(genome_style == "NCBI")) {
    if (is_true(genome == "GRCh37")) {
      frag_seqinfo <- bedtorch::get_seqinfo(genome = "GRCh37")
    } else {
      frag_seqinfo <- bedtorch::get_seqinfo(genome = "GRCh38")
    }
  } else if (is_true(genome_style == "UCSC")) {
    if (is_true(genome == "GRCh37")) {
      frag_seqinfo <- bedtorch::get_seqinfo(genome = "hg19")
    } else {
      frag_seqinfo <- bedtorch::get_seqinfo(genome = "hg38")
    }
  } else {
    stop(paste0("Unknown style: ", genome_style))
  }

  seqlevels(gr) <- seqlevels(frag_seqinfo)
  seqinfo(gr) <- frag_seqinfo
  gr
}

#' Read fragment BED file
#' @param genome Should be GRCh37 or GRCh38
#' @export
read_fragments <- function(file_path, range = NULL, genome = NULL, verbose = FALSE) {
  logging::logdebug("Reading BED data")
  frag <- bedtorch::read_bed(
    file_path = file_path,
    range = range
  )

  frag <- assign_seqinfo(frag, genome)

  GenomicRanges::strand(frag) <- "*"
  logging::logdebug("Done reading BED data")
  frag_metadata <- mcols(frag)

  if ("mapq" %in% colnames(frag_metadata)) {
    if (!is.integer(frag_metadata$mapq)) {
      logging::logerror("MAPQ should be integers")
      print(frag)
      stop()
    }
  } else if (ncol(frag_metadata) >= 2 && is.integer(frag_metadata[[2]])) {
    frag_metadata$mapq <- frag_metadata[[2]]
  } else {
    logging::logwarn("Cannot find MAPQ in the data frame")
  }

  # Only keep necessary columns
  meta_cols <- base::intersect(names(frag_metadata), c("mapq", "cigar1", "cigar2"))
  mcols(frag) <- frag_metadata[meta_cols]
  if (verbose) {
    print(frag)
  }
  frag
}

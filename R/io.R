#' @export
read_fragments <- function(file_path, range = NULL, genome = NULL) {
  logging::logdebug("Reading BED data")
  frag <- bedtorch::read_bed(file_path = file_path,
                             range = range,
                             genome = genome)
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
  print(frag)
  frag
}

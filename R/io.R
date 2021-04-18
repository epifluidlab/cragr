#' @export
read_fragments <- function(file_path, range = NULL, genome = NULL) {
  logging::logdebug("Reading BED data")
  frag <- bedtorch::read_bed(file_path = file_path,
                             range = range,
                             genome = genome)
  logging::logdebug("Done reading BED data")
  frag_metadata <- mcols(frag)

  mapq_guessed <- FALSE
  if (!"mapq" %in% colnames(frag_metadata)) {
    logging::logwarn("mapq not in column names. Need to guess which column is mapq")

    for (col_idx in seq_along(colnames(frag_metadata))) {
      if (is.integer(frag_metadata[[col_idx]])) {
        logging::logwarn(str_interp("Column ${colnames(frag_metadata)[col_idx]} seems to be mapq"))
        colnames(frag_metadata)[col_idx] <- "mapq"
        mapq_guessed <- TRUE
        break
      }
    }

    if (!"mapq" %in% colnames(frag_metadata))
      stop("Cannot find mapq in the data frame")
  }

  # Only need mapq. Drop other columns
  mcols(frag) <- frag_metadata["mapq"]
  if (mapq_guessed)
    print(frag)
  frag
}

#' Selectively load fragments using tabix
#'
#' @param file_path Path to the tsv file.
#' @param index_path Path to the .tbi index file.
#' @param regions Which genomic region to load. Can be either a character vector
#'   following the UCSC convention, e.g. 21:1-9411192, or a
#'   GenomicRanges::GRanges object.
#' @return A data.table object
.load_tabix <-
  function(file_path,
           index_path = paste0(file_path, ".tbi"),
           regions = NULL) {
    assertthat::assert_that(check_binaries(binaries = "tabix", verbose = FALSE))

    file_path <- normalizePath(file_path)
    index_path <- normalizePath(index_path)

    if (!file.exists(index_path)) {
      warning(str_interp("Failed to locat the index fine: ${index_path}"))
    }

    # Try to get column names from the header
    fields <-
      colnames(data.table::fread(file = file_path, nrows = 1000))
    # Remove the leading #-symbol from the 1st column name
    fields[1] <-  str_replace(fields[1], "^#[ ]*", "")

    na_strings <- c("NA", "na", "NaN", "nan", ".", "")
    if (is.null(regions) || str_trim(regions) == "") {
      results <-
        data.table::fread(file = file_path, na.strings = na_strings)
    } else if (class(regions) %in% c("GRanges", "character")) {
      tabix_regions <- paste0(as.character(regions), collapse = " ")
      tabix_cmd <- str_interp("tabix ${file_path} ${tabix_regions}")
      results <-
        data.table::fread(cmd = tabix_cmd, na.strings = na_strings)
    } else {
      stop(
        str_interp(
          "regions should be either a character vector or a GenomicRanges::GRanges object."
        )
      )
    }

    colnames(results) <- fields
    results[, chrom := as.character(chrom)][,]
  }

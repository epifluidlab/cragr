#' Load fragment data files or BED data files
#'
#' Fragment data files are essentially BED data files. `load_fragments` and
#' `load_bed` will load these data files into `data.table` objects.
#'
#' If `regions` is not `NULL`, `load_fragments` will invoke `tabix` to load the data.
#'
#' @param file_path Path to the data file (in tab-separated values format, a.k.a "tsv").
#' @param region Which genomic region to load. Can be either a character vector
#'   following the UCSC convention, e.g. `21:1-9411192`, or a
#'   GenomicRanges::GRanges object. Default is `NULL`.
#' @param chrom_sizes Path to the chromosome size file.
#' @return A `data.table` object
#' @export
load_bed <-
  function(file_path,
           region = NULL,
           chrom_sizes = system.file("extdata", "human_g1k_v37.chrom.sizes", package = "cragr")) {
    file_path <- normalizePath(file_path)

    # Try to get column names from the header
    fields <-
      colnames(data.table::fread(file = file_path, nrows = 1000))
    # Remove the leading #-symbol from the 1st column name
    fields[1] <- str_replace(fields[1], "^#[ ]*", "")

    filter_flag <- FALSE

    na_strings <- c("NA", "na", "NaN", "nan", ".", "")
    if (is.null(region) || str_trim(region) == "") {
      results <-
        data.table::fread(file = file_path, na.strings = na_strings)
    } else if (class(region) %in% c("GRanges", "character", "factor")) {
      # Check if we can use tabix to extract the reads
      tabix_exist <- check_binaries(binaries = "tabix", verbose = FALSE)
      index_exist <- file.exists(paste0(file_path, ".tbi"))
      use_tabix <- tabix_exist && index_exist

      if (use_tabix) {
        tabix_regions <- paste0(as.character(region), collapse = " ")
        tabix_cmd <- str_interp("tabix ${file_path} ${tabix_regions}")
        results <-
          data.table::fread(cmd = tabix_cmd, na.strings = na_strings)
      } else {
        if (endsWith(file_path, ".gz")) {
          warning("Either tabix or index does not exist. Resource to sequential scanning, which may negatively impact the performance.")
        }
        results <- data.table::fread(file = file_path, na.strings = na_strings)

        # Postpone the filtering to after assigning column names
        filter_flag <- TRUE
      }
    } else {
      stop(
        str_interp(
          "regions should be either a character vector or a GenomicRanges::GRanges object."
        )
      )
    }

    if (nrow(results) == 0) {
      NULL
    } else {
      # Rule: first 3 columns: chrom, start, end. 4th column: score if BEDGRAPH, name if BED
      if (all(str_detect(fields, "V[0-9]+"))) {
        # field names are default values, e.g. V1, V2, etc.
        data.table::setnames(results, old = 1:3, new = c("chrom", "start", "end"))
        if (typeof(results[[4]]) %in% c("integer", "double"))
          # BEDGRAPH
          data.table::setnames(results, old = 4, new = "score")
        else
          # BED
          data.table::setnames(results, old = 4, new = "name")
      } else {
        data.table::setnames(results, new = fields)
      }

      chrom_list <- data.table::fread(chrom_sizes)[[1]]
      results[, chrom := factor(chrom, levels = chrom_list)]
      # 4th column should be either character or numeric
      col_types <- sapply(results, typeof)
      if (!(col_types[4] %in% c("integer", "double")))
        data.table::set(results, j = 4L, value = as.character(results[[4L]]))

      if (filter_flag)
        results <- .filter_by_regions(results, region = region)

      data.table::setkeyv(results, c("chrom", "start", "end"))
      results[]
    }
  }


#' @rdname load_bed
#' @export
load_fragments <-
  function(file_path,
           region = NULL) {
    frags <- load_bed(
      file_path = file_path,
      region = region
    )
    fields <- colnames(frags)
    assertthat::assert_that(length(fields) >=6)
    fields[1:6] <- c("chrom", "start", "end", "name", "mapq", "strand")
    if (length(fields) >= 8)
      fields[7:8] <- c("cigar1", "cigar2")
    data.table::setnames(frags, fields)
    frags
  }


#' Filter the dataset by chromosomes
#'
#' @param dt A `data.table` object.
#' @param chrom Character vector containing chromosome names.
#' @param negate Indicate whether to include (`FALSE`) or exclude (`TRUE`) these
#'   fragments.
#' @return A filtered `data.table` object.
.filter_by_chrom <- function(dt, chrom, negate = FALSE) {
  if (nrow(dt) == 0)
    return(dt)

  if (length(chrom) == 0) {
    if (negate)
      return(dt)
    else {
      # Return an empty data.table
      return(dt[1][start == start + 1])
    }
  }

  chrom_list <- chrom
  if (negate)
    dt[!(chrom %in% chrom_list)]
  else
    dt[(chrom %in% chrom_list)]
}


#' Filter the dataset by genomic ranges
#'
#' @param dt A `data.table` object.
#' @param granges A `GenomicRanges::GRanges` object, or a compatible chracter vector.
#' @param negate Indicate whether to include (`FALSE`) or exclude (`TRUE`) these
#'   fragments.
#' @return A filtered `data.table` object.
.filter_by_granges <- function(dt, granges, negate = FALSE) {
  if (nrow(dt) == 0)
    return(dt)

  if (length(granges) == 0) {
    if (negate)
      return(dt)
    else {
      # Return an empty data.table
      return(dt[1][start == start + 1])
    }
  }

  assertthat::assert_that(check_binaries(binaries = "bedtools"))

  if (typeof(granges) == "character")
    granges <- GenomicRanges::GRanges(granges)

  # Convert granges to a data frame so that we can apply bedtools intersect on it
  granges <-
    granges %>%
    as_tibble() %>%
    rename(chrom = seqnames, start = start, end = end) %>%
    select(chrom:end)

  # After bedtools, the data types of some columns will change to character.
  # We need to convert them back to their original types.
  col_types <- map_chr(dt, typeof)

  params <- ifelse(negate, "-v -sorted", "-u -sorted")

  dt <- bedr::bedr(
    method = "intersect",
    params = params,
    input = list(a = dt, b = granges),
    check.chr = FALSE,
    check.zero.based = FALSE,
    check.valid = FALSE,
    check.sort = FALSE,
    check.merge = FALSE,
    verbose = FALSE,
    capture.output = "disk"
  )
  data.table::setDT(dt)

  # Convert back to numeric column types
  col_types2 <- map_chr(dt, typeof)
  fields <- colnames(dt)
  for (idx in which((col_types == "integer") & (col_types != col_types2))) {
    data.table::set(dt, j = idx, value = suppressWarnings(as.integer(dt[[idx]])))
  }
  for (idx in which((col_types == "double") & (col_types != col_types2))) {
    data.table::set(dt, j = idx, value = suppressWarnings(as.numeric(dt[[idx]])))
  }
  for (idx in which((col_types == "character") &
                    (col_types != col_types2))) {
    data.table::set(dt, j = idx, value = suppressWarnings(as.character(dt[[idx]])))
  }
  for (idx in which((col_types == "logical") & (col_types != col_types2))) {
    data.table::set(dt, j = idx, value = suppressWarnings(as.logical(dt[[idx]])))
  }
  dt
}


#' Filter the dataset by regions
#'
#' @param dt A `data.table` object.
#' @param region Can be either `GenomicRanges::GRanges` or character vector
#' @param negate Indicate whether to include (`FALSE`) or exclude (`TRUE`) these
#'   fragments.
#' @return A filtered `data.table` object.
.filter_by_regions <- function(dt, region, negate = FALSE) {
  # Separate region to chromosome names and genomic ranges
  if (is.character(region)) {
    granges <-
      GenomicRanges::GRanges(region[str_detect(region, "^[^:]+:[0-9]+-[0-9]+$")])
    chrom <- region[str_detect(region, "^[^:]+$")]
  } else {
    granges <- region
    chrom <- NULL
  }

  if (negate) {
    dt <- .filter_by_chrom(dt, chrom = chrom, negate = negate)
    dt <- .filter_by_granges(dt, granges = granges, negate = negate)
  } else {
    # In this case, .filter_by_chrom and .filter_by_granges may have common
    # returns. Need to deduplicate them.
    dt[, id := seq_along(chrom)]
    dt1 <- .filter_by_chrom(dt, chrom = chrom)
    dt2 <- .filter_by_granges(dt, granges = granges)
    dt[, id := NULL]
    dt <- unique(data.table::rbindlist(list(dt1, dt2)), by = "id")
    dt[, id := NULL]
  }
  dt
}


#' Write data frame to disk in BED format
#'
#' If the target file indicates gzip-compression, will invoke bgzip rather than
#' standard gzip. Another feature is, if the output file should contain header
#' lines, they will start with `#`, as specified in BED format specification.
#' @param dt Data frame to write.
#' @param file_path Output file name.
#' @param col_names A logical value indicating whether to output colum names.
#' @param create_index A logical value indicating whether to create an index
#'   file. Ignored if the output is not gzipped.
#' @export
write_bed <- function(dt, file_path, col_names = TRUE, create_index = FALSE, ...) {
  is_gzipped <- endsWith(file_path, ".gz")
  scipen <- 999

  if (col_names == TRUE) {
    fields <- colnames(dt)
    old_fields <- data.table::copy(fields)
    data.table::setnames(dt, old = fields[1], new = paste0("#", fields[1]))
    on.exit(data.table::setnames(dt, old_fields), add = TRUE)
  }

  if (is_gzipped) {
    # bgzip should be present
    check_binaries(binaries = "bgzip", stop_on_fail = TRUE)
    if (create_index)
      check_binaries(binaries = "tabix", stop_on_fail = TRUE)

    temp_bed <- tempfile(fileext = ".bed")
    on.exit(file.remove(temp_bed), add = TRUE)

    data.table::fwrite(
      dt,
      file = temp_bed,
      col.names = col_names,
      sep = "\t",
      quote = FALSE,
      scipen = scipen,
      na = ".",
      ...
    )

    system(str_interp("bgzip < ${temp_bed} > ${file_path}"))
    if (create_index)
      system(str_interp("tabix -p bed ${file_path}"))
  } else {
    data.table::fwrite(
      dt,
      file = file_path,
      col.names = col_names,
      sep = "\t",
      quote = FALSE,
      scipen = scipen,
      na = ".",
      ...
    )
  }
}


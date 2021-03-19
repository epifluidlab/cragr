#' Calculate IFS scores
#'
#' @param fragment_data A `data.table` containing cfDNA fragment data.
#' @param chrom Calculate IFS scores for certain chromosomes. If `NULL`, all
#'   chromosomes in `fragment_data` will be used.
#' @param window_size Width of the sliding window.
#' @param excluded_regions Those fragments whose midpoint fall within any of the
#'   excluded regions will not be used in the analysis. `excluded_regions` can
#'   be either a character vector that contains the names of the files defining
#'   the regions (should be in BED format), or a list of `data.table` objects.
#' @param chrom_sizes Path to the chromosome size file.
#' @param min_mapq Exclude fragments whose `MAPQ` scores are smaller than `min_mapq`.
#' @param min_fraglen Exclude fragments shorter than `min_fraglen`.
#' @param max_fraglen Exclude fragments longer than `max_fraglen`.
#' @param exclude_soft_clipping For fragments with leading soft clipping, their
#'   inferred lengths are less reliable. If `TRUE`, such fragments will be
#'   excluded from the analysis.
#' @export
ifs_score <-
  function(fragment_data,
           chrom = NULL,
           window_size = 20L,
           excluded_regions = NULL,
           chrom_sizes = system.file("extdata", "human_g1k_v37.chrom.sizes", package = "cragr"),
           min_mapq = 30L,
           min_fraglen = 50L,
           max_fraglen = 1000L,
           exclude_soft_clipping = FALSE) {
    if (is.null(chrom)) {
      chrom <- unique(fragment_data$chrom)
    }
    if (length(chrom) > 1) {
      # For multiple chromosomes, sequentially
      ifs <-
        data.table::rbindlist(lapply(chrom, function(chrom) {
          ifs_score(
            fragment_data,
            chrom = chrom,
            excluded_regions = excluded_regions,
            min_mapq = min_mapq,
            exclude_soft_clipping = exclude_soft_clipping
          )
        }))
      data.table::setkeyv(ifs, c("chrom", "start", "end"))
      return(ifs)
    }

    chrom0 <- chrom
    fragment_data <- fragment_data[chrom == chrom0]
    rm(chrom0)

    # Apply filters
    if (min_mapq > 0)
      fragment_data <- fragment_data[mapq >= min_mapq]

    fragment_data[, length := end - start]
    min_fraglen <- min_fraglen %||% 0
    max_fraglen <- max_fraglen %||% Inf
    fragment_data <- fragment_data[between(length, min_fraglen, max_fraglen)]

    if (exclude_soft_clipping && all(c("cigar1", "cigar2") %in% colnames(fragment_data)))
      fragment_data <- fragment_data[!str_detect(cigar1, "^[0-9]+S") & !str_detect(cigar2, "[0-9]+S$")]

    # Each fragment is characterized by the midpoint rather than the start-end pair
    #
    # After this operation:
    #        chrom    start      end mapq length   frag_s
    # 1:        22 16050320 16050321   48    145 16050248
    # 2:        22 16050324 16050325   48    147 16050251
    # 3:        22 16050320 16050321   48    135 16050252
    # 4:        22 16050453 16050454   49    118 16050394
    # 5:        22 16050552 16050553   48    143 16050481
    # ---
    # 26122:    22 51242582 51242583   46     91 51242537
    # 26123:    22 51243136 51243137   46    112 51243080
    # 26124:    22 51243190 51243191   46    116 51243132
    # 26125:    22 51243212 51243213   45    133 51243146
    # 26126:    22 51243198 51243199   47    103 51243147
    fragment_data <-
      fragment_data[
        , .(chrom, start, end, mapq, length)
      ][, `:=`(
        midpoint = round((start + end) / 2),
        frag_s = start)][
          , `:=`(start = midpoint,
                 end = midpoint + 1,
                 midpoint = NULL)
        ]
    data.table::setkeyv(fragment_data, c("chrom", "start", "end"))

    if (!is.null(excluded_regions)) {
      fragment_data <-
        .exclude_region(
          fragment_data,
          region = excluded_regions,
          check_sort = TRUE,
          chrom_sizes = chrom_sizes
        )
    }

    fragment_data[, window_id := start %/% window_size]
    avg_len <- mean(fragment_data$length)
    ifs <-
      fragment_data[,
                    .(chrom = chrom[1],
                      cov = .N,
                      len_sum = sum(length, na.rm = TRUE)),
                    by = .(window_id)]
    ifs[, `:=`(
      score = cov + len_sum / avg_len,
      start = window_id * window_size,
      end = (window_id + 1) * window_size
    )]
    ifs <- ifs[, .(chrom, start, end, score, cov)]
    data.table::setkeyv(ifs, c("chrom", "start", "end"))
    ifs
  }


#' Exclude certain regions from the dataset, using `bedtools`
#'
#' @param dt A `data.table` object containing fragment data.
#' @param region A `data.table` defining excluded regions.
#' @param check_sort `region` may not be correctly sorted. If `TRUE`, will sort
#'   `region` by both chromosomes and positions, based on `chrom_sizes`
#' @param chrom_sizes Path to the chromosome size file.
#' @return A filtered `data.table` object.
.exclude_region <-
  function(dt,
           region,
           check_sort = FALSE,
           chrom_sizes = system.file("extdata", "human_g1k_v37.chrom.sizes", package = "cragr")) {
    if (nrow(dt) == 0 || is.null(region) || nrow(region) == 0)
      return(dt)

    region <- region[, .(chrom, start, end)]

    if (check_sort) {
      # Make sure region is correctly sorted
      chrom_list <- data.table::fread(chrom_sizes)[[1]]
      region[, chrom := factor(chrom, levels = chrom_list)]
      region <- region[order(chrom)]
    }

    # After bedtools, the data types of some columns will change to character.
    # We need to convert them back to their original types.
    col_types <- map_chr(dt, typeof)

    dt <- bedr::bedr(
      method = "intersect",
      params = paste0("-v -sorted -g ", chrom_sizes),
      input = list(a = dt, b = region[,1:3]),
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
    for (idx in which((col_types == "integer") &
                      (col_types != col_types2))) {
      data.table::set(dt, j = idx, value = suppressWarnings(as.integer(dt[[idx]])))
    }
    for (idx in which((col_types == "double") &
                      (col_types != col_types2))) {
      data.table::set(dt, j = idx, value = suppressWarnings(as.numeric(dt[[idx]])))
    }
    for (idx in which((col_types == "character") &
                      (col_types != col_types2))) {
      data.table::set(dt, j = idx, value = suppressWarnings(as.character(dt[[idx]])))
    }
    for (idx in which((col_types == "logical") &
                      (col_types != col_types2))) {
      data.table::set(dt, j = idx, value = suppressWarnings(as.logical(dt[[idx]])))
    }

    chrom_list <- data.table::fread(chrom_sizes)[[1]]
    dt[, chrom := factor(as.character(chrom), levels = chrom_list)]
    data.table::setkeyv(dt, c("chrom", "start", "end"))
    dt
  }

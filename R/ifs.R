#' Calculate IFS scores
#'
#' @param fragment_data A `data.table` containing cfDNA fragment data.
#' @param chrom Calculate IFS scores for certain chromosomes. If `NULL`, all
#'   chromosomes in `fragment_data` will be used.
#' @param window_size Width of the sliding window. Must be multiples of `step_size`.
#' @param step_size Incremental steps of the sliding window.
#' @param gc A `data.table` or character vector indicating the G+C% track. If
#'   `NULL`, do not perform GC-correction.
#' @param blacklist_region Those fragments whose midpoint fall within any of the
#'   excluded regions will not be used in the analysis. `blacklist_region` can
#'   be either a character vector that contains the names of the files defining
#'   the regions (should be in BED format), a `data.table` object, or a list of
#'   `data.table` objects.
#' @param mappability_region A `data.table` object of mappability scores, or a
#'   character vector of the mappability score file path.
#' @param mappability_threshold The mappability score threshold, below which the
#'   region is considered as low-mappability and excluded from analysis.
#' @param chrom_sizes Path to the chromosome size file.
#' @param min_mapq Exclude fragments whose `MAPQ` scores are smaller than `min_mapq`.
#' @param min_fraglen Exclude fragments shorter than `min_fraglen`.
#' @param max_fraglen Exclude fragments longer than `max_fraglen`.
#' @param exclude_soft_clipping For fragments with leading soft clipping, their
#'   inferred lengths are less reliable. If `TRUE`, such fragments will be
#'   excluded from the analysis. `cigar1` and `cigar2` should be present in the
#'   dataset.
#' @export
ifs_score <-
  function(fragment_data,
           chrom = NULL,
           window_size = 200L,
           step_size = 20L,
           gc = NULL,
           blacklist_region = NULL,
           mappability_region = NULL,
           mappability_threshold = 0.9,
           chrom_sizes = system.file("extdata", "human_g1k_v37.chrom.sizes", package = "cragr"),
           min_mapq = 30L,
           min_fraglen = 50L,
           max_fraglen = 1000L,
           exclude_soft_clipping = FALSE) {
    assertthat::are_equal(window_size %% step_size, 0)
    if (is.null(chrom)) {
      chrom <- as.character(unique(fragment_data$chrom))
    } else {
      chrom <- as.character(chrom)
    }

    if (length(chrom) > 1) {
      # For multiple chromosomes, sequentially
      ifs <-
        data.table::rbindlist(lapply(chrom, function(chrom) {
          ifs_score(
            fragment_data,
            chrom = chrom,
            window_size = window_size,
            step_size = step_size,
            gc = gc,
            blacklist_region = blacklist_region,
            mappability_region = mappability_region,
            mappability_threshold = mappability_threshold,
            chrom_sizes = chrom_sizes,
            min_mapq = min_mapq,
            min_fraglen = min_fraglen,
            max_fraglen = max_fraglen,
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

    if (exclude_soft_clipping)
      fragment_data <- fragment_data[!str_detect(cigar1, "^[0-9]+S") & !str_detect(cigar2, "[0-9]+S$")]

    # Exclude fragments from certain regions
    excluded_regions <-
      .build_exclude_region(
        chrom = chrom,
        blacklist_region = blacklist_region
        # mappability_region = mappability_region,
        # mappability_threshold = mappability_threshold
      )

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


    # Exclude dark regions
    if (!is.null(excluded_regions)) {
      fragment_data <-
        .exclude_region(
          fragment_data,
          region = excluded_regions,
          # check_sort = TRUE,
          chrom_sizes = chrom_sizes
        )
    }


    # Get coverage and calcuate IFS score for each 20bp (step-size) window
    fragment_data[, window_id := start %/% step_size]
    avg_len <- mean(fragment_data$length)
    ifs <-
      fragment_data[,
                    .(chrom = chrom[1],
                      cov = .N,
                      len_sum = sum(length, na.rm = TRUE)),
                    by = .(window_id)]
    ifs[, `:=`(
      score = cov + len_sum / avg_len,
      start = window_id * step_size,
      end = (window_id + 1) * step_size
    )]
    ifs <- ifs[, .(chrom, start, end, score, cov)]
    data.table::setkeyv(ifs, c("chrom", "start", "end"))

    # Perform rolling sum over the sliding windows, therefore we have results
    # for rolling windows (200bp, window_size) at step size of (20bp, step_size)
    ifs <- .slide_window(ifs, window_size = window_size, step_size = step_size)

    # Exclude low-mappability regions
    if (!is_null(mappability_region)) {
      if (is_character(mappability_region)) {
        mappability_region <- load_bed(mappability_region, region = chrom)
      } else {
        chrom0 <- chrom
        mappability_region <- mappability_region[chrom == chrom0]
        rm(chrom0)
      }

      ifs <- .exclude_low_mappability(ifs, mappability = mappability_region, mappability_threshold = mappability_threshold)
      rm(mappability_region)
    }

    # Don't perform GC correctoion
    if (!is_null(gc)) {
      if (is_character(gc)) {
        gc <- load_bed(gc, region = chrom)
      } else {
        chrom0 <- chrom
        gc <- gc[chrom == chrom0]
        rm(chrom0)
      }

      ifs <- .gc_correct(ifs, gc, span = 0.75)
      rm(gc)
    }

    ifs[]
  }


#' Perform GC-correction for IFS scores
#'
#' @param ifs A `data.table` of original IFS scores.
#'   in-place.
#' @param gc A `data.table` of GC contents.
#' @param span `span` used in LOESS fit.
#' @return A `data.table` containing GC-corrected IFS scores/
.gc_correct <- function(ifs, gc, span = 0.75) {
  chrom <- as.character(unique(ifs$chrom))
  assertthat::are_equal(length(chrom), 1)
  assertthat::are_equal(chrom, as.character(unique(gc$chrom)))

  check_binaries(stop_on_fail = TRUE)
  # Map GC to fragment BEDs
  gc <- bedr::bedr(method = "map",
                   input = list(a = ifs[, .(chrom, start, end)],
                                b = gc[, .(chrom, start, end, score)]),
                   params = "-c 4 -o mean",
                   check.chr = FALSE,
                   check.zero.based = FALSE,
                   check.valid = FALSE,
                   check.sort = FALSE,
                   check.merge = FALSE,
                   verbose = FALSE,
                   capture.output = "disk")

  data.table::setDT(gc)
  data.table::setnames(gc, c("chrom", "start", "end", "score"))
  gc[, chrom := factor(chrom, levels = levels(ifs$chrom))]
  data.table::setkey(gc, "chrom", "start", "end")

  ifs <- ifs[gc, nomatch = 0][, `:=`(gc = i.score, i.score = NULL)]
  ifs[, score0 := score]
  ifs[score0 > 0, score := score0 - lowess(x = gc, y = score0, f = span)$y + mean(score0)]

  # Also perform GC-correction for coverage
  ifs[score > 0, cov_corrected := cov - lowess(x = gc, y = cov, f = span)$y + mean(cov)]
  ifs[is.na(cov_corrected), cov_corrected := 0]

  # ifs[, .(chrom, start, end, score, cov, gc, score0)]
  ifs[]
}


.exclude_low_mappability <- function(ifs, mappability, mappability_threshold) {
  chrom <- as.character(unique(ifs$chrom))
  assertthat::are_equal(length(chrom), 1)
  assertthat::are_equal(chrom, as.character(unique(mappability$chrom)))

  check_binaries(stop_on_fail = TRUE)
  # Map mappability to fragment BEDs
  mappability <- bedr::bedr(method = "map",
                   input = list(a = ifs[, .(chrom, start, end)],
                                b = mappability[, .(chrom, start, end, score)]),
                   params = "-c 4 -o mean",
                   check.chr = FALSE,
                   check.zero.based = FALSE,
                   check.valid = FALSE,
                   check.sort = FALSE,
                   check.merge = FALSE,
                   verbose = FALSE,
                   capture.output = "disk")

  data.table::setDT(mappability)
  data.table::setnames(mappability, c("chrom", "start", "end", "score"))
  mappability[, chrom := factor(chrom, levels = levels(ifs$chrom))]
  data.table::setkey(mappability, "chrom", "start", "end")

  ifs <- ifs[mappability, nomatch = 0][, `:=`(mappability = i.score, i.score = NULL)]
  ifs <- ifs[mappability >= mappability_threshold]
  # ifs[, .(chrom, start, end, score, cov, gc, mappability, score0)]
  ifs[]
}


#' Perform a sliding window over raw IFS scores
#'
#' Raw IFS scores are calculated for each `step_size` bin. This function
#' performs a sliding window with the width of `window_size`. Each window
#' contains a whole number of `step_size` bins.
#' @param window_size Width of the sliding window. Must be multiples of `step_size`.
#' @param step_size Incremental steps of the sliding window.
.slide_window <- function(ifs, window_size, step_size) {
  assertthat::are_equal(window_size %% step_size, 0)

  n <- window_size %/% step_size

  # Construct a continuous IFS score track, with filled 0s.
  scaffold <-
    data.table::data.table(chrom = ifs$chrom[1],
                           start = seq(min(ifs$start), max(ifs$start),
                                       by = step_size))
  scaffold[, end := start + step_size]
  data.table::setkey(scaffold, "chrom", "start", "end")
  ifs <- ifs[scaffold][is.na(score), score := 0][is.na(cov), cov := 0]

  ifs <- ifs[, `:=`(
    end = start + window_size,
    score = zoo::rollsum(
      score,
      k = n,
      fill = NA,
      align = "left"
    ),
    cov = zoo::rollsum(
      cov,
      k = n,
      fill = NA,
      align = "left"
    )
  )][!is.na(score) & !is.na(cov)]

  # Eliminate float-point errors
  ifs[abs(score) < 1e-9, score := 0]

  data.table::setkey(ifs, "chrom", "start", "end")
  ifs[score > 0]



  # # Build a matrix: N by n, where N is the number of IFS rows. If window_size is
  # # 200 and step_size is 20, n is 10, means that at most we need combine 10
  # # step_size bins into a window bin.
  # # Each column is the IFS scores, but incrementally shifted.
  # m_ifs <- 0:(n-1) %>% map(function(s) {
  #   if (s == 0) return(ifs$score)
  #   data.table::shift(ifs$score, n = -s)
  # }) %>% unlist() %>% matrix(ncol = n)
  #
  # m_cov <- 0:(n-1) %>% map(function(s) {
  #   if (s == 0) return(ifs$cov)
  #   data.table::shift(ifs$cov, n = -s)
  # }) %>% unlist() %>% matrix(ncol = n)
  #
  # # A matrix of the same dimension. Each column is the shifted start position.
  # m_pos <- 0:(n-1) %>% map(function(s) {
  #   if (s == 0) return(ifs$start)
  #   data.table::shift(ifs$start, n = -s)
  # }) %>% unlist() %>% matrix(ncol = n)
  #
  # # Indicate whether the element should be included in the window. If an element
  # # is too far way (further than the window_size), it should be excluded.
  # m_include <- (m_pos <= m_pos[,1] + window_size - step_size)
  #
  # ifs <-
  #   ifs[, `:=`(
  #     end = start + window_size,
  #     cov = rowSums(m_cov * m_include),
  #     score = rowSums(m_ifs * m_include)
  #   )][!is.na(score)]
  # data.table::setkey(ifs, "chrom", "start", "end")
  # ifs[]
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

    check_binaries(stop_on_fail = TRUE)
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



#' Construct the overall excluded region
#'
#' @param chrom A character vector. If not `NULL`, only regions in these
#'   chromosomes will be considered.
#' @param blacklist_region Can be a `data.table` object, a list of `data.table`
#'   object, or a character vector of the data file paths.
#' @param mappability_region A `data.table` object of mappability scores, or a
#'   character vector of the mappability score file path.
#' @param mappability_threshold The mappability score threshold, below which the
#'   region is considered as low-mappability and excluded from analysis.
.build_exclude_region <-
  function(chrom = NULL,
           blacklist_region = NULL) {
           # mappability_region = NULL,
           # mappability_threshold = 0.9) {
    # if (is_null(blacklist_region) && is_null(mappability_region)) return(NULL)
    if (is_null(blacklist_region)) return(NULL)

    if (is_list(blacklist_region)) {
      blacklist_region %<>% data.table::rbindlist()
    } else if (is_character(blacklist_region)) {
      blacklist_region %<>% map(~ load_bed(., region = chrom)) %>% data.table::rbindlist()
    }

    # if (is_character(mappability_region)) {
    #   mappability_region %<>% load_bed(region = chrom)
    # }
    # if (!is_null(mappability_region)) {
    #   mappability_region <-
    #     mappability_region[score < mappability_threshold]
    # }
    #
    # merged_region <-
    #   data.table::rbindlist(list(blacklist_region, mappability_region),
    #                         fill = TRUE)
    # if (is_null(merged_region) || nrow(merged_region) == 0) return(NULL)

    merged_region <- blacklist_region
    merged_region <- merged_region[, .(chrom, start, end)]
    data.table::setkey(merged_region, "chrom", "start", "end")

    chrom_list <- levels(merged_region$chrom)
    merged_region <- bedr::bedr(
      method = "merge",
      input = list(i = merged_region),
      check.chr = FALSE,
      check.zero.based = FALSE,
      check.merge = FALSE,
      check.valid = FALSE,
      check.sort = FALSE,
      verbose = FALSE
    )
    data.table::setDT(merged_region)
    merged_region[, chrom := factor(chrom, levels = chrom_list)]
    data.table::setkey(merged_region, "chrom", "start", "end")
    merged_region[]
  }


#' Build a continuous counterpart of Poisson distribution
.build_cpois <- function(lambda) {
  # Determine the normalization factor
  if (lambda > 10)
    factor <- 1
  else {
    f <- function(x) {
      exp(-lambda) * lambda ^ x / gamma(x + 1)
    }
    factor <- 1 / integrate(f, lower = 0, upper = 100)$value
  }

  # Continuous-Poisson density
  cpois_density <- function(x) {
    factor * exp(-lambda) * lambda ^ x / gamma(x + 1)
  }

  function(x) {
    sapply(X = x,
           FUN = function(v) integrate(cpois_density, lower = 0, upper = v)$value)
  }
}


#'
#'
#' @param cpois A logical value indicating whether to calculate p-values based
#'   on continuous Poisson model.
#' @param local_layout Example: list(`5k` = 5000L, `10k` = 10000L)
#' @export
call_peak <- function(ifs, cpois = FALSE) {
  # Standard Poisson model
  ifs[, pval := ppois(score, lambda = mean(score)), by = chrom][, pval_adjust := p.adjust(pval, method = "BH"), by = chrom]
  if (cpois) {
    # Fit the continuous Poisson model
    d <-
      ifs[, {
        cpois = .build_cpois(mean(score))
        list(pval = cpois(score))
      }, by = chrom]
    ifs[, pval_cpois := d$pval][, pval_cpois_adjust := p.adjust(pval_cpois, method = "BH"), by = chrom]
  }
  ifs
}


#'
#' @export
call_peak_local <- function(ifs, window_size, step_size, local_layout) {
  assertthat::are_equal(window_size %% step_size, 0)
  assertthat::assert_that(all(unlist(local_layout) %% step_size == 0))

  # Construct a continuous IFS score track, with filled 0s.
  scaffold <-
    ifs[, .(start = seq(min(start), max(start), by = step_size)), by = chrom][, end := start + window_size]
  data.table::setkey(scaffold, "chrom", "start", "end")
  ifs <- ifs[scaffold]

  for (local_suffix in names(local_layout)) {
    local_width <- local_layout[[local_suffix]]

    ifs[, score_rollmean := zoo::rollmean(score, k = local_width %/% step_size, fill = NA, align = "center", na.rm = TRUE), by = chrom]
    ifs[, (paste0("pval_", local_suffix)) := ppois(score, score_rollmean)][, score_rollmean := NULL]
  }
  ifs[!is.na(score)]
}


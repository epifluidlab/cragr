preprocess_fragment <-
  function(fragment_data,
           chrom,
           min_mapq,
           min_fraglen,
           max_fraglen,
           exclude_soft_clipping,
           blacklist_region) {
    stopifnot(length(chrom) == 1)

    fragment_data <-
      keepSeqlevels(fragment_data, chrom, pruning.mode = "coarse")

    # Apply filters
    logging::loginfo("Applying MAPQ filter ...")
    if (min_mapq > 0)
      fragment_data <- fragment_data[fragment_data$mapq >= min_mapq]
    logging::loginfo(str_interp("${length(fragment_data)} fragments have passed filter"))

    logging::loginfo("Applying fragment length filter ...")
    min_fraglen <- min_fraglen %||% 0
    max_fraglen <- max_fraglen %||% Inf
    fragment_data <-
      fragment_data[between(width(fragment_data), min_fraglen, max_fraglen)]
    logging::loginfo(str_interp("${length(fragment_data)} fragments have passed filter"))


    if (exclude_soft_clipping &&
        all(c("cigar1", "cigar2") %in% colnames(mcols(fragment_data))))
      fragment_data <-
      fragment_data[!str_detect(fragment_data$cigar1, "^[0-9]+S") &
                      !str_detect(fragment_data$cigar2, "[0-9]+S$")]

    # Exclude fragments from certain regions
    excluded_regions <-
      build_exclude_region(chrom = chrom,
                            blacklist_region = blacklist_region)

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
    fragment_data$frag_s <- start(fragment_data)
    fragment_data$length <- width(fragment_data)
    midpoint = round((start(fragment_data) + end(fragment_data)) / 2)
    start(fragment_data) <- midpoint
    width(fragment_data) <- 1

    logging::loginfo("Applying blacklist region filter ...")
    # Exclude dark regions
    logging::loginfo("Excluding blacklist regions ...")
    if (!is.null(excluded_regions)) {
      fragment_data <-
        fragment_data[-queryHits(findOverlaps(fragment_data, excluded_regions))]
    }
    logging::loginfo(str_interp("${length(fragment_data)} fragments have passed filter"))

    fragment_data
  }


calculate_raw_ifs <- function(fragment_data, window_size, step_size) {
  # Get coverage and calcuate IFS score for each 20bp (step-size) window
  logging::loginfo("Calculating raw IFS scores ...")

  # Here we need do grouping and aggregating.
  # It seems data.table is way faster
  fragment_data$window_id <- (start(fragment_data) - 1) %/% step_size

  fragment_data <- bedtorch::as.bedtorch_table(fragment_data)
  avg_len <- mean(fragment_data$length)
  ifs <- fragment_data[,
                       .(chrom = chrom[1],
                         start = window_id * step_size,
                         end = (window_id + 1) * step_size,
                         score = .N + sum(length, na.rm = TRUE) / avg_len,
                         cov = .N),
                       by = window_id][, window_id := NULL]
  ifs %<>%
    bedtorch::as.bedtorch_table(genome = attr(fragment_data, "genome"))
  data.table::setkeyv(ifs, c("chrom", "start", "end"))

  # Perform rolling sum over the sliding windows, therefore we have results
  # for rolling windows (200bp, window_size) at step size of (20bp, step_size)
  ifs <-
    slide_window(
      ifs,
      window_size = window_size,
      step_size = step_size
    )
}


#' Calculate IFS scores
#'
#' @param fragment_data A `data.table` containing cfDNA fragment data.
#' @param chrom Calculate IFS scores for certain chromosomes. If `NULL`, all
#'   chromosomes in `fragment_data` will be used.
#' @param window_size Width of the sliding window. Must be multiples of `step_size`.
#' @param step_size Incremental steps of the sliding window.
#' @param gc A `data.table` or character vector indicating the G+C% track. If
#'   `NULL`, do not perform GC-correction.
#' @param gc_correct Logical value. Whether to perform GC correction.
#' @param blacklist_region Those fragments whose midpoint fall within any of the
#'   excluded regions will not be used in the analysis. `blacklist_region` can
#'   be either a character vector that contains the names of the files defining
#'   the regions (should be in BED format), a `data.table` object, or a list of
#'   `data.table` objects.
#' @param high_mappability_region A `data.table` object of mappability scores, or a
#'   character vector of the mappability score file path.
#' @param mappability_threshold The mappability score threshold, below which the
#'   region is considered as low-mappability and excluded from analysis.
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
           gc_correct = TRUE,
           blacklist_region = NULL,
           high_mappability_region = NULL,
           min_mapq = 30L,
           min_fraglen = 50L,
           max_fraglen = 1000L,
           exclude_soft_clipping = FALSE) {
    assertthat::are_equal(window_size %% step_size, 0)
    assertthat::assert_that(!is.null(high_mappability_region))

    if (is.null(chrom)) {
      chrom <- as.character(unique(seqnames(fragment_data)))
    } else {
      chrom <- as.character(chrom)
    }

    if (length(chrom) > 1) {
      # For multiple chromosomes, sequentially
      ifs <-
        chrom %>% map(function(chrom) {
          ifs_score(
            fragment_data,
            chrom = chrom,
            window_size = window_size,
            step_size = step_size,
            gc_correct = gc_correct,
            blacklist_region = blacklist_region,
            high_mappability_region = high_mappability_region,
            min_mapq = min_mapq,
            min_fraglen = min_fraglen,
            max_fraglen = max_fraglen,
            exclude_soft_clipping = exclude_soft_clipping
          )
        }) %>%
        do.call(c, args = .)
      return(ifs)
    }

    fragment_data <-
      preprocess_fragment(
        fragment_data,
        chrom = chrom,
        min_mapq = min_mapq,
        min_fraglen = min_fraglen,
        max_fraglen = max_fraglen,
        exclude_soft_clipping = exclude_soft_clipping,
        blacklist_region = blacklist_region
      )

    log_mem("Done preprocessing fragment")

    ifs <- calculate_raw_ifs(fragment_data, window_size = window_size, step_size = step_size)
    ifs <- bedtorch::as.GenomicRanges(ifs)

    log_mem("Done calculating raw IFS")

    # Exclude low-mappability regions
    logging::loginfo("Excluding low mappability regions ...")
    if (!is_null(high_mappability_region)) {
      if (is.character(high_mappability_region)) {
        high_mappability_region <- bedtorch::read_bed(high_mappability_region, range = chrom)
      } else {
        high_mappability_region <-
          keepSeqlevels(high_mappability_region, chrom, pruning.mode = "coarse")
      }
      log_mem("Done loading high mappability regions")

      # Intervals in ifs and high_mappability_region should be exactly matched (200bp window)
      ifs <- ifs[queryHits(findOverlaps(ifs, high_mappability_region, type = "equal"))]
      rm(high_mappability_region)
    }

    log_mem("Done mappability filtering")

    # Don't perform GC correctoion
    if (gc_correct) {
      logging::loginfo("Performing GC correction ...")
      ifs <- gc_correct(ifs, span = 0.75)
    }

    log_mem("Done GC correction")

    # Calculate z-scores
    logging::loginfo("Calculating z-scores ...")
    calc_ifs_z_score(ifs)
  }


#' Calculate z-scores for the IFS
calc_ifs_z_score <- function(ifs) {
  assertthat::are_equal(length(unique(seqnames(ifs))), 1)

  mad_factor <- 1.28

  score_median <- median(ifs$score)
  score_mad <- mad(ifs$score, constant = mad_factor)
  ifs$z_score <- (ifs$score - score_median) / score_mad

  ifs
}


calc_gc <- function(ifs) {
  gcs <- as.numeric(biovizBase::GCcontent(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, unstrand(frags)))
  frags$gc <- gcs
}


# Perform GC-correction for IFS scores
#
# @param ifs A `data.table` of original IFS scores.
#   in-place.
# @param gc A `data.table` of GC contents.
# @param span `span` used in LOESS fit.
# @return A `data.table` containing GC-corrected IFS scores/
# @export
gc_correct <- function(ifs, span = 0.75) {
  assertthat::are_equal(length(unique(seqnames(ifs))), 1)
  # chrom <- as.character(unique(ifs$chrom))
  # assertthat::are_equal(length(chrom), 1)
  # assertthat::are_equal(chrom, as.character(unique(gc$chrom)))

  genome <- GenomeInfoDb::genome(ifs) %>% unique()
  # Only works for hs37-1kg
  stopifnot(genome == "hs37-1kg")
  bsgenome <- BSgenome.Hsapiens.1000genomes.hs37d5::hs37d5

  # Dirty hack: need to adjust the seqinfo
  logging::loginfo("Calculating GC content for each fragment...")
  seqinfo_hs37_1kg <- seqinfo(ifs)

  levels_hs37_1kg <- seqlevels(ifs)
  levels_bsgenome <- seqlevels(bsgenome)

  seqinfo(ifs, new2old = match(levels_bsgenome, levels_hs37_1kg)) <- seqinfo(bsgenome)
  ifs$gc <- biovizBase::GCcontent(bsgenome, ifs) %>% as.numeric()
  # Restore the seqinfo
  seqinfo(ifs, new2old = match(levels_hs37_1kg, levels_bsgenome)) <- seqinfo_hs37_1kg

  ifs$score0 <- ifs$score

  # Use all points
  logging::loginfo(str_interp("Training using LOESS mode, total n = ${length(ifs)}"))
  train_data <- data.table::data.table(gc = ifs$gc, score0 = ifs$score0) %>% slice_sample(n = 2e6L)
  model <-
    loess(
      formula = score0 ~ gc,
      data = train_data,
      span = span,
      control = loess.control(
        surface = "interpolate",
        statistics = "approximate",
        trace.hat = "approximate"
      )
    )
  rm(train_data)

  logging::loginfo("Finished training")
  log_mem("Done GC training")

  logging::loginfo("Applying the model for GC correction ...")
  ifs$score <- ifs$score0 - predict(model, newdata = ifs$gc) + mean(ifs$score0)
  ifs$score <- ifelse(ifs$score < 0, 0, ifs$score)

  # In rare cases, GC-corrected scores can be NA. Don't know why.
  # Here we drop out any of those
  na_score_idx <- is.na(ifs$score)
  if (any(na_score_idx)) {
    logging::logwarn(str_interp("Excluded ${sum(na_score_idx)} records with NA scores"))
    ifs <- ifs[!na_score_idx]
  }
  logging::loginfo("Finished performing GC correction")
  ifs

  # #
  # #
  # # # Implement this using foverlaps
  # # ifs[, start := start + 1]
  # # gc[, start := start + 1]
  # # data.table::setkey(ifs, "chrom", "start", "end")
  # # data.table::setkey(gc, "chrom", "start", "end")
  # # overlap_idx <-
  # #   data.table::foverlaps(gc,
  # #                         ifs,
  # #                         type = "any",
  # #                         which = TRUE,
  # #                         nomatch = NULL)
  # # # x: gc, y: 200bp window
  # # overlap_idx[, gc := gc[xid]$score]
  # # overlap_idx[, gc := mean(gc), by = yid]
  # #
  # # ifs[, start := start - 1]
  # # gc[, start := start - 1]
  # # data.table::setkey(ifs, "chrom", "start", "end")
  # # data.table::setkey(gc, "chrom", "start", "end")
  # # ifs <- ifs[overlap_idx$yid, gc := overlap_idx$gc]
  # #
  # # data.table::setnames(ifs, "score", "score0")
  #
  # # Use all points
  # logging::loginfo(str_interp("Training using LOESS mode, total n = ${nrow(ifs)}"))
  # train_data <- ifs[, .(gc, score0)] %>% slice_sample(n = 2e6L)
  # model <-
  #   loess(
  #     formula = score0 ~ gc,
  #     data = train_data,
  #     control = loess.control(
  #       surface = "interpolate",
  #       statistics = "approximate",
  #       trace.hat = "approximate"
  #     )
  #   )
  # rm(train_data)
  # logging::loginfo("Finished training")
  # logging::loginfo("Applying the model for GC correction ...")
  # ifs[, score := score0 - predict(model, newdata = gc) + mean(score0)]
  # # ifs[, score := model$residuals + mean(score0)]
  #
  # # train_data <- ifs[score0 > 0][1:20e3L, .(gc, score0)]
  # # model <- loess(formula = score0 ~ gc, data = train_data)
  # # ifs[score0 > 0, score := score0 - predict(model, newdata = gc) + mean(score0)]
  #
  # ifs[score < 0, score := 0]
  #
  #
  # # In rare cases, GC-corrected scores can be NA. Don't know why.
  # # Here we drop out any of those
  # na_score_idx <- ifs[, is.na(score)]
  # if (any(na_score_idx)) {
  #   logging::logwarn(str_interp("Excluded ${sum(na_score_idx)} records with NA scores"))
  #   ifs <- ifs[!na_score_idx]
  # }
  # logging::loginfo("Finished performing GC correction")
  #
  # # # Also perform GC-correction for coverage
  # # ifs[score > 0, cov_corrected := cov - lowess(x = gc, y = cov, f = span)$y + mean(cov)]
  # # ifs[is.na(cov_corrected), cov_corrected := 0]
  #
  # # ifs[, .(chrom, start, end, score, cov, gc, score0)]
  # ifs[]
}


#' Perform a sliding window over raw IFS scores
#'
#' Raw IFS scores are calculated for each `step_size` bin. This function
#' performs a sliding window with the width of `window_size`. Each window
#' contains a whole number of `step_size` bins.
#' @param window_size Width of the sliding window. Must be multiples of `step_size`.
#' @param step_size Incremental steps of the sliding window.
slide_window <- function(ifs, window_size, step_size) {
  assertthat::are_equal(window_size %% step_size, 0)
  assertthat::are_equal(length(unique(ifs$chrom)), 1)

  genome <- attr(ifs, "genome")
  stopifnot(is.character(genome) || length(genome) == 1)

  chrom_sizes <- bedtorch::get_seqinfo(genome)
  chrom_sizes <-
    data.table(
      chrom = factor(seqnames(chrom_sizes), levels = as.character(seqnames(chrom_sizes))),
      start = 0L,
      end = GenomeInfoDb::seqlengths(chrom_sizes)
    ) %>%
    bedtorch::as.bedtorch_table(genome = genome)
  chrom_sizes[, size := end - start]

  # chrom_sizes <-
  #   data.table::fread(chrom_sizes,
  #                     col.names = c("chrom", "size"),
  #                     header = FALSE)
  # chrom_sizes[, chrom := factor(chrom)]

  n <- window_size %/% step_size

  # Construct a continuous IFS score track, with filled 0s.

  scaffold <-
    data.table::data.table(chrom = ifs$chrom[1])[chrom_sizes, nomatch = 0, on = "chrom"]
  scaffold <- scaffold[, {
    gid <- 0:(size %/% step_size)
    start <- as.integer(gid * step_size)
    end <- as.integer((gid + 1) * step_size)
    list(chrom, start, end)
  }]

  data.table::setkey(scaffold, "chrom", "start", "end")
  ifs <- ifs[scaffold][is.na(score), score := 0][is.na(cov), cov := 0]
  rm(scaffold)

  ifs <- ifs[, `:=`(
    end = start + window_size,
    score = rollsum(
      score,
      k = n,
      na_pad = TRUE,
      align = "left"
    ),
    cov = rollsum(
      cov,
      k = n,
      na_pad = TRUE,
      align = "left"
    )
  )][!is.na(score) & !is.na(cov)]

  # Eliminate float-point errors
  ifs[abs(score) < 1e-9, score := 0]

  data.table::setkey(ifs, "chrom", "start", "end")
  ifs <- bedtorch::as.bedtorch_table(ifs, genome = genome)
}


# Construct the overall excluded region
#
# @param chrom A character vector. If not `NULL`, only regions in these
#   chromosomes will be considered.
# @param blacklist_region Can be a `data.table` object, a list of `data.table`
#   object, or a character vector of the data file paths.
build_exclude_region <-
  function(chrom = NULL,
           blacklist_region = NULL) {
    stopifnot(length(chrom) == 1)
    if (is.null(blacklist_region)) return(NULL)

    if (is.list(blacklist_region)) {
      blacklist_region <- do.call(c, args = blacklist_region)
    } else if (is_character(blacklist_region)) {
      blacklist_region %<>%
        map(function(region) {
          if (file.exists(region))
            bedtorch::read_bed(region, range = chrom)
          else if (region == "encode.blacklist.hs37-1kg")
            bedtorch::read_bed(system.file("extdata", "wgEncodeDacMapabilityConsensusExcludable.hs37-1kg.bed", package = "cragr"))
          else
            stop(paste0("Invalid region: ", region))
        }) %>%
        do.call(c, args = .)
    }

    bedtorch::merge_bed(blacklist_region)
  }


.build_pcpois <- function() {
  lut <- NULL
  lut_factor <- NULL
  res <- 3
  res_lambda <- 2

  # Calculate the normalization factors
  calc_factor <- function(lambda) {
    results <- rep(1, length(lambda))
    idx <- lambda <= 10
    results[idx] <- lambda[idx] %>%
      map_dbl(function(l) {
        v <- integrate(function(x) exp(x * log(l) - l - lgamma(x + 1)), lower = 0, upper = Inf)$value
        1 / v
      })
    results
  }

  function(x, lambda) {
    if (length(x) != 1 && length(lambda) != 1)
      assertthat::are_equal(length(x), length(lambda))

    x <- round(x, digits = res)
    lambda <- round(lambda, digits = res_lambda)

    # Update LUT factor and calculate factors as needed
    if (is_null(lut_factor)) {
      lut_factor <<-
        data.table::as.data.table(list(lambda = unique(lambda), factor = calc_factor(unique(lambda))))
      data.table::setkey(lut_factor, "lambda")
    } else{
      lut_factor <<- merge(lut_factor,
                           data.table::as.data.table(list(lambda = unique(lambda))),
                           all = TRUE)
      # For lambda values that don't have factors associated yet
      lut_factor[is.na(factor), factor := calc_factor(lambda)]
    }

    # Upate LUT
    if (is_null(lut)) {
      lut <<- data.table::as.data.table(list(lambda = 0.1, x = 0.1, p = 0.1))[lambda == 0]
      data.table::setkey(lut, "lambda", "x")
    }

    dt <- data.table::as.data.table(list(lambda = lambda, x = x))
    data.table::setkey(dt, "lambda", "x")

    dt <- dt[lut_factor, nomatch=0]
    # All factors have been calculated
    assertthat::assert_that(!any(is.na(dt$factor)))
    dt <- lut[dt]

    dcpois <- function(x, lambda, factor) {
      factor * exp(x * log(lambda) - lambda - lgamma(x + 1))
    }

    # Calculate probabilities as needed
    dt <- unique(dt)[
      is.na(p),
      p := {
        # Calculate the CDF using numeric integral
        cdf_cpois <- function(v) {
          if (is.na(v))
            return(NA)
          else if (v <= 0)
            return(0)
          else {
            integrate(function(z)
              dcpois(z, lambda[1], factor[1]),
              lower = 0,
              upper = v)$value
          }
        }
        sapply(x, FUN = cdf_cpois)
      },
      by = lambda][, .(lambda, x, p)]

    # data.table::setkey(dt, "lambda", "x")
    full_lut <- merge(dt, lut, all = TRUE)[, p := ifelse(is.na(p.y), p.x, p.y)][, .(lambda, x, p)]

    results <- merge(data.table::as.data.table(list(lambda = lambda, x = x)),
          full_lut,
          by = c("lambda", "x"),
          all.x = TRUE,
          sort = FALSE)$p

    lut <<- full_lut[!is.na(x) & x > 0]
    results
  }
}

pcpois <- .build_pcpois()


#' Calculate p-values based on Poisson model for each window
#'
#' @param cpois A logical value indicating whether to calculate p-values based
#'   on continuous Poisson model.
#' @export
calc_pois_pval <- function(ifs, cpois = FALSE, threshold = 1e-5) {
  # browser()
  # Standard Poisson model
  grl <- GenomicRanges::split(ifs, seqnames(ifs))
  log_mem("Done spliting IFS data frame")
  grl %>%
    lapply(function(gr) {
      gr$pval <- ppois(gr$score, lambda = mean(gr$score))
      gr$pval_adjust <- p.adjust(gr$pval, method = "BH")
      gr
    }) %>%
    GenomicRanges::GRangesList() %>% unlist()

#
#   # # Alternatively, use MAD-outlier detection
#   # ifs[, pval := {
#   #   score_median <- median_score
#   #   mad <- median(abs(score - score_median))
#   #   score_filtered <- score %>% discard(function(x) abs(x - score_median) > 3 * mad)
#   #   lambda <- mean(score_filtered)
#   #   ppois(score, lambda = lambda)
#   # }]
#
#   # Only those pval <= threshld are used in p-value adjustment
#   ifs[, pval_adjust := {
#     # v <- ifelse(pval > threshold, 1, pval)
#     p.adjust(pval, method = "BH")
#   },
#   by = chrom]
#
#   if (cpois) {
#     ifs[, pval_cpois := pcpois(score, lambda = mean(score)), by = chrom]
#
#     # Only those pval <= threshld are used in p-value adjustment
#     ifs[, pval_cpois_adjust := {
#       # v <- ifelse(pval_cpois > threshold, 1, pval_cpois)
#       p.adjust(pval, method = "BH")
#     },
#     by = chrom]
#   }
#   ifs
}


#' Calculate local p-values based on Poisson model for each window
#'
#' @param local_layout Example: list(`5k` = 5000L, `10k` = 10000L)
#' @export
calc_pois_pval_local <- function(ifs, window_size, step_size, local_layout, cpois = FALSE) {
  assertthat::are_equal(window_size %% step_size, 0)
  assertthat::assert_that(all(unlist(local_layout) %% step_size == 0))

  # Use data.table for the calculation
  seqinfo_original <- seqinfo(ifs)
  ifs <- bedtorch::as.bedtorch_table(ifs)
  # Construct a continuous IFS score track, with filled 0s.
  scaffold <-
    ifs[, .(start = seq(min(start), max(start), by = step_size)), by = chrom][, end := start + window_size]
  data.table::setkey(scaffold, "chrom", "start", "end")
  ifs_expanded <- ifs[scaffold]
  rm(scaffold)
  log_mem("Done expanding IFS data frame")

  # for (local_suffix in names(local_layout)) {
  #   local_width <- local_layout[[local_suffix]]
  #
  # For each row, calculate the rollmean over the "local" region, i.e. 5k or 10k
  # bp around the center
  local_peaks <-
    ifs_expanded[,
        {
          valid_score_idx <- !is.na(score)

          results <- local_layout %>% map(function(local_width) {
            # local_width: may be 5000L or 10000L
            # Local means over the 5k/10k region
            outer_sum <- bedtorch::rollsum(
              score,
              k = local_width %/% step_size,
              na_pad = TRUE,
              align = "center",
              na.rm = TRUE
            )
            inner_sum <- bedtorch::rollsum(
              score,
              k = window_size %/% step_size * 2,
              na_pad = TRUE,
              align = "center",
              na.rm = TRUE
            )

            # How many valid 200-bp windows in the 5k/10 rolling region?
            outer_count <- bedtorch::rollsum(
              as.integer(!is.na(score)),
              k = local_width %/% step_size,
              na_pad = TRUE,
              align = "center"
            )
            inner_count <- bedtorch::rollsum(
              as.integer(!is.na(score)),
              k = window_size %/% step_size * 2,
              na_pad = TRUE,
              align = "center"
            )
            valid_count <- outer_count - inner_count
            score_rollmean <- ifelse(!is.na(valid_count) & valid_count > 0, (outer_sum - inner_sum) / valid_count, NA)
            score_rollmean[score_rollmean < 0] <- 0


            # After calculating rollmean and valid_count, we can safely remove
            # rows with NA scores, which are dummy rows only for rolling
            # purposes
            score_rollmean <- score_rollmean[valid_score_idx]
            valid_count <- valid_count[valid_score_idx]
            score <- score[valid_score_idx]

            # Only consider regions with sufficient ...
            valid_roll_idx <- !is.na(valid_count) & (valid_count >= 30) & !is.na(score_rollmean)

            results <- list()
            # By default, the p-value is 1
            results$pval <-
              rep(1, length(score_rollmean))
            results$pval[valid_roll_idx] <-
              ppois(score[valid_roll_idx], score_rollmean[valid_roll_idx])
            results$pval[is.na(valid_count) | is.na(score_rollmean)] <- NA

            if (cpois) {
              results$pval_cpois <-
                rep(1, length(score_rollmean))
              results$pval_cpois[valid_roll_idx] <-
                pcpois(score[valid_roll_idx], score_rollmean[valid_roll_idx])
              results$pval_cpois[is.na(valid_count) | is.na(score_rollmean)] <- NA
            }

            results
          })

          # results: a list. Each element: 5kb and 10k, in which: pval_* and pval_cpois_*
          # Compare 5k and 10k local probabilities and find the maximum
          pval_local <- local({
            # Column: `5k` or `10k`
            mat <- results %>% map(~ .$pval) %>% unlist() %>%
              matrix(ncol = length(results))
            apply(mat, MARGIN = 1, FUN = function(v) ifelse(all(is.na(v)), NA, max(v, na.rm = TRUE)))
          })

          pval_results <- list(
            start = start[valid_score_idx],
            end = end[valid_score_idx],
            pval_local = pval_local
          )

          if (cpois) {
            pval_cpois_local <- local({
              # Column: `5k` or `10k`
              mat <- results %>% map( ~ .$pval_cpois) %>% unlist() %>%
                matrix(ncol = length(results))
              apply(
                mat,
                MARGIN = 1,
                FUN = function(v)
                  ifelse(all(is.na(v)), NA, max(v, na.rm = TRUE))
              )
            })
            pval_results$pval_cpois_local <- pval_cpois_local
          }

          pval_results
        },
        by = chrom]
  data.table::setkey(local_peaks, "chrom", "start", "end")

  ifs <- ifs[local_peaks]
  ifs <- bedtorch::as.GenomicRanges(ifs)

  levels_original <- seqlevels(seqinfo_original)
  levels_now <- seqlevels(ifs)
  seqinfo(ifs, new2old = match(levels_original, levels_now)) <- seqinfo_original

  ifs
}


#' Call hotspots based on pvalues (and merge)
#' @export
call_hotspot <- function(ifs, use_cpois = FALSE, fdr_cutoff = 0.01, pval_cutoff = 1e-5, local_pval_cutoff = 1e-5, merge_distance = 200) {
  browser()
  ifs <- bedtorch::as.bedtorch_table(ifs)
  ifs <-
    ifs[, .(chrom,
            start,
            end,
            z_score,
            score,
            pval,
            pval_adjust,
            pval_local,
            cov,
            score0)]
  fields <- colnames(ifs)

  if (use_cpois) {
    hotspot <-
      ifs[pval <= pval_cutoff &
            pval_cpois_adjust <= fdr_cutoff &
            pval_cpois_local <= local_pval_cutoff]
  } else {
    hotspot <-
      ifs[pval <= pval_cutoff &
            pval_adjust <= fdr_cutoff & pval_local <= local_pval_cutoff]
  }
  if (nrow(hotspot) == 0) {
    # hotspot is empty
    return(NULL)
  }

  bedtorch::merge_bed(
    bedtorch::as.GenomicRanges(hotspot),
    operation = list(
      z_score = list(
        on = "z_score",
        func = function(x)
          mean(x, na.rm = TRUE)
      ),
      score = list(
        on = "score",
        func = function(x)
          mean(x, na.rm = TRUE)
      ),
      pval = list(
        on = "pval",
        func = function(x)
          max(x, na.rm = TRUE)
      ),
      pval_adjust = list(
        on = "pval_adjust",
        func = function(x)
          max(x, na.rm = TRUE)
      ),
      pval_local = list(
        on = "pval_local",
        func = function(x)
          max(x, na.rm = TRUE)
      ),
      cov = list(
        on = "cov",
        func = function(x)
          mean(x, na.rm = TRUE)
      ),
      score0 = list(
        on = "score0",
        func = function(x)
          mean(x, na.rm = TRUE)
      )
    )
  )
}




# Prepross fragment data
# 1. MAPQ filter
# 2. Fragment length filter
# 3. Blacklist region filter
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
    if (min_mapq > 0 && "mapq" %in% colnames(mcols(fragment_data))) {
      fragment_data <- fragment_data[fragment_data$mapq >= min_mapq]
    }
    logging::loginfo(str_interp("${length(fragment_data)} fragments have passed filter"))

    logging::loginfo("Applying fragment length filter ...")
    min_fraglen <- min_fraglen %||% 0
    max_fraglen <- max_fraglen %||% Inf
    fragment_data <-
      fragment_data[between(width(fragment_data), min_fraglen, max_fraglen)]
    logging::loginfo(str_interp("${length(fragment_data)} fragments have passed filter"))


    if (exclude_soft_clipping &&
        all(c("cigar1", "cigar2") %in% colnames(mcols(fragment_data)))) {
      logging::loginfo("Removing soft-clipping fragments ...")
      fragment_data <-
        fragment_data[!grepl(pattern = "^[0-9]+S", x = fragment_data$cigar1) &
                        !grepl(pattern = "[0-9]+S$", x = fragment_data$cigar2)]
      logging::loginfo(str_interp("${length(fragment_data)} fragments have passed filter"))
    }

    # Exclude fragments from certain regions
    if (!is.null(blacklist_region)) {
      excluded_regions <-
        build_exclude_region(chrom = chrom,
                             blacklist_region = blacklist_region)
    } else
      excluded_regions <- NULL

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
    midpoint <-
      round((start(fragment_data) + end(fragment_data)) / 2)
    start(fragment_data) <- midpoint
    width(fragment_data) <- 1
    fragment_data <- sort(fragment_data)

    # Exclude dark regions
    logging::loginfo("Applying blacklist region filter ...")
    if (!is.null(excluded_regions)) {
      logging::loginfo("Excluding blacklist regions ...")
      hits <- queryHits(findOverlaps(fragment_data, excluded_regions))
      # Caution: hits may be empty!
      if (length(hits) > 0) {
        fragment_data <- fragment_data[-hits]
      }
    }
    logging::loginfo(str_interp("${length(fragment_data)} fragments have passed filter"))

    fragment_data
  }


# Calculate raw IFS scores, i.e. the ones before GC-correction
calculate_raw_ifs <-
  function(fragment_data, window_size, step_size) {
    # Get coverage and calcuate IFS score for each 20bp (step-size) window
    logging::loginfo("Calculating raw IFS scores ...")

    # Here we need do grouping and aggregating.
    # It seems data.table is way faster
    fragment_data$window_id <-
      (start(fragment_data) - 1) %/% step_size

    fragment_data <- bedtorch::as.bedtorch_table(fragment_data)
    avg_len <- mean(fragment_data$length)
    ifs <- fragment_data[,
      .(
        chrom = chrom[1],
        start = as.integer(window_id * step_size),
        end = as.integer((window_id + 1) * step_size),
        # score = .N + sum(length, na.rm = TRUE) / avg_len,
        sum_len = sum(length, na.rm = TRUE),
        cov = as.integer(.N)
      ),
      by = window_id
    ][, window_id := NULL]
    ifs <-
      ifs[, .(chrom,
              start,
              end,
              score = cov + sum_len / avg_len,
              cov,
              frag_len = sum_len / cov)]

    ifs %<>% bedtorch::as.bedtorch_table(genome = attr(fragment_data, "genome"))
    data.table::setkeyv(ifs, c("chrom", "start", "end"))

    logging::logdebug("Applying rollmean")
    # Perform rolling sum over the sliding windows, therefore we have results
    # for rolling windows (200bp, window_size) at step size of (20bp, step_size)
    ifs <-
      slide_window(ifs,
                   window_size = window_size,
                   step_size = step_size,
                   avg_len = avg_len
      )

    # ifs[, avg_len := avg_len]
    return(ifs)
  }


#' Calculate raw IFS scores
#'
#' Calculate raw IFS scores from coverage and fragment length, and do the
#' mappability filtering.
#'
#' @param fragment_data A `GRanges` data frame containing cfDNA fragment data.
#' @param chrom Calculate IFS scores for certain chromosomes. If `NULL`, all
#'   chromosomes in `fragment_data` will be used.
#' @param window_size Width of the sliding window. Must be multiples of
#'   `step_size`.
#' @param step_size Incremental steps of the sliding window.
#' @param gc_correct Logical value. Whether to perform GC correction. If `TRUE`,
#'   the corresponding `BSgenome` should be present. Currently, only support
#'   `BSgenome.Hsapiens.1000genomes.hs37d5`.
#' @param blacklist_region Those fragments whose midpoint fall within any of the
#'   excluded regions will not be used in the analysis. `blacklist_region` can
#'   be either a character vector that contains the names of the files defining
#'   the regions (should be in BED format), a `GRanges` object, or a list of
#'   `GRanges` objects.
#' @param high_mappability_region Regions with high mappability scores. Only the
#'   windows within such regions with be used in the analysis. Must be a
#'   `GRanges`, or path of the BED file path.
#' @param min_mapq Exclude fragments whose `MAPQ` scores are smaller than
#'   `min_mapq`.
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
    stopifnot(is(fragment_data, "GRanges"))
    assertthat::are_equal(window_size %% step_size, 0)
    # assertthat::assert_that(!is.null(high_mappability_region))

    if (is.null(chrom)) {
      chrom <- as.character(unique(seqnames(fragment_data)))
    } else {
      chrom <- as.character(chrom)
    }

    assertthat::assert_that(length(chrom) == 1)
    # if (length(chrom) > 1) {
    #   # For multiple chromosomes, sequentially
    #   ifs <-
    #     chrom %>%
    #     map(function(chrom) {
    #       ifs_score(
    #         fragment_data,
    #         chrom = chrom,
    #         window_size = window_size,
    #         step_size = step_size,
    #         gc_correct = gc_correct,
    #         blacklist_region = blacklist_region,
    #         high_mappability_region = high_mappability_region,
    #         min_mapq = min_mapq,
    #         min_fraglen = min_fraglen,
    #         max_fraglen = max_fraglen,
    #         exclude_soft_clipping = exclude_soft_clipping
    #       )
    #     }) %>%
    #     do.call(c, args = .)
    #   return(ifs)
    # }

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

    avg_len <- mean(fragment_data$length)
    frag_cnt <- length(fragment_data$length)
    ifs <-
      calculate_raw_ifs(fragment_data,
                        window_size = window_size,
                        step_size = step_size
      )
    ifs <- bedtorch::as.GenomicRanges(ifs)

    log_mem("Done calculating raw IFS")

    # Exclude low-mappability regions
    if (!is_null(high_mappability_region)) {
      logging::loginfo("Excluding low mappability regions ...")
      if (is.character(high_mappability_region)) {
        high_mappability_region <-
          bedtorch::read_bed(high_mappability_region, range = chrom)
      } else {
        high_mappability_region <-
          keepSeqlevels(high_mappability_region, chrom, pruning.mode = "coarse")
      }
      log_mem("Done loading high mappability regions")

      # Intervals in ifs and high_mappability_region should be exactly matched (200bp window)
      ifs <-
        ifs[queryHits(findOverlaps(ifs, high_mappability_region, type = "equal"))]
      rm(high_mappability_region)
      log_mem("Done mappability filtering")
    }

    return(list(ifs = ifs, avg_len = avg_len, frag_cnt = frag_cnt))
  }


# postprocess_ifs <- function(ifs, gc_correct = TRUE) {
#   # Don't perform GC correctoion
#   if (gc_correct) {
#     logging::loginfo("Performing GC correction ...")
#     ifs <- gc_correct(ifs, span = 0.75)
#   }
#
#   log_mem("Done GC correction")
#
#   # Calculate z-scores
#   logging::loginfo("Calculating z-scores ...")
#   ifs <- calc_ifs_z_score(ifs)
#   ifs
# }


# Calculate z-scores for the IFS
#' @export
calc_ifs_z_score <- function(ifs) {
  assertthat::are_equal(length(unique(seqnames(ifs))), 1)

  score <- exclude_outlier(ifs$score)
  score_mean <- mean(score)
  score_sd <- sd(score)
  ifs$z_score <- (ifs$score - score_mean) / score_sd

  ifs
}


#' Calculate GC content for each fragment
#'
#' @param ifs A `GRanges` object. In order to calculate GC content, `ifs` must
#'   have a valid genome. Currently, only hs37-1kg is supported. Furthermore, you need
#'
#' @export
calc_gc <- function(ifs) {
  genome_name <- GenomeInfoDb::genome(ifs) %>% unique()
  assertthat::assert_that(is_scalar_character(genome_name) &&
                (genome_name %in% c("GRCh37", "hs37-1kg", "GRCh38")))

  bsgenome <- switch(
    genome_name,
    "GRCh37" = BSgenome::getBSgenome(genome = "BSgenome.Hsapiens.1000genomes.hs37d5", load.only = TRUE),
    "hs37-1kg" = BSgenome::getBSgenome(genome = "BSgenome.Hsapiens.1000genomes.hs37d5", load.only = TRUE),
    "GRCh38" = BSgenome::getBSgenome(genome = "BSgenome.Hsapiens.NCBI.GRCh38", load.only = TRUE),
    stop(paste0("Invalid genome: ", genome_name))
  )

  seqlevels(ifs) <- seqlevels(bsgenome)
  seqinfo(ifs) <- seqinfo(bsgenome)

  # ifs may contain intervals out of the chromosome, make sure to exclude them
  ifs <- GenomicRanges::trim(ifs)
  common_width <- median(GenomicRanges::width(ifs))
  ifs <- ifs[GenomicRanges::width(ifs) == common_width]

  genome_seq <- BSgenome::getSeq(bsgenome, ifs)
  gc <- Biostrings::letterFrequency(genome_seq, letters = "CG")[, 1]
  fourbases <- Biostrings::letterFrequency(genome_seq, letters = "ATCG")[, 1]
  gc <- gc / fourbases

  # Discard a bin if it contains too many Ns
  mask_threshold <- 0.05
  base_masked <- Biostrings::letterFrequency(genome_seq, letters = "N", as.prob = TRUE)[, 1] > mask_threshold
  gc[base_masked] <- NA

  ifs$gc <- gc
  return(ifs)
}

#' Perform GC-correction for IFS scores
#'
#' @param span `span` used in LOESS fit.
#' @param max_training_dataset LOESS can be slow for a large dataset. If it
#'   contains more data points than this parameter, we will sample a subset and
#'   do the training.
#' @param method If `standard`, the built-in `loess` function will be used for the calculation. If `caret`, will use `gamLOESS`
#' @export
gc_correct <- function(ifs, span = 0.75, max_training_dataset = 1e6L, method = "standard", ...) {
  if (method == "standard") {
    return(gc_correct_standard(ifs, span, max_training_dataset, ...))
  } else if (method == "caret") {
    return(gc_correct_caret(ifs, span, max_training_dataset, ...))
  } else {
    stop(paste0("Unsupported method: ", method))
  }
}


gc_correct_caret <-
  function(ifs,
           span = 0.75,
           max_training_dataset = 1e6L,
           thread = 1L,
           ...) {
    assertthat::are_equal(length(unique(seqnames(ifs))), 1)
    assertthat::assert_that("gc" %in% colnames(mcols(ifs)))
    
    logging::loginfo("Performin GC correction through caret/gamLoess ...")
    
    ifs$score0 <- ifs$score
    
    sel_idx <- 1:length(ifs)
    # Exclude GC/score with NAs, otherwise the model training may fail
    sel_idx <- sel_idx[!is.na(ifs$gc) & !is.na(ifs$score0)]
    
    # Use all points
    if (length(sel_idx) > max_training_dataset) {
      logging::loginfo(
        str_interp(
          "Training using LOESS mode, total n = ${length(sel_idx)}, subsample n = ${max_training_dataset}"
        )
      )
      sel_idx <- sample(sel_idx, size = max_training_dataset)
    } else {
      logging::loginfo(str_interp("Training using LOESS mode, total n = ${length(sel_idx)}"))
    }
    
    library(caret)
    library(doParallel)
    
    # LOESS model ----
    logging::loginfo("Building LOESS model ...")
    
    cl <- makeForkCluster(max(1, thread))
    registerDoParallel(cl)
    
    model <- train(
      x = data.frame(gc = ifs$gc[sel_idx]),
      y = ifs$score0[sel_idx],
      method = "gamLoess"
    )
    stopCluster(cl)
    
    logging::loginfo("Finished training")
    gc()
    log_mem("Done GC training")
    
    logging::loginfo("Applying the model for GC correction ...")
    
    # If GC is NA, model prediction may fail. Let's exclude them from the analysis
    na_idx <- is.na(ifs$gc)
    pred <- predict(model, newdata = data.frame(gc = ifs$gc[!na_idx]))
    ifs$score <- NA
    ifs$score[!na_idx] <-
      pmax(0, ifs$score0[!na_idx] - pred + mean(ifs$score0, na.rm = TRUE))
    
    # In rare cases, GC-corrected scores can be NA. Don't know why.
    # Here we drop out any of those
    na_score_idx <- is.na(ifs$score)
    if (any(na_score_idx)) {
      logging::logwarn(str_interp("Excluded ${sum(na_score_idx)} records with NA scores"))
      ifs <- ifs[!na_score_idx]
    }
    logging::loginfo("Finished performing GC correction")
    ifs
  }


gc_correct_standard <-
  function(ifs,
           span = 0.75,
           max_training_dataset = 1e6L,
           ...) {
    assertthat::are_equal(length(unique(seqnames(ifs))), 1)
    assertthat::assert_that("gc" %in% colnames(mcols(ifs)))
    
    logging::loginfo("Performin GC correction through standard LOESS regression ...")

    ifs$score0 <- ifs$score
    sel_idx <- 1:length(ifs)

    # Use all points
    if (length(sel_idx) > max_training_dataset) {
      logging::loginfo(
        str_interp(
          "Training using LOESS mode, total n = ${length(sel_idx)}, subsample n = ${max_training_dataset}"
        )
      )
      sel_idx <- sample(sel_idx, size = max_training_dataset)
    } else {
      logging::loginfo(str_interp("Training using LOESS mode, total n = ${length(sel_idx)}"))
    }

    # A subset of IFS for model training
    ifs_train <- ifs[sel_idx]
    
    # Exclude outliers from the training dataset
    outlier_flag <- exclude_outlier(ifs_train$score0, mark = TRUE)
    train_data <-
      data.table::data.table(gc = ifs_train$gc[!outlier_flag], score0 = ifs_train$score0[!outlier_flag]) %>%
      slice_sample(n = max_training_dataset)
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
    ifs$score <-
      ifs$score0 - predict(model, newdata = ifs$gc) + mean(ifs$score0)
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
  }


# Perform a sliding window over raw IFS scores
#
# Raw IFS scores are calculated for each `step_size` bin. This function
# performs a sliding window with the width of `window_size`. Each window
# contains a whole number of `step_size` bins.
# The result is a collection of overlapping `window_size` bins.
# @param window_size Width of the sliding window. Must be multiples of `step_size`.
# @param step_size Incremental steps of the sliding window.
slide_window <- function(ifs, window_size, step_size, avg_len) {
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
  ifs <-
    ifs[scaffold][is.na(score), score := 0][is.na(cov), cov := 0]
  rm(scaffold)

  ifs <- ifs[, `:=`(
    end = as.integer(start + window_size),
    score = rollsum(
      score,
      k = n,
      na_pad = TRUE,
      align = "left"
    ),
    cov = as.integer(rollsum(
      cov,
      k = n,
      na_pad = TRUE,
      align = "left"
    ))
  )][!is.na(score) & !is.na(cov)]
  ifs <- ifs[cov > 0]
  ifs[, frag_len := (score - cov) / cov * avg_len]

  # Eliminate float-point errors
  ifs[abs(score) < 1e-9, score := 0]

  data.table::setkey(ifs, "chrom", "start", "end")
  bedtorch::as.bedtorch_table(ifs, genome = genome)
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
    if (is.null(blacklist_region)) {
      return(NULL)
    }

    if (is.list(blacklist_region)) {
      blacklist_region <- do.call(c, args = blacklist_region)
    } else if (is_character(blacklist_region)) {
      blacklist_region %<>%
        map(function(region) {
          if (file.exists(region)) {
            gr <- bedtorch::read_bed(region)
            gr[GenomicRanges::seqnames(gr) == chrom]
          } else if (region == "encode.blacklist.hs37-1kg") {
            bedtorch::read_bed(
              system.file(
                "extdata",
                "wgEncodeDacMapabilityConsensusExcludable.hs37-1kg.bed",
                package = "cragr"
              )
            )
          } else {
            stop(paste0("Invalid region: ", region))
          }
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
        v <-
          integrate(function(x) {
            exp(x * log(l) - l - lgamma(x + 1))
          },
          lower = 0,
          upper = Inf
          )$value
        1 / v
      })
    results
  }

  function(x, lambda) {
    if (length(x) != 1 && length(lambda) != 1) {
      assertthat::are_equal(length(x), length(lambda))
    }

    x <- round(x, digits = res)
    lambda <- round(lambda, digits = res_lambda)

    # Update LUT factor and calculate factors as needed
    if (is_null(lut_factor)) {
      lut_factor <<-
        data.table::as.data.table(list(lambda = unique(lambda), factor = calc_factor(unique(lambda))))
      data.table::setkey(lut_factor, "lambda")
    } else {
      lut_factor <<- merge(lut_factor,
                           data.table::as.data.table(list(lambda = unique(lambda))),
                           all = TRUE
      )
      # For lambda values that don't have factors associated yet
      lut_factor[is.na(factor), factor := calc_factor(lambda)]
    }

    # Upate LUT
    if (is_null(lut)) {
      lut <<-
        data.table::as.data.table(list(
          lambda = 0.1,
          x = 0.1,
          p = 0.1
        ))[lambda == 0]
      data.table::setkey(lut, "lambda", "x")
    }

    dt <- data.table::as.data.table(list(lambda = lambda, x = x))
    data.table::setkey(dt, "lambda", "x")

    dt <- dt[lut_factor, nomatch = 0]
    # All factors have been calculated
    assertthat::assert_that(!any(is.na(dt$factor)))
    dt <- lut[dt]

    dcpois <- function(x, lambda, factor) {
      factor * exp(x * log(lambda) - lambda - lgamma(x + 1))
    }

    # Calculate probabilities as needed
    dt <- unique(dt)[is.na(p),
                     p := {
                       # Calculate the CDF using numeric integral
                       cdf_cpois <- function(v) {
                         if (is.na(v)) {
                           return(NA)
                         } else if (v <= 0) {
                           return(0)
                         } else {
                           integrate(function(z) {
                             dcpois(z, lambda[1], factor[1])
                           },
                           lower = 0,
                           upper = v
                           )$value
                         }
                       }
                       sapply(x, FUN = cdf_cpois)
                     },
                     by = lambda
    ][, .(lambda, x, p)]

    full_lut <-
      merge(dt, lut, all = TRUE)[, p := ifelse(is.na(p.y), p.x, p.y)][, .(lambda, x, p)]

    results <-
      merge(
        data.table::as.data.table(list(lambda = lambda, x = x)),
        full_lut,
        by = c("lambda", "x"),
        all.x = TRUE,
        sort = FALSE
      )$p

    lut <<- full_lut[!is.na(x) & x > 0]
    results
  }
}

pcpois <- .build_pcpois()


# @param mark if TRUE, return a logical vector indicating positions of outliers
exclude_outlier <- function(x, mark = FALSE, threshold = 5) {
  n1 <- sum(!is.na(x))

  if (n1 < 10) {
    if (mark)
      return(rep(FALSE, length(x)))
    else
      return(x)
  }

  mad_x <- mad(x, na.rm = TRUE)
  outlier_flag <- !is.na(x) & (abs(x - median(x, na.rm = TRUE)) >= threshold * mad_x)
  n2 <- sum(!is.na(x) & !outlier_flag)

  if (n2 / n1 < 0.9) {
    logging::logwarn(str_interp("Tentative outlier exclusion: #${n1} -> #${n2}. Pass rate: ${n2/n1}"))
    if (mark)
      return(rep(FALSE, length(x)))
    else
      return(x)
  }

  if (mark)
    return(outlier_flag)
  else
    return(x[!outlier_flag])
}


#' Calculate p-values based on Poisson model for each window
#'
#' @param cpois A logical value indicating whether to calculate p-values based
#'   on continuous Poisson model. Currently, must be `FALSE`.
#' @export
calc_pois_pval <- function(ifs,
                           cpois = FALSE,
                           # model = c("poisson", "negbinom"),
                           threshold = 1e-5) {
  assertthat::assert_that(!cpois)
  # model <- match.arg(model)
  # model <- c("poisson", "negbinom")

  # Calculate p-values using Poisson model
  pval_ppois <- function(x) {
    x_no_outlier <- exclude_outlier(x)
    stats::ppois(x, lambda = mean(x_no_outlier))
  }

  # Calculate p-values using negative binomial model
  pval_pnbinom <- function(x) {
    x_no_outlier <- exclude_outlier(x)
    prob <- mean(x_no_outlier) / var(x_no_outlier)
    size <-
      (mean(x_no_outlier) ** 2) / var(x_no_outlier) / (1 - prob)
    stats::pnbinom(x, size = size, prob = prob)
  }

  # Standard Poisson model
  grl <- GenomicRanges::split(ifs, seqnames(ifs))
  log_mem("Done spliting IFS data frame")
  grl %>%
    lapply(function(gr) {
      if (length(gr) == 0)
        return(gr)

      # Poisson-based global p-values
      gr$pval_pois <- pval_ppois(gr$score)
      # Poisson-based global p-values with FDR correction
      gr$pval_pois_adjust <- p.adjust(gr$pval_pois, method = "BH")
      # NB-based global p-values
      gr$pval_nbinom <- pval_pnbinom(gr$score)
      gr
    }) %>%
    GenomicRanges::GRangesList() %>%
    unlist()
}


#' Calculate local p-values based on Poisson model for each window
#'
#' @param local_layout Example: list(`5k` = 5000L, `10k` = 10000L)
#' @param cpois Currently must be `FALSE`.
#' @export
calc_pois_pval_local <-
  function(ifs,
           window_size,
           step_size,
           local_layout,
           cpois = FALSE) {
    assertthat::are_equal(window_size %% step_size, 0)
    assertthat::assert_that(all(unlist(local_layout) %% step_size == 0))

    # Use data.table for the calculation
    seqinfo_original <- seqinfo(ifs)
    # Must have seqinfo
    assertthat::assert_that(all(!is.na(GenomeInfoDb::genome(seqinfo_original))))

    ifs <- bedtorch::as.bedtorch_table(ifs)
    # Mark outliers for each chromosome
    ifs[, is_outlier := exclude_outlier(score, mark = TRUE), by = chrom]

    # Construct a continuous IFS score track, with filled 0s.
    scaffold <-
      ifs[, .(start = seq(min(start), max(start), by = step_size)), by = chrom][, end := start + window_size]
    data.table::setkey(ifs, "chrom", "start", "end")
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
                     # Because we just expanded the IFS data table, some rows will contain
                     # NA scores
                     valid_score_idx <- !is.na(score)

                     # Used for var calculation
                     score_sq <- score ** 2

                     results <-
                       names(local_layout) %>% map(function(local_suffix) {
                         # local_suffix: for example: _5k
                         # local_width: 5000L, or 10000L, for example
                         local_width <- local_layout[[local_suffix]]

                         # Outer rim: center -       outer_shift -> center -> center      +       outer_shift
                         # Inner rim:       center - inner_shift -> center -> center + inner_shift
                         #
                         # Only use outer_rim - inner_rim to determin the background
                         outer_shift <- as.integer(local_width / step_size / 2)
                         inner_shift <- as.integer(window_size / step_size)

                         # local_width: may be 5000L or 10000L
                         # Local means over the 5k/10k region
                         outer_sum <- bedtorch::rollsum(
                           ifelse(is_outlier, NA, score),
                           k = 2 * outer_shift + 1,
                           na_pad = TRUE,
                           align = "center",
                           na.rm = TRUE
                         )
                         inner_sum <- bedtorch::rollsum(
                           ifelse(is_outlier, NA, score),
                           k = 2 * inner_shift + 1,
                           na_pad = TRUE,
                           align = "center",
                           na.rm = TRUE
                         )
                         outer_sum_sq <- bedtorch::rollsum(
                           ifelse(is_outlier, NA, score_sq),
                           k = 2 * outer_shift + 1,
                           na_pad = TRUE,
                           align = "center",
                           na.rm = TRUE
                         )
                         inner_sum_sq <- bedtorch::rollsum(
                           ifelse(is_outlier, NA, score_sq),
                           k = 2 * inner_shift + 1,
                           na_pad = TRUE,
                           align = "center",
                           na.rm = TRUE
                         )

                         # How many valid 200-bp windows in the 5k/10 rolling region?
                         # Valid: sore is non-NA, is_outlier is FALSE
                         outer_count <- bedtorch::rollsum(
                           as.integer(!is.na(score) & !is.na(is_outlier) & !is_outlier),
                           k = 2 * outer_shift + 1,
                           na_pad = TRUE,
                           align = "center"
                         )
                         inner_count <- bedtorch::rollsum(
                           as.integer(!is.na(score) & !is.na(is_outlier) & !is_outlier),
                           k = 2 * inner_shift + 1,
                           na_pad = TRUE,
                           align = "center"
                         )
                         valid_count <- outer_count - inner_count

                         score_rollmean <-
                           ifelse(!is.na(valid_count) &
                                    valid_count > 1,
                                  (outer_sum - inner_sum) / valid_count,
                                  NA
                           )
                         score_rollmean[score_rollmean < 0] <- 0

                         score_sq_rollmean <-
                           ifelse(!is.na(valid_count) &
                                    valid_count > 1,
                                  (outer_sum_sq - inner_sum_sq) / valid_count, NA)
                         score_sq_rollmean[score_sq_rollmean < 0] <- 0
                         score_rollvar <- valid_count/(valid_count - 1) * (score_sq_rollmean - score_rollmean**2)


                         # After calculating rollmean and valid_count, we can safely remove
                         # rows with NA scores, which are dummy rows only for rolling
                         # purposes
                         score_rollmean <-
                           score_rollmean[valid_score_idx]
                         score_rollvar <-
                           score_rollvar[valid_score_idx]
                         valid_count <- valid_count[valid_score_idx]
                         score <- score[valid_score_idx]

                         # Only consider regions with sufficient ...
                         valid_roll_idx <-
                           !is.na(valid_count) &
                           (valid_count >= 30) &
                           !is.na(score_rollmean) &
                           !is.na(score_rollvar)

                         results <- list()

                         # Local pois p-values
                         v <- rep(NA, length(score_rollmean))
                         v[valid_roll_idx] <- ppois(score[valid_roll_idx],
                                                    lambda = score_rollmean[valid_roll_idx])
                         results[[paste0("pval_pois_", local_suffix)]] <- v

                         # Local NB p-values
                         v <- rep(NA, length(score_rollmean))
                         nbinom_mu <- score_rollmean[valid_roll_idx]
                         nbinom_var <- score_rollvar[valid_roll_idx]
                         nbinom_size <- ifelse(
                           # in these cases we can't infer the nbinom distribution
                           nbinom_var == nbinom_mu | nbinom_var < nbinom_mu,
                           NA,
                           nbinom_mu ** 2 / (nbinom_var - nbinom_mu)
                         )

                         v[valid_roll_idx] <- pnbinom(score[valid_roll_idx],
                                                      mu = nbinom_mu,
                                                      size = nbinom_size)

                         results[[paste0("pval_nbinom_", local_suffix)]] <- v

                         results[[paste0("mu_", local_suffix)]] <- score_rollmean
                         # results[[paste0("var_", local_suffix)]] <- score_rollvar

                         results
                       })

                     # Example:
                     #
                     # List of 8
                     # $ pval_pois_5k   : num [1:1323099] NA NA NA 1 1 1 1 1 1 1 ...
                     # $ pval_nbinom_5k : num [1:1323099] NA NA NA NaN NaN NaN NaN NaN NaN NaN ...
                     # $ mu_5k          : num [1:1323099] NA NA NA 0 0 0 0 0 0 0 ...
                     # $ var_5k         : num [1:1323099] NA NA NA 0 0 0 0 0 0 0 ...
                     # $ pval_pois_10k  : num [1:1323099] NA NA NA 0.897 0.889 ...
                     # $ pval_nbinom_10k: num [1:1323099] NA NA NA 0.905 0.897 ...
                     # $ mu_10k         : num [1:1323099] NA NA NA 0.108 0.118 ...
                     # $ var_10k        : num [1:1323099] NA NA NA 0.127 0.139 ...
                     results <- do.call(c, args = results)
                     c(
                       list(start = start[valid_score_idx],
                            end = end[valid_score_idx]),
                       names(results) %>%
                         map(function(name)
                           results[[name]]) %>%
                         set_names(names(results))
                     )
                   },
                   by = chrom
      ]
    data.table::setkey(local_peaks, "chrom", "start", "end")

    ifs <- ifs[local_peaks]

    ifs[, is_outlier := NULL]

    # Aggregate all local pois p-values by taking the maximum
    # Example: pval_pois_5k, pval_pois_10k, ...
    pval_pois_local_cols <- paste0("pval_pois_", names(local_layout))
    v <- do.call(pmax, args = ifs[, ..pval_pois_local_cols])
    ifs[, pval_pois_local := v]
    ifs[, (pval_pois_local_cols) := NULL]

    # Aggregate all NB p-values by taking the maximum
    # Example: pval_nbinom, pval_nbinom_5k, pval_nbinom_10k, ...
    pval_nb_cols <- c("pval_nbinom", paste0("pval_nbinom_", names(local_layout)))

    # Here we only use local background p-values to calculate the results
    # pval_nb_cols <- paste0("pval_nbinom_", names(local_layout))
    v <- do.call(pmax, args = ifs[, ..pval_nb_cols])
    ifs[, pval := v]
    ifs[, pval_adjust := p.adjust(pval, method = "BH")]
    # ifs[, (pval_nb_cols) := NULL]
    # # Remove global background p-values
    # ifs[, pval_nbinom := NULL]

    ifs <- bedtorch::as.GenomicRanges(ifs)

    levels_original <- seqlevels(seqinfo_original)
    levels_now <- seqlevels(ifs)
    seqinfo(ifs, new2old = match(levels_original, levels_now)) <-
      seqinfo_original

    ifs
  }

#' Call hotspots using FDR only
#' @export
call_hotspot <- function(
  ifs,
  fdr_cutoff = 0.2,
  pval_cutoff = 1e-5,
  local_pval_cutoff = 1e-5,
  method = c("pois", "nb")
) {
  method <- match.arg(method)

  # Relocate column orders
  ifs_md <- mcols(ifs) %>%
    as_tibble() %>%
    relocate(z_score, 1)
  mcols(ifs) <- ifs_md

  assertthat::assert_that(method %in% c("pois", "nb"))
  if (method == "pois") {
    # All local p-values (pval_pois_5k and alike) should be filtered
    pval_local_filter_result <-
      ifs_md %>%
      select(starts_with("pval_pois_")) %>%
      select(-"pval_pois_adjust") %>%
      mutate(id = seq.int(n())) %>%
      mutate(flag = if_all(starts_with("pval_pois_"),
                           ~ !is.na(.) & . <= local_pval_cutoff)) %>%
      .$flag

    other_filter_result <- with(
      ifs_md,
      !is.na(pval_pois) &
        pval_pois <= pval_cutoff &
        !is.na(pval_pois_adjust) &
        pval_pois_adjust <= fdr_cutoff
    )

    hotspot <- ifs[pval_local_filter_result & other_filter_result]
  } else if (method == "nb") {
    hotspot_idx <-
      with(ifs_md,!is.na(pval_adjust) & pval_adjust <= fdr_cutoff)
    hotspot <- ifs[hotspot_idx]
  } else {
    stop()
  }

  if (length(hotspot) == 0) {
    # hotspot is empty
    return(NULL)
  }

  # mcols(hotspot) <- NULL
  hotspot
}

